#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


def is_project_source(source: Path, repo_root: Path) -> bool:
    """Exclude third-party deps and generated build files."""
    try:
        rel = source.relative_to(repo_root)
    except ValueError:
        return False
    rel_text = rel.as_posix()
    if rel_text.startswith("build/") or "_deps/" in rel_text:
        return False
    return True


KNOWN_BENCHMARK_TARGETS = {
    "benchmarks",
    "bench_rmm",
    "bench_rmm_sdsl",
    "louds_tree_benchmarks",
    "alignment_comparison",
}

HEADER_EXTENSIONS = {
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".inc",
    ".ipp",
    ".tcc",
}

BUILD_INFRA_FILES = {
    "CMakeLists.txt",
    "CMakePresets.json",
}

DIFF_HUNK_RE = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")

CPP_FUNCTION_START_RE = re.compile(
    r"^\s*"
    r"(?:template\s*<[^>]*>\s*)?"
    r"(?:(?:inline|constexpr|consteval|constinit|static|friend|virtual|explicit)\s+)*"
    r"[A-Za-z_~][\w:<>,\s\*&\[\]]*\s+"
    r"([~A-Za-z_][A-Za-z0-9_]*)\s*"
    r"\([^;{}]*\)\s*"
    r"(?:const\s*)?"
    r"(?:noexcept\s*)?"
    r"(?:->\s*[^\{]+)?\{"
)


@dataclass
class CompileCommandEntry:
    directory: Path
    source: Path
    arguments: list[str]
    output: Path | None
    target: str | None
    dependencies: set[Path] = field(default_factory=set)
    dep_error: str | None = None


@dataclass
class AstImpactResult:
    benchmark_names: set[str] = field(default_factory=set)
    affected_names: set[str] = field(default_factory=set)
    ast_error: str | None = None


def run_command(
    args: list[str],
    cwd: Path,
    check: bool = True,
    timeout: float | None = 60.0,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=str(cwd),
        text=True,
        capture_output=True,
        check=check,
        timeout=timeout,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze benchmark impact between baseline and HEAD via "
            "compile_commands dependency mapping and clang AST analysis."
        )
    )
    parser.add_argument(
        "--baseline",
        default="main",
        help="Baseline ref used as <baseline>...HEAD (default: main).",
    )
    parser.add_argument(
        "--head",
        default="HEAD",
        help="Contender ref (default: HEAD).",
    )
    parser.add_argument(
        "--compile-commands",
        default=None,
        help=(
            "Path to compile_commands.json. If omitted, auto-discovers most "
            "recent build/**/compile_commands.json."
        ),
    )
    parser.add_argument(
        "--clangxx",
        default=None,
        help="clang++ executable for AST dump (auto-detected if omitted).",
    )
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text).",
    )
    parser.add_argument(
        "--include-working-tree",
        dest="include_working_tree",
        action="store_true",
        default=True,
        help=(
            "Include local unstaged/staged changes in changed-files set, "
            "in addition to <baseline>...<head> (default: enabled)."
        ),
    )
    parser.add_argument(
        "--no-include-working-tree",
        dest="include_working_tree",
        action="store_false",
        help="Disable working-tree inclusion and only analyze <baseline>...<head>.",
    )
    return parser.parse_args()


def git_repo_root() -> Path:
    proc = run_command(["git", "rev-parse", "--show-toplevel"], cwd=Path.cwd())
    return Path(proc.stdout.strip()).resolve()


def resolve_compile_commands(repo_root: Path, explicit_path: str | None) -> Path:
    if explicit_path:
        compile_path = Path(explicit_path)
        if not compile_path.is_absolute():
            compile_path = (repo_root / compile_path).resolve()
        if not compile_path.exists():
            raise FileNotFoundError(f"compile_commands.json not found: {compile_path}")
        return compile_path

    candidates = sorted(
        repo_root.glob("build/**/compile_commands.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            "No compile_commands.json found under build/**. "
            "Configure with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON first."
        )
    return candidates[0].resolve()


def load_compile_commands(
    compile_commands_path: Path,
    repo_root: Path,
) -> list[CompileCommandEntry]:
    entries: list[CompileCommandEntry] = []
    data = json.loads(compile_commands_path.read_text(encoding="utf-8"))
    for raw_entry in data:
        directory = Path(raw_entry["directory"]).resolve()

        raw_source = Path(raw_entry["file"])
        if raw_source.is_absolute():
            source = raw_source.resolve()
        else:
            source = (directory / raw_source).resolve()

        if not is_project_source(source, repo_root):
            continue

        if "arguments" in raw_entry:
            arguments = [str(arg) for arg in raw_entry["arguments"]]
        else:
            arguments = shlex.split(raw_entry["command"])

        output = infer_output_path(arguments, directory)
        target = infer_cmake_target_from_output(output)

        entries.append(
            CompileCommandEntry(
                directory=directory,
                source=source,
                arguments=arguments,
                output=output,
                target=target,
            )
        )
    return entries


def infer_output_path(arguments: list[str], directory: Path) -> Path | None:
    output_token: str | None = None
    for idx, arg in enumerate(arguments):
        if arg == "-o" and idx + 1 < len(arguments):
            output_token = arguments[idx + 1]
        elif arg.startswith("-o") and len(arg) > 2:
            output_token = arg[2:]
        elif arg.startswith("/Fo") and len(arg) > 3:
            output_token = arg[3:]

    if output_token is None:
        return None

    out_path = Path(output_token)
    if out_path.is_absolute():
        return out_path.resolve()
    return (directory / out_path).resolve()


def infer_cmake_target_from_output(output: Path | None) -> str | None:
    if output is None:
        return None
    parts = output.parts
    for index, part in enumerate(parts):
        if part == "CMakeFiles" and index + 1 < len(parts):
            target_part = parts[index + 1]
            if target_part.endswith(".dir"):
                return target_part[: -len(".dir")]
            return target_part
    return None


def git_changed_files(repo_root: Path, baseline: str, head: str) -> set[Path]:
    diff_range = f"{baseline}...{head}"
    proc = run_command(["git", "diff", "--name-only", diff_range], cwd=repo_root)
    changed_files: set[Path] = set()
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        changed_files.add((repo_root / line).resolve())
    return changed_files


def git_working_tree_changed_files(repo_root: Path) -> set[Path]:
    changed_files: set[Path] = set()
    commands = [
        ["git", "diff", "--name-only"],
        ["git", "diff", "--name-only", "--cached"],
    ]
    for cmd in commands:
        proc = run_command(cmd, cwd=repo_root)
        for line in proc.stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            changed_files.add((repo_root / line).resolve())
    return changed_files


def parse_changed_lines_from_diff_text(
    diff_text: str,
    repo_root: Path,
) -> dict[Path, set[int]]:
    changed_lines: dict[Path, set[int]] = defaultdict(set)

    current_file: Path | None = None
    in_hunk = False
    new_line = 0

    for raw_line in diff_text.splitlines():
        if raw_line.startswith("+++ "):
            file_token = raw_line[4:].strip()
            if file_token == "/dev/null":
                current_file = None
                in_hunk = False
                continue
            if file_token.startswith("b/"):
                file_token = file_token[2:]
            current_file = (repo_root / file_token).resolve()
            in_hunk = False
            continue

        hunk_match = DIFF_HUNK_RE.match(raw_line)
        if hunk_match:
            in_hunk = current_file is not None
            new_line = int(hunk_match.group(1))
            continue

        if not in_hunk or current_file is None:
            continue

        if raw_line.startswith("+") and not raw_line.startswith("+++"):
            changed_lines[current_file].add(new_line)
            new_line += 1
            continue

        if raw_line.startswith("-") and not raw_line.startswith("---"):
            continue

        if raw_line.startswith(" "):
            new_line += 1
            continue

    return changed_lines


def git_changed_line_map(
    repo_root: Path,
    baseline: str,
    head: str,
    include_working_tree: bool,
) -> dict[Path, set[int]]:
    changed_lines: dict[Path, set[int]] = defaultdict(set)

    proc = run_command(
        ["git", "diff", "--unified=0", f"{baseline}...{head}"],
        cwd=repo_root,
    )
    baseline_map = parse_changed_lines_from_diff_text(proc.stdout, repo_root)
    for path, lines in baseline_map.items():
        changed_lines[path].update(lines)

    if include_working_tree:
        for cmd in (
            ["git", "diff", "--unified=0"],
            ["git", "diff", "--cached", "--unified=0"],
        ):
            wt_proc = run_command(cmd, cwd=repo_root)
            wt_map = parse_changed_lines_from_diff_text(wt_proc.stdout, repo_root)
            for path, lines in wt_map.items():
                changed_lines[path].update(lines)

    return changed_lines


def extract_changed_symbol_names_from_file(
    file_path: Path,
    changed_lines: set[int],
) -> set[str]:
    if not changed_lines or not file_path.exists():
        return set()

    lines = file_path.read_text(encoding="utf-8", errors="replace").splitlines()
    symbols: set[str] = set()

    line_index = 1
    max_line = len(lines)
    while line_index <= max_line:
        line = lines[line_index - 1]
        match = CPP_FUNCTION_START_RE.match(line)
        if not match:
            line_index += 1
            continue

        symbol_name = match.group(1)
        start_line = line_index
        brace_depth = line.count("{") - line.count("}")
        end_line = start_line

        while brace_depth > 0 and end_line < max_line:
            end_line += 1
            body_line = lines[end_line - 1]
            brace_depth += body_line.count("{") - body_line.count("}")

        if any(start_line <= line_no <= end_line for line_no in changed_lines):
            symbols.add(symbol_name)

        line_index = end_line + 1

    return symbols


def collect_changed_symbol_names(
    changed_line_map: dict[Path, set[int]],
) -> set[str]:
    symbol_names: set[str] = set()
    for file_path, changed_lines in changed_line_map.items():
        symbol_names.update(
            extract_changed_symbol_names_from_file(file_path, changed_lines)
        )
    return symbol_names


def clean_command_for_dependency_scan(arguments: list[str]) -> list[str]:
    cleaned: list[str] = []
    skip_next = False
    flags_with_value = {
        "-o",
        "-MF",
        "-MT",
        "-MQ",
        "-MJ",
        "-Xclang",
    }
    standalone_drop = {
        "-c",
        "-MD",
        "-MMD",
        "-MP",
        "-MM",
        "-M",
        "-E",
        "-S",
    }

    index = 0
    while index < len(arguments):
        arg = arguments[index]
        if skip_next:
            skip_next = False
            index += 1
            continue

        if arg in flags_with_value:
            skip_next = True
            index += 1
            continue
        if arg in standalone_drop:
            index += 1
            continue
        if arg.startswith("-o") and len(arg) > 2:
            index += 1
            continue
        if arg.startswith("-MF") and len(arg) > 3:
            index += 1
            continue
        if arg.startswith("-MT") and len(arg) > 3:
            index += 1
            continue
        if arg.startswith("-MQ") and len(arg) > 3:
            index += 1
            continue
        if arg.startswith("-MJ") and len(arg) > 3:
            index += 1
            continue

        cleaned.append(arg)
        index += 1

    return cleaned


def parse_makefile_dependencies(stdout_text: str) -> list[str]:
    flattened = stdout_text.replace("\\\n", " ").replace("\n", " ")
    if ":" not in flattened:
        return []
    dep_payload = flattened.split(":", 1)[1].strip()
    if not dep_payload:
        return []
    return shlex.split(dep_payload)


def compute_tu_dependencies(entry: CompileCommandEntry) -> None:
    dep_cmd = clean_command_for_dependency_scan(entry.arguments)
    if not dep_cmd:
        entry.dep_error = "Empty compile command after sanitization"
        entry.dependencies = {entry.source}
        return

    dep_cmd.extend(["-MM", "-MF", "-", "-MT", "__pixie_tu__"])
    source_arg = str(entry.source)
    if source_arg not in dep_cmd:
        dep_cmd.append(source_arg)

    try:
        proc = run_command(dep_cmd, cwd=entry.directory, check=False)
    except FileNotFoundError as exc:
        entry.dep_error = str(exc)
        entry.dependencies = {entry.source}
        return

    dependencies: set[Path] = {entry.source}
    if proc.returncode != 0:
        stderr = proc.stderr.strip()
        entry.dep_error = (
            stderr if stderr else f"Dependency scan failed ({proc.returncode})"
        )
        entry.dependencies = dependencies
        return

    for dep in parse_makefile_dependencies(proc.stdout):
        dep_path = Path(dep)
        resolved = (
            dep_path.resolve()
            if dep_path.is_absolute()
            else (entry.directory / dep_path).resolve()
        )
        dependencies.add(resolved)

    entry.dependencies = dependencies


def is_build_infra_change(repo_root: Path, changed: set[Path]) -> bool:
    for path in changed:
        if path.name in BUILD_INFRA_FILES:
            return True
        try:
            rel = path.relative_to(repo_root)
        except ValueError:
            continue
        rel_text = rel.as_posix()
        if rel_text.startswith("cmake/"):
            return True
    return False


def identify_benchmark_targets(
    entries: list[CompileCommandEntry], repo_root: Path
) -> set[str]:
    benchmark_targets: set[str] = set()
    targets_present = {entry.target for entry in entries if entry.target}
    for entry in entries:
        if entry.target is None:
            continue
        try:
            rel = entry.source.relative_to(repo_root)
            rel_text = rel.as_posix()
        except ValueError:
            rel_text = entry.source.as_posix()

        if rel_text.startswith("src/benchmarks/"):
            benchmark_targets.add(entry.target)

    benchmark_targets.update(targets_present.intersection(KNOWN_BENCHMARK_TARGETS))
    return benchmark_targets


def is_benchmark_source(source: Path, repo_root: Path) -> bool:
    try:
        rel_text = source.relative_to(repo_root).as_posix()
    except ValueError:
        return False
    return rel_text.startswith("src/benchmarks/")


def dedupe_entries_by_target_source(
    entries: list[CompileCommandEntry],
) -> list[CompileCommandEntry]:
    deduped: list[CompileCommandEntry] = []
    seen: set[tuple[str | None, Path]] = set()
    for entry in entries:
        key = (entry.target, entry.source)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(entry)
    return deduped


def discover_clangxx(explicit: str | None) -> str:
    if explicit:
        return explicit

    candidates = [
        "clang++",
        "clang++-19",
        "clang++-18",
        "clang++-17",
        "clang++-16",
    ]
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise FileNotFoundError(
        "clang++ was not found on PATH. Provide --clangxx to select a clang compiler."
    )


def clean_command_for_ast(arguments: list[str], clangxx: str) -> list[str]:
    cleaned = clean_command_for_dependency_scan(arguments)
    if not cleaned:
        return []
    cleaned[0] = clangxx
    cleaned.extend(["-Xclang", "-ast-dump=json", "-fsyntax-only"])
    return cleaned


def normalize_path_candidate(path_text: str | None, working_dir: Path) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text)
    if path.is_absolute():
        return path.resolve()
    return (working_dir / path).resolve()


def file_from_loc(loc: dict[str, Any] | None, working_dir: Path) -> Path | None:
    if not isinstance(loc, dict):
        return None
    if "file" in loc:
        return normalize_path_candidate(str(loc["file"]), working_dir)
    for nested_key in ("spellingLoc", "expansionLoc", "includedFrom"):
        nested_loc = loc.get(nested_key)
        if isinstance(nested_loc, dict):
            resolved = file_from_loc(nested_loc, working_dir)
            if resolved is not None:
                return resolved
    return None


def iter_ast_nodes(node: Any):
    if isinstance(node, dict):
        yield node
        inner = node.get("inner", [])
        if isinstance(inner, list):
            for child in inner:
                yield from iter_ast_nodes(child)
    elif isinstance(node, list):
        for item in node:
            yield from iter_ast_nodes(item)


def referenced_decl_file(node: dict[str, Any], working_dir: Path) -> Path | None:
    referenced = node.get("referencedDecl")
    if not isinstance(referenced, dict):
        return None
    return file_from_loc(referenced.get("loc"), working_dir)


def node_references_changed_symbol(
    node: dict[str, Any],
    changed_symbol_names: set[str],
) -> bool:
    if not changed_symbol_names:
        return False

    for subnode in iter_ast_nodes(node):
        if not isinstance(subnode, dict):
            continue

        kind = subnode.get("kind")
        if kind == "MemberExpr":
            member_name = subnode.get("name")
            if isinstance(member_name, str) and member_name in changed_symbol_names:
                return True

        if kind == "DeclRefExpr":
            ref_decl = subnode.get("referencedDecl")
            if not isinstance(ref_decl, dict):
                continue
            ref_name = ref_decl.get("name")
            if isinstance(ref_name, str) and ref_name in changed_symbol_names:
                return True

    return False


def call_expr_callee_name(call_expr: dict[str, Any]) -> str | None:
    for node in iter_ast_nodes(call_expr):
        if not isinstance(node, dict):
            continue
        if node.get("kind") != "DeclRefExpr":
            continue
        referenced = node.get("referencedDecl")
        if isinstance(referenced, dict) and isinstance(referenced.get("name"), str):
            return referenced["name"]
    return None


def string_literals_in_node(node: dict[str, Any]) -> list[str]:
    values: list[str] = []
    for cur in iter_ast_nodes(node):
        if not isinstance(cur, dict):
            continue
        if cur.get("kind") != "StringLiteral":
            continue
        value = cur.get("value")
        if isinstance(value, str):
            if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
                value = value[1:-1]
            values.append(value)
    return values


def benchmark_names_from_source(source: Path) -> set[str]:
    names: set[str] = set()
    if not source.exists():
        return names
    text = source.read_text(encoding="utf-8", errors="replace")
    for match in re.finditer(r"BENCHMARK\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)", text):
        names.add(match.group(1))
    for match in re.finditer(r"register_op\(\s*\"([^\"]+)\"", text):
        names.add(match.group(1))
    return names


def ast_analyze_entry(
    entry: CompileCommandEntry,
    changed_files: set[Path],
    changed_symbol_names: set[str],
    clangxx: str,
) -> AstImpactResult:
    result = AstImpactResult()

    ast_cmd = clean_command_for_ast(entry.arguments, clangxx)
    if not ast_cmd:
        result.ast_error = "Failed to build AST command"
        return result

    try:
        proc = run_command(ast_cmd, cwd=entry.directory, check=False)
    except FileNotFoundError as exc:
        result.ast_error = str(exc)
        return result

    if proc.returncode != 0:
        stderr = proc.stderr.strip()
        result.ast_error = (
            stderr if stderr else f"AST command failed ({proc.returncode})"
        )
        return result

    try:
        ast_root = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        result.ast_error = f"Invalid AST JSON: {exc}"
        return result

    function_callees: dict[str, set[str]] = defaultdict(set)
    direct_impacted_functions: set[str] = set()
    dynamic_benchmarks_by_function: dict[str, set[str]] = defaultdict(set)

    for node in iter_ast_nodes(ast_root):
        if not isinstance(node, dict):
            continue

        if node.get("kind") not in {"FunctionDecl", "CXXMethodDecl"}:
            continue

        function_name = node.get("name")
        if not isinstance(function_name, str) or not function_name:
            continue

        if function_name.startswith("BM_"):
            result.benchmark_names.add(function_name)

        function_callees.setdefault(function_name, set())

        function_loc = file_from_loc(node.get("loc"), entry.directory)
        is_directly_impacted = function_loc in changed_files
        if not is_directly_impacted:
            is_directly_impacted = node_references_changed_symbol(
                node, changed_symbol_names
            )

        for subnode in iter_ast_nodes(node):
            if not isinstance(subnode, dict):
                continue

            sub_kind = subnode.get("kind")
            if sub_kind in {"CallExpr", "CXXMemberCallExpr", "CXXOperatorCallExpr"}:
                callee = call_expr_callee_name(subnode)
                if callee:
                    function_callees[function_name].add(callee)

                if callee == "register_op":
                    literal_values = string_literals_in_node(subnode)
                    if literal_values:
                        dynamic_benchmarks_by_function[function_name].add(
                            literal_values[0]
                        )

            if not is_directly_impacted:
                ref_file = referenced_decl_file(subnode, entry.directory)
                if ref_file in changed_files:
                    is_directly_impacted = True

        if is_directly_impacted:
            direct_impacted_functions.add(function_name)

    # Reverse call-graph propagation: if a function is directly impacted,
    # every caller in this TU is impacted as well (fixed-point DFS/BFS).
    callers_of: dict[str, set[str]] = defaultdict(set)
    for caller, callees in function_callees.items():
        for callee in callees:
            callers_of[callee].add(caller)

    impacted_functions = set(direct_impacted_functions)
    stack = list(direct_impacted_functions)
    while stack:
        callee_name = stack.pop()
        for caller_name in callers_of.get(callee_name, set()):
            if caller_name in impacted_functions:
                continue
            impacted_functions.add(caller_name)
            stack.append(caller_name)

    for function_name in impacted_functions:
        if function_name.startswith("BM_"):
            result.affected_names.add(function_name)

    for function_name, names in dynamic_benchmarks_by_function.items():
        result.benchmark_names.update(names)
        if function_name in impacted_functions:
            result.affected_names.update(names)

    return result


def regex_for_benchmarks(names: set[str]) -> str | None:
    if not names:
        return None
    ordered = sorted(names)
    body = "|".join(re.escape(name) for name in ordered)
    return rf"^({body})(/|$)"


def relpath_or_abs(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def main() -> int:
    cli = parse_args()

    try:
        repo_root = git_repo_root()
        changed_files = git_changed_files(repo_root, cli.baseline, cli.head)
        if cli.include_working_tree:
            changed_files.update(git_working_tree_changed_files(repo_root))
        changed_line_map = git_changed_line_map(
            repo_root,
            cli.baseline,
            cli.head,
            cli.include_working_tree,
        )
        changed_symbol_names = collect_changed_symbol_names(changed_line_map)
        compile_commands_path = resolve_compile_commands(
            repo_root, cli.compile_commands
        )
        entries = load_compile_commands(compile_commands_path, repo_root)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        if stderr:
            print(f"error: {stderr}", file=sys.stderr)
        else:
            print(f"error: command failed: {' '.join(exc.cmd)}", file=sys.stderr)
        return 2

    target_to_entries: dict[str, list[CompileCommandEntry]] = defaultdict(list)
    source_to_entries: dict[Path, list[CompileCommandEntry]] = defaultdict(list)
    for entry in entries:
        source_to_entries[entry.source].append(entry)
        if entry.target:
            target_to_entries[entry.target].append(entry)

    benchmark_targets = identify_benchmark_targets(entries, repo_root)
    all_targets = {entry.target for entry in entries if entry.target}
    benchmark_entries = dedupe_entries_by_target_source(
        [entry for entry in entries if entry.target in benchmark_targets]
    )

    infra_change = is_build_infra_change(repo_root, changed_files)
    relevant_changed_files = {
        path
        for path in changed_files
        if is_project_source(path, repo_root)
        or path.name in BUILD_INFRA_FILES
        or relpath_or_abs(path, repo_root).startswith("cmake/")
    }
    has_header_changes = any(
        path.suffix.lower() in HEADER_EXTENSIONS for path in relevant_changed_files
    )
    benchmark_source_extensions = {".c", ".cc", ".cpp", ".cxx"}
    only_benchmark_source_changes = bool(relevant_changed_files) and all(
        is_benchmark_source(path, repo_root)
        and path.suffix.lower() in benchmark_source_extensions
        for path in relevant_changed_files
    )

    directly_affected_targets: set[str] = set()
    for changed_path in changed_files:
        for entry in source_to_entries.get(changed_path, []):
            if entry.target:
                directly_affected_targets.add(entry.target)

    dependency_scan_entries: list[CompileCommandEntry] = []
    if not infra_change and not only_benchmark_source_changes:
        if has_header_changes:
            dependency_scan_entries = dedupe_entries_by_target_source(entries)
        else:
            dependency_scan_entries = benchmark_entries

    if dependency_scan_entries:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(8, (os.cpu_count() or 4))
        ) as pool:
            list(pool.map(compute_tu_dependencies, dependency_scan_entries))

    affected_targets: set[str] = set(directly_affected_targets)
    for entry in dependency_scan_entries:
        has_changed_dependency = any(dep in changed_files for dep in entry.dependencies)
        if has_changed_dependency and entry.target:
            affected_targets.add(entry.target)

    if infra_change:
        affected_targets.update(all_targets)

    dependency_impacted_benchmark_targets = affected_targets.intersection(
        benchmark_targets
    )
    impacted_benchmark_entries = [
        entry
        for entry in benchmark_entries
        if entry.target in dependency_impacted_benchmark_targets
    ]

    ast_errors: dict[str, str] = {}
    benchmark_target_to_names: dict[str, set[str]] = defaultdict(set)
    benchmark_target_to_affected: dict[str, set[str]] = defaultdict(set)
    warnings: list[str] = []
    ast_fallback_used = False
    ast_entries_scanned = 0

    if impacted_benchmark_entries:
        try:
            clangxx = discover_clangxx(cli.clangxx)
        except FileNotFoundError as exc:
            clangxx = ""
            warnings.append(str(exc))

        if not clangxx:
            ast_fallback_used = True
            for entry in impacted_benchmark_entries:
                target_name = entry.target or "<unknown-target>"
                fallback_names = benchmark_names_from_source(entry.source)
                benchmark_target_to_names[target_name].update(fallback_names)
                benchmark_target_to_affected[target_name].update(fallback_names)
        else:
            max_ast_workers = min(2, (os.cpu_count() or 2))
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=max_ast_workers
            ) as pool:
                futures = {
                    pool.submit(
                        ast_analyze_entry,
                        entry,
                        changed_files,
                        changed_symbol_names,
                        clangxx,
                    ): entry
                    for entry in impacted_benchmark_entries
                }
                ast_entries_scanned = len(futures)
                for future in concurrent.futures.as_completed(futures):
                    entry = futures[future]
                    target_name = entry.target or "<unknown-target>"
                    source_path = entry.source
                    source_is_changed = source_path in changed_files

                    try:
                        ast_result = future.result(timeout=120)
                    except Exception as exc:
                        ast_result = AstImpactResult(
                            ast_error=f"AST worker failed: {exc}"
                        )

                    if ast_result.ast_error:
                        ast_errors[relpath_or_abs(source_path, repo_root)] = (
                            ast_result.ast_error
                        )

                    benchmark_names = ast_result.benchmark_names
                    if not benchmark_names:
                        benchmark_names = benchmark_names_from_source(source_path)
                    benchmark_target_to_names[target_name].update(benchmark_names)

                    if ast_result.affected_names:
                        benchmark_target_to_affected[target_name].update(
                            ast_result.affected_names
                        )
                    elif source_is_changed or ast_result.ast_error:
                        benchmark_target_to_affected[target_name].update(
                            benchmark_names
                        )
                        if benchmark_names:
                            ast_fallback_used = True

    if infra_change and benchmark_targets:
        for target_name in sorted(benchmark_targets):
            for entry in target_to_entries.get(target_name, []):
                names = benchmark_names_from_source(entry.source)
                benchmark_target_to_names[target_name].update(names)
                benchmark_target_to_affected[target_name].update(names)

    if infra_change:
        affected_benchmark_targets = sorted(benchmark_targets)
    else:
        affected_benchmark_targets = sorted(
            target for target, names in benchmark_target_to_affected.items() if names
        )

    all_affected_benchmarks: set[str] = set()
    for names in benchmark_target_to_affected.values():
        all_affected_benchmarks.update(names)

    dep_scan_failures = {
        relpath_or_abs(entry.source, repo_root): entry.dep_error
        for entry in dependency_scan_entries
        if entry.dep_error
    }

    scope_mode = "normal"
    if infra_change:
        scope_mode = "infra_fallback"
    elif ast_fallback_used:
        scope_mode = "ast_fallback"

    report: dict[str, Any] = {
        "baseline": cli.baseline,
        "head": cli.head,
        "include_working_tree": cli.include_working_tree,
        "changed_symbols": sorted(changed_symbol_names),
        "compile_commands": relpath_or_abs(compile_commands_path, repo_root),
        "changed_files": sorted(
            relpath_or_abs(path, repo_root) for path in changed_files
        ),
        "affected_targets": sorted(affected_targets),
        "affected_benchmark_targets": affected_benchmark_targets,
        "affected_benchmarks": {
            target: sorted(names)
            for target, names in sorted(benchmark_target_to_affected.items())
            if names
        },
        "suggested_filter_regex": regex_for_benchmarks(all_affected_benchmarks),
        "dependency_entries_scanned": len(dependency_scan_entries),
        "benchmark_entries_scanned": len(benchmark_entries),
        "ast_entries_scanned": ast_entries_scanned,
        "scope_mode": scope_mode,
        "dependency_scan_failures": dep_scan_failures,
        "ast_failures": ast_errors,
        "warnings": warnings,
    }

    if cli.format == "json":
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    print(f"Baseline: {cli.baseline}")
    print(f"Head: {cli.head}")
    print(f"Compile commands: {report['compile_commands']}")
    print(f"Scope mode: {report['scope_mode']}")
    print(
        "Scan counts: "
        f"dependency={report['dependency_entries_scanned']}, "
        f"benchmark={report['benchmark_entries_scanned']}, "
        f"ast={report['ast_entries_scanned']}"
    )
    print("")

    print(f"Changed files ({len(report['changed_files'])}):")
    for item in report["changed_files"]:
        print(f"- {item}")
    if not report["changed_files"]:
        print("- none")

    print("")
    print(f"Affected targets ({len(report['affected_targets'])}):")
    for item in report["affected_targets"]:
        print(f"- {item}")
    if not report["affected_targets"]:
        print("- none")

    print("")
    print(f"Affected benchmark targets ({len(report['affected_benchmark_targets'])}):")
    for item in report["affected_benchmark_targets"]:
        print(f"- {item}")
    if not report["affected_benchmark_targets"]:
        print("- none")

    print("")
    print("Affected benchmark functions:")
    if report["affected_benchmarks"]:
        for target, names in report["affected_benchmarks"].items():
            print(f"- {target}:")
            for name in names:
                print(f"  - {name}")
    else:
        print("- none")

    print("")
    print("Suggested --benchmark_filter regex:")
    print(report["suggested_filter_regex"] or "none")

    if dep_scan_failures:
        print("")
        print("Dependency scan failures:")
        for source, error in dep_scan_failures.items():
            print(f"- {source}: {error}")

    if ast_errors:
        print("")
        print("AST failures:")
        for source, error in ast_errors.items():
            print(f"- {source}: {error}")

    if warnings:
        print("")
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
