#!/usr/bin/env python3
"""Draw a rooted tree from a succinct tree encoding using pydot.

Supported modes (--mode):
- bp    : balanced parentheses in DFS order, one pair per node.
- louds : level-order unary degree sequence with a leading sentinel '('.
          Each node contributes d opening parentheses + ')'. BFS order.
- dfuds : depth-first unary degree sequence with a leading sentinel '('.
          Each node contributes d opening parentheses + ')'. DFS order.

The --sequence-dir flag produces one frame per node, highlighting that node.
The --highlight-node flag highlights a single node by 0-based id.

Examples:
  uv run --no-project --with pydot python scripts/draw_bp_tree.py \\
    "((()()())(()()))" --mode bp --output bp_tree.svg --format svg
  uv run --no-project --with pydot python scripts/draw_bp_tree.py \\
    "((()((()(())))))" --mode louds --sequence-dir out/louds --sequence-prefix tree
  uv run --no-project --with pydot python scripts/draw_bp_tree.py \\
    "((()((())))(()))" --mode dfuds --output dfuds_tree.png
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

import pydot


@dataclass
class Node:
    node_id: int
    parent_id: int | None
    depth: int
    children: list[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Normalisation helpers
# ---------------------------------------------------------------------------

def _normalize_parens(text: str) -> str:
    filtered = "".join(ch for ch in text if not ch.isspace())
    invalid = sorted({ch for ch in filtered if ch not in "()"})
    if invalid:
        chars = ", ".join(repr(ch) for ch in invalid)
        raise ValueError(f"sequence contains invalid characters: {chars}")
    if not filtered:
        raise ValueError("sequence is empty")
    return filtered


def _normalize_louds(text: str) -> tuple[str, str]:
    """Return (bit_sequence, display_sequence).

    Accepts '0'/'1' bits or '('/')' parentheses (where '(' = 1, ')' = 0).
    """
    filtered = "".join(ch for ch in text if not ch.isspace())
    if not filtered:
        raise ValueError("sequence is empty")
    if all(ch in "01" for ch in filtered):
        return filtered, filtered
    if all(ch in "()" for ch in filtered):
        bits = "".join("1" if ch == "(" else "0" for ch in filtered)
        return bits, filtered
    invalid = sorted({ch for ch in filtered if ch not in "01()"})
    if invalid:
        chars = ", ".join(repr(ch) for ch in invalid)
        raise ValueError(f"LOUDS sequence contains invalid characters: {chars}")
    raise ValueError("LOUDS sequence must use either only '0/1' or only '(/' and ')' tokens")


def _make_node(parent_id: int | None, depth: int, nodes: list[Node]) -> int:
    node_id = len(nodes)
    nodes.append(Node(node_id=node_id, parent_id=parent_id, depth=depth))
    if parent_id is not None:
        nodes[parent_id].children.append(node_id)
    return node_id


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_bp(sequence: str) -> list[Node]:
    """Parse a BP string into a node list."""
    seq = _normalize_parens(sequence)
    nodes: list[Node] = []
    stack: list[int] = []

    for pos, ch in enumerate(seq):
        if ch == "(":
            parent_id = stack[-1] if stack else None
            node_id = _make_node(parent_id, len(stack), nodes)
            stack.append(node_id)
        else:
            if not stack:
                raise ValueError(f"Unmatched ')' at position {pos}")
            stack.pop()

    if stack:
        raise ValueError(f"Unmatched '(' — {len(stack)} unclosed node(s)")
    return nodes


def parse_louds(sequence: str) -> list[Node]:
    """Parse a LOUDS parenthesis sequence with a leading sentinel '('."""
    bit_seq, _ = _normalize_louds(sequence)
    n = len(bit_seq)
    if n < 1 or bit_seq[0] != "1":
        raise ValueError("LOUDS sequence must start with a sentinel '1' (or '(' in paren form)")

    nodes: list[Node] = []
    pos = 1  # skip sentinel
    root_id = _make_node(None, 0, nodes)
    queue: deque[int] = deque([root_id])

    while queue:
        node_id = queue.popleft()
        while pos < n and bit_seq[pos] == "1":
            child_id = _make_node(node_id, nodes[node_id].depth + 1, nodes)
            queue.append(child_id)
            pos += 1
        if pos >= n or bit_seq[pos] != "0":
            raise ValueError(f"LOUDS node {node_id} missing terminating '0' at position {pos}")
        pos += 1  # consume '0'

    if pos != n:
        raise ValueError(f"LOUDS sequence has trailing data starting at position {pos}")
    return nodes


def parse_dfuds(sequence: str) -> list[Node]:
    """Parse a DFUDS parenthesis sequence with a leading sentinel '('."""
    seq = _normalize_parens(sequence)
    n = len(seq)
    if n < 1 or seq[0] != "(":
        raise ValueError("DFUDS sequence must start with a sentinel '('")

    nodes: list[Node] = []
    pos = 1  # skip sentinel
    root_id = _make_node(None, 0, nodes)
    queue: deque[int] = deque([root_id])

    while queue:
        node_id = queue.popleft()
        while pos < n and seq[pos] == "(":
            child_id = _make_node(node_id, nodes[node_id].depth + 1, nodes)
            queue.append(child_id)
            pos += 1
        if pos >= n or seq[pos] != ")":
            raise ValueError(f"DFUDS node {node_id} missing terminating ')' at position {pos}")
        pos += 1  # consume ')'

    if pos != n:
        raise ValueError(f"DFUDS sequence has trailing data starting at position {pos}")
    return nodes


def parse_nodes(sequence: str, mode: str) -> list[Node]:
    if mode == "bp":
        return parse_bp(sequence)
    if mode == "louds":
        return parse_louds(sequence)
    if mode == "dfuds":
        return parse_dfuds(sequence)
    raise ValueError(f"Unsupported mode: {mode}")


# ---------------------------------------------------------------------------
# Graph building
# ---------------------------------------------------------------------------

def _validate_highlight_node(nodes: list[Node], highlight_node: int | None) -> None:
    if highlight_node is None:
        return
    if highlight_node < 0 or highlight_node >= len(nodes):
        raise ValueError(
            f"highlight node {highlight_node} out of range [0, {len(nodes) - 1}]"
        )


def make_graph(
    nodes: list[Node],
    label_mode: str,
    highlight_node: int | None = None,
) -> pydot.Dot:
    _validate_highlight_node(nodes, highlight_node)

    graph = pydot.Dot(
        graph_type="digraph",
        rankdir="TB",
        dpi="220",
        nodesep="0.3",
        ranksep="0.45",
        splines="true",
        outputorder="edgesfirst",
    )
    graph.set_node_defaults(
        shape="circle",
        fontname="Helvetica",
        fontsize="20",
        width="0.7",
        height="0.7",
        penwidth="2.0",
        margin="0.06,0.04",
    )
    graph.set_edge_defaults(color="#444444", penwidth="1.8", arrowsize="0.6", arrowhead="none")

    for node in nodes:
        if label_mode == "id":
            label = str(node.node_id)
        elif label_mode == "depth":
            label = f"{node.node_id}\\nd={node.depth}"
        else:
            label = ""

        node_style: dict[str, str] = {}
        if highlight_node == node.node_id:
            node_style = {
                "color": "#c62828",
                "fontcolor": "#c62828",
                "fillcolor": "#fde7e7",
                "style": "filled,bold",
                "penwidth": "3.2",
            }

        graph.add_node(pydot.Node(str(node.node_id), label=label, **node_style))

    for node in nodes:
        if node.parent_id is not None:
            graph.add_edge(pydot.Edge(str(node.parent_id), str(node.node_id)))

    return graph


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_graph(graph: pydot.Dot, output: Path, fmt: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "dot":
        graph.write_raw(str(output))
    elif fmt == "png":
        graph.write_png(str(output))
    elif fmt == "svg":
        graph.write_svg(str(output))
    else:
        raise ValueError(f"Unsupported output format: {fmt}")


def infer_format_from_output(output: Path, explicit_format: str | None) -> str:
    if explicit_format is not None:
        return explicit_format
    ext = output.suffix.lower().lstrip(".")
    if ext in {"dot", "png", "svg"}:
        return ext
    return "png"


def write_highlight_sequence(
    nodes: list[Node],
    label_mode: str,
    output_dir: Path,
    prefix: str,
    fmt: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    digits = max(2, len(str(len(nodes) - 1)))

    for node in nodes:
        graph = make_graph(nodes, label_mode, highlight_node=node.node_id)
        output = output_dir / f"{prefix}_{node.node_id:0{digits}d}.{fmt}"
        write_graph(graph, output, fmt)
        print(f"[saved] {output} ({fmt})")

    print(f"[info] generated_frames={len(nodes)}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Draw a rooted tree from a succinct tree encoding."
    )
    parser.add_argument(
        "sequence",
        help="Input sequence interpreted according to --mode.",
    )
    parser.add_argument(
        "--mode",
        choices=["bp", "louds", "dfuds"],
        default="bp",
        help="Input encoding mode. Default: %(default)s",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tree.png"),
        help="Output file path. Default: %(default)s",
    )
    parser.add_argument(
        "--format",
        choices=["dot", "png", "svg"],
        default=None,
        help="Output format. If omitted, inferred from --output extension.",
    )
    parser.add_argument(
        "--label-mode",
        choices=["none", "id", "depth"],
        default="id",
        help="Node labels: none, numeric id, or id+depth. Default: %(default)s",
    )
    parser.add_argument(
        "--highlight-node",
        type=int,
        default=None,
        help="0-based node id to emphasise in red.",
    )
    parser.add_argument(
        "--sequence-dir",
        type=Path,
        default=None,
        help=(
            "Output directory for one-frame-per-node highlighting sequence. "
            "If set, --output is ignored."
        ),
    )
    parser.add_argument(
        "--sequence-prefix",
        default="tree",
        help="Filename prefix for sequence frames. Default: %(default)s",
    )
    args = parser.parse_args()

    nodes = parse_nodes(args.sequence, args.mode)

    if args.sequence_dir is not None:
        fmt = args.format or "png"
        write_highlight_sequence(
            nodes=nodes,
            label_mode=args.label_mode,
            output_dir=args.sequence_dir,
            prefix=args.sequence_prefix,
            fmt=fmt,
        )
        print(f"[info] mode={args.mode} nodes={len(nodes)}")
        return 0

    graph = make_graph(nodes, args.label_mode, args.highlight_node)
    fmt = infer_format_from_output(args.output, args.format)
    write_graph(graph, args.output, fmt)

    print(f"[saved] {args.output} ({fmt})")
    if args.highlight_node is not None:
        print(f"[info] highlighted_node={args.highlight_node}")
    print(f"[info] mode={args.mode} nodes={len(nodes)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
