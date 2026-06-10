#!/usr/bin/env python3
"""Draw standalone succinct tree encodings with optional token highlighting.

Supported modes:
- BP: balanced parentheses in DFS order, one parenthesis pair per node.
- LOUDS: level-order unary degree sequence, accepted as bits (1/0) or
  equivalent parentheses (() form where '(' = 1 and ')' = 0).
- DFUDS: depth-first unary degree sequence over parentheses with a single leading
  sentinel '('. Each node with d children contributes d opening parentheses
  followed by one closing parenthesis. The shown interval for each node spans
  its entire d+1-token block, anchored at the block's first token.

Examples:
  uv run --no-project --with pydot python scripts/draw_bp_representation.py \
    "((()()())(()()))" --mode bp --output report/bp_representation.png
  uv run --no-project --with pydot python scripts/draw_bp_representation.py \
    "((()((()(())))))" --mode louds --output report/louds_representation.svg
  uv run --no-project --with pydot python scripts/draw_bp_representation.py \
    "((()((())))(()))" --mode dfuds --output report/dfuds_representation.svg
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

import pydot


@dataclass
class TreeNode:
    node_id: int
    parent_id: int | None
    depth: int
    children: list[int] = field(default_factory=list)


@dataclass
class EncodedTree:
    mode: str
    sequence: str
    nodes: list[TreeNode]
    node_token_positions: dict[int, int]
    interval_targets: dict[int, int] = field(default_factory=dict)


def _normalize_parens(text: str) -> str:
    filtered = "".join(ch for ch in text if not ch.isspace())
    invalid = sorted({ch for ch in filtered if ch not in "()"})
    if invalid:
        chars = ", ".join(repr(ch) for ch in invalid)
        raise ValueError(f"sequence contains invalid parenthesis characters: {chars}")
    if not filtered:
        raise ValueError("sequence is empty")
    return filtered


def _normalize_louds(text: str) -> tuple[str, str]:
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
    raise ValueError("LOUDS sequence must use either only '0/1' or only '(' and ')' tokens")


def _validate_highlight_index(length: int, highlight_index: int | None) -> None:
    if highlight_index is None:
        return
    if highlight_index < 0 or highlight_index >= length:
        raise ValueError(
            f"highlight index {highlight_index} out of range [0, {length - 1}]"
        )


def _make_node(parent_id: int | None, depth: int, nodes: list[TreeNode]) -> int:
    node_id = len(nodes)
    nodes.append(TreeNode(node_id=node_id, parent_id=parent_id, depth=depth))
    if parent_id is not None:
        nodes[parent_id].children.append(node_id)
    return node_id


def parse_bp_encoded(bp: str) -> EncodedTree:
    sequence = _normalize_parens(bp)
    stack: list[tuple[int, int]] = []
    nodes: list[TreeNode] = []
    node_token_positions: dict[int, int] = {}
    interval_targets: dict[int, int] = {}

    for pos, ch in enumerate(sequence):
        if ch == "(":
            parent_id = stack[-1][0] if stack else None
            depth = len(stack)
            node_id = _make_node(parent_id, depth, nodes)
            node_token_positions[node_id] = pos
            stack.append((node_id, pos))
        else:
            if not stack:
                raise ValueError(f"Unmatched ')' at position {pos}")
            node_id, _ = stack.pop()
            interval_targets[node_id] = pos

    if stack:
        _, pos = stack[-1]
        raise ValueError(f"Unmatched '(' at position {pos}")

    return EncodedTree(
        mode="bp",
        sequence=sequence,
        nodes=nodes,
        node_token_positions=node_token_positions,
        interval_targets=interval_targets,
    )


def parse_louds_encoded(louds: str) -> EncodedTree:
    """Parse a LOUDS sequence.

    Convention: the sequence starts with a single sentinel '(' followed by the
    node blocks in BFS order. Each node v with d children contributes a block of
    d+1 tokens: d opening parentheses followed by one closing parenthesis. A leaf
    contributes just ')'. The sentinel is NOT part of any node's block.

    Accepts either parentheses ('(' = 1, ')' = 0) or a bit string ('1'/'0');
    in both cases the display sequence preserves the input form.

    The interval shown for node v spans its entire block [block_start, block_end].

    Example for tree 0->[1,2], 1->[3,4,5], 2->[6,7], leaves 3-7:
        ((()((()(()))))
        ^               sentinel
         ^^  ^          node0: block [1,3]
            ^^^  ^      node1: block [4,7]
    """
    bit_sequence, display_sequence = _normalize_louds(louds)
    n = len(bit_sequence)
    if n < 1 or bit_sequence[0] != "1":
        raise ValueError("LOUDS sequence must start with a sentinel '1' (or '(' in paren form)")

    nodes: list[TreeNode] = []
    node_token_positions: dict[int, int] = {}
    interval_targets: dict[int, int] = {}

    # Blocks appear in BFS order; use a FIFO.
    pos = 1  # skip sentinel at position 0
    root_id = _make_node(parent_id=None, depth=0, nodes=nodes)
    nodes_queue: deque[int] = deque([root_id])

    while nodes_queue:
        node_id = nodes_queue.popleft()
        block_start = pos

        # Read d '1' tokens — each creates one child
        while pos < n and bit_sequence[pos] == "1":
            child_id = _make_node(
                parent_id=node_id,
                depth=nodes[node_id].depth + 1,
                nodes=nodes,
            )
            nodes_queue.append(child_id)
            pos += 1

        # Expect the closing '0'
        if pos >= n or bit_sequence[pos] != "0":
            raise ValueError(
                f"LOUDS node {node_id} block starting at {block_start} "
                f"is missing a terminating '0' at position {pos}"
            )
        block_end = pos
        pos += 1

        # Shift by +1 for display because we insert a sentinel ')' at index 1
        node_token_positions[node_id] = block_start + 1
        interval_targets[node_id] = block_end + 1

    if pos != n:
        raise ValueError(f"LOUDS sequence has trailing data starting at position {pos}")

    # Build display sequence with an explicit sentinel ')' inserted after the opening '('
    # so the sentinel appears as a '()' pair rather than a bare '('.
    sentinel_ch = display_sequence[0]      # '(' or '1'
    sentinel_close = ")" if sentinel_ch == "(" else "0"
    augmented_display = sentinel_ch + sentinel_close + display_sequence[1:]

    return EncodedTree(
        mode="louds",
        sequence=augmented_display,
        nodes=nodes,
        node_token_positions=node_token_positions,
        interval_targets=interval_targets,
    )


def parse_dfuds_encoded(dfuds: str) -> EncodedTree:
    """Parse a DFUDS sequence.

    Convention: the sequence starts with a single sentinel '(' followed by the
    node blocks in DFS pre-order. Each node v with d children contributes a block
    of d+1 tokens: d opening parentheses followed by one closing parenthesis. A
    leaf contributes just ')'. The sentinel is NOT part of any node's block.

    The interval shown for node v spans its entire block [block_start, block_end],
    where block_start is the first token of the block and block_end is the closing ')'.
    For a leaf this is a single-token interval.

    Example for tree 0->[1,5], 1->[2,3,4], 5->[6,7], leaves 2,3,4,6,7:
        ((()((())))(()))
        ^               sentinel
         ^^  ^          node0: block [1,3]
            ^^^  ^      node1: block [4,7]
                  ^     node2: block [8,8]  ...
    """
    sequence = _normalize_parens(dfuds)
    n = len(sequence)
    if n < 1 or sequence[0] != "(":
        raise ValueError("DFUDS sequence must start with a sentinel '('")

    nodes: list[TreeNode] = []
    node_token_positions: dict[int, int] = {}
    interval_targets: dict[int, int] = {}

    # Blocks appear in DFS pre-order (BFS ordering of block positions in the sequence).
    # Use a FIFO so we process nodes in the same order their blocks appear.
    pos = 1  # skip sentinel at position 0
    root_id = _make_node(parent_id=None, depth=0, nodes=nodes)
    nodes_queue: deque[int] = deque([root_id])

    while nodes_queue:
        node_id = nodes_queue.popleft()
        block_start = pos

        # Read d opening parens — each creates one child
        while pos < n and sequence[pos] == "(":
            child_id = _make_node(
                parent_id=node_id,
                depth=nodes[node_id].depth + 1,
                nodes=nodes,
            )
            nodes_queue.append(child_id)
            pos += 1

        # Expect the closing ')'
        if pos >= n or sequence[pos] != ")":
            raise ValueError(
                f"DFUDS node {node_id} block starting at {block_start} "
                f"is missing a terminating ')' at position {pos}"
            )
        block_end = pos
        pos += 1

        node_token_positions[node_id] = block_start
        interval_targets[node_id] = block_end

    if pos != n:
        raise ValueError(
            f"DFUDS sequence has trailing data starting at position {pos}"
        )

    return EncodedTree(
        mode="dfuds",
        sequence=sequence,
        nodes=nodes,
        node_token_positions=node_token_positions,
        interval_targets=interval_targets,
    )


def parse_encoded_tree(sequence: str, mode: str) -> EncodedTree:
    if mode == "bp":
        return parse_bp_encoded(sequence)
    if mode == "louds":
        return parse_louds_encoded(sequence)
    if mode == "dfuds":
        return parse_dfuds_encoded(sequence)
    raise ValueError(f"Unsupported mode: {mode}")


def _mode_label(mode: str) -> str:
    return mode.upper()


def _mode_token_font(mode: str) -> str:
    return "Times New Roman" if mode != "louds" else "Helvetica"


def _add_token_row(
    graph: pydot.Dot,
    sequence: str,
    mode: str,
    with_label: bool,
    highlight_index: int | None,
) -> list[str]:
    token_subgraph = pydot.Subgraph(f"{mode}_tokens", rank="sink")
    token_names: list[str] = []
    fontname = _mode_token_font(mode)

    for i, ch in enumerate(sequence):
        is_highlighted = i == highlight_index
        token_name = f"token_{i}"
        token_names.append(token_name)
        token_subgraph.add_node(
            pydot.Node(
                token_name,
                shape="plaintext",
                label=ch,
                fontname=(f"{fontname} Bold" if is_highlighted else fontname),
                fontsize="34" if is_highlighted else "30",
                fontcolor="#c62828" if is_highlighted else "black",
                group=f"g{i}",
            )
        )

    if with_label:
        token_subgraph.add_node(
            pydot.Node(
                "encoding_label",
                shape="plaintext",
                label=_mode_label(mode),
                fontname="Times New Roman",
                fontsize="26",
            )
        )

    graph.add_subgraph(token_subgraph)

    if with_label and token_names:
        graph.add_edge(
            pydot.Edge(
                "encoding_label",
                token_names[0],
                style="invis",
                weight="100",
            )
        )

    for i in range(len(token_names) - 1):
        graph.add_edge(
            pydot.Edge(
                token_names[i],
                token_names[i + 1],
                style="invis",
                weight="100",
            )
        )

    return token_names


def _make_graph_base(splines: str = "ortho") -> pydot.Dot:
    return pydot.Dot(
        graph_type="digraph",
        rankdir="TB",
        splines=splines,
        nodesep="0.12",
        ranksep="0.32",
        dpi="220",
    )


def _add_single_level_intervals(
    graph: pydot.Dot,
    encoded: EncodedTree,
    highlight_index: int | None,
    node_shape: str,
) -> None:
    ids_subgraph = pydot.Subgraph(f"{encoded.mode}_ids", rank="same")

    for node in encoded.nodes:
        start_pos = encoded.node_token_positions[node.node_id]
        end_pos = encoded.interval_targets.get(node.node_id)
        touches_highlight = highlight_index in {start_pos, end_pos}
        id_name = f"id_{node.node_id}"

        ids_subgraph.add_node(
            pydot.Node(
                id_name,
                shape=node_shape,
                fixedsize="true",
                width="0.34",
                height="0.34",
                label=str(node.node_id),
                fontname="Times New Roman",
                fontsize="16",
                penwidth="2.2" if touches_highlight else "1.2",
                color="#c62828" if touches_highlight else "black",
                fontcolor="#c62828" if touches_highlight else "black",
                group=f"g{start_pos}",
            )
        )

        graph.add_edge(
            pydot.Edge(
                id_name,
                f"token_{start_pos}",
                style="invis",
                weight="130",
            )
        )

        if end_pos is None:
            continue

        graph.add_edge(
            pydot.Edge(
                id_name,
                f"token_{end_pos}",
                arrowhead="normal",
                arrowsize="0.45",
                penwidth="1.1",
                color="#c62828" if touches_highlight else "#444444",
                constraint="false",
            )
        )

    graph.add_subgraph(ids_subgraph)


def make_bp_graph(
    encoded: EncodedTree,
    with_label: bool,
    highlight_index: int | None = None,
) -> pydot.Dot:
    _validate_highlight_index(len(encoded.sequence), highlight_index)
    graph = _make_graph_base(splines="ortho")
    _add_token_row(graph, encoded.sequence, encoded.mode, with_label, highlight_index)
    _add_single_level_intervals(
        graph,
        encoded=encoded,
        highlight_index=highlight_index,
        node_shape="circle",
    )

    return graph


def make_louds_graph(
    encoded: EncodedTree,
    with_label: bool,
    highlight_index: int | None = None,
) -> pydot.Dot:
    _validate_highlight_index(len(encoded.sequence), highlight_index)
    graph = _make_graph_base(splines="ortho")
    _add_token_row(graph, encoded.sequence, encoded.mode, with_label, highlight_index)
    _add_single_level_intervals(
        graph,
        encoded=encoded,
        highlight_index=highlight_index,
        node_shape="box",
    )

    return graph


def make_dfuds_graph(
    encoded: EncodedTree,
    with_label: bool,
    highlight_index: int | None = None,
) -> pydot.Dot:
    _validate_highlight_index(len(encoded.sequence), highlight_index)
    graph = _make_graph_base(splines="ortho")
    _add_token_row(graph, encoded.sequence, encoded.mode, with_label, highlight_index)
    _add_single_level_intervals(
        graph,
        encoded=encoded,
        highlight_index=highlight_index,
        node_shape="box",
    )

    return graph


def make_graph(
    encoded: EncodedTree,
    with_label: bool,
    highlight_index: int | None = None,
) -> pydot.Dot:
    if encoded.mode == "bp":
        return make_bp_graph(encoded, with_label, highlight_index)
    if encoded.mode == "louds":
        return make_louds_graph(encoded, with_label, highlight_index)
    if encoded.mode == "dfuds":
        return make_dfuds_graph(encoded, with_label, highlight_index)
    raise ValueError(f"Unsupported mode: {encoded.mode}")


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


def infer_format(output: Path, explicit_format: str | None) -> str:
    if explicit_format:
        return explicit_format
    ext = output.suffix.lower().lstrip(".")
    if ext in {"dot", "png", "svg"}:
        return ext
    return "png"


def write_highlight_sequence(
    encoded: EncodedTree,
    with_label: bool,
    output_dir: Path,
    prefix: str,
    fmt: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    digits = max(2, len(str(len(encoded.sequence) - 1)))

    for i in range(len(encoded.sequence)):
        graph = make_graph(encoded, with_label, highlight_index=i)
        frame_path = output_dir / f"{prefix}_{i:0{digits}d}.{fmt}"
        write_graph(graph, frame_path, fmt)
        print(f"[saved] {frame_path} ({fmt})")

    print(f"[info] generated_frames={len(encoded.sequence)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Draw a standalone BP, LOUDS, or DFUDS representation."
    )
    parser.add_argument(
        "sequence",
        help=(
            "Input sequence interpreted according to --mode. "
            "For --mode louds, accepts either 0/1 or equivalent parentheses tokens."
        ),
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
        default=Path("tree_representation.png"),
        help="Output file path. Default: %(default)s",
    )
    parser.add_argument(
        "--format",
        choices=["dot", "png", "svg"],
        default=None,
        help="Output format. If omitted, inferred from --output extension.",
    )
    parser.add_argument(
        "--no-label",
        dest="with_label",
        action="store_false",
        default=True,
        help="Hide the encoding label on the left.",
    )
    parser.add_argument(
        "--highlight-index",
        type=int,
        default=None,
        help=(
            "0-based token index to emphasize in red. "
            "Example: --highlight-index 5"
        ),
    )
    parser.add_argument(
        "--sequence-dir",
        type=Path,
        default=None,
        help=(
            "Output directory for a full highlight sequence (one frame per token). "
            "If set, --output is ignored."
        ),
    )
    parser.add_argument(
        "--sequence-prefix",
        default="encoding_step",
        help="Filename prefix for sequence frames. Default: %(default)s",
    )
    args = parser.parse_args()

    encoded = parse_encoded_tree(args.sequence, args.mode)

    if args.sequence_dir is not None:
        fmt = args.format or "png"
        write_highlight_sequence(
            encoded=encoded,
            with_label=args.with_label,
            output_dir=args.sequence_dir,
            prefix=args.sequence_prefix,
            fmt=fmt,
        )
        print(f"[info] mode={encoded.mode} length={len(encoded.sequence)} nodes={len(encoded.nodes)}")
        return 0

    graph = make_graph(encoded, args.with_label, args.highlight_index)
    fmt = infer_format(args.output, args.format)
    write_graph(graph, args.output, fmt)

    print(f"[saved] {args.output} ({fmt})")
    if args.highlight_index is not None:
        print(f"[info] highlighted_index={args.highlight_index}")
    print(f"[info] mode={encoded.mode} length={len(encoded.sequence)} nodes={len(encoded.nodes)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
