#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMPLATE_DIR="$SCRIPT_DIR/scripts/templates"

usage() {
    cat <<EOF
Usage: $(basename "$0") [--style <style>] <paper-name>

Scaffolds a new paper and installs the appropriate journal extension.

Options:
  --style STYLE   Journal style: acm (default) or ieee

Examples:
  ./scripts/init.sh my-paper
  ./scripts/init.sh --style ieee my-paper
EOF
}

style="acm"
name=""

while [ $# -gt 0 ]; do
    case "$1" in
        --style)
            style="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Error: unknown option '$1'" >&2
            usage
            exit 1
            ;;
        *)
            name="$1"
            shift
            ;;
    esac
done

if [ -z "$name" ]; then
    echo "Error: paper name required" >&2
    usage
    exit 1
fi

case "$style" in
    acm)
        ext="quarto-journals/acm"
        ;;
    ieee)
        ext="dfolio/quarto-ieee"
        ;;
    *)
        echo "Error: unknown style '$style' (supported: acm, ieee)" >&2
        exit 1
        ;;
esac

if [ ! -f "$TEMPLATE_DIR/$style.qmd" ]; then
    echo "Error: template not found: $TEMPLATE_DIR/$style.qmd" >&2
    exit 1
fi

dest="$SCRIPT_DIR/papers/$name"

if [ -d "$dest" ]; then
    echo "Error: directory '$dest' already exists" >&2
    exit 1
fi

echo "Scaffolding new paper: $name (style: $style)"
cp -r "$SCRIPT_DIR/papers/sample-paper" "$dest"
cp "$TEMPLATE_DIR/$style.qmd" "$dest/paper.qmd"

echo "Installing $ext extension..."
cd "$dest"
quarto add "$ext" --no-prompt

echo "Created $dest"
echo "  - Edit $dest/paper.qmd to start writing"
echo "  - Update $dest/references.bib with your references"
echo "  - Render with: cd papers/$name && quarto render"
