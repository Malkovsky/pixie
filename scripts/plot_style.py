"""Shared Matplotlib defaults for Pixie's local plotting scripts.

Call apply_style() once, and use line_style() for consistent series colors.
Tune this file to change the appearance of plots using these helpers.
"""

import matplotlib as mpl

PANEL_SIZE = (4.6, 6.6)
PANEL_COLUMNS = 3
LEGEND_COLUMNS = 4
DPI = 180
RC_PARAMS = {
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "lines.linewidth": 1.4,
    "lines.markersize": 4,
    "legend.fontsize": 10,
    "legend.frameon": False,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
}
# Stable identifiers can also be used by RMQ, RmM, and other local plots.
SERIES = {
    "baseline": ("#0072B2", "o"),
    "compressed": ("#009E73", "s"),
    "pasta": ("#E69F00", "^"),
    "3star": ("#CC79A7", "D"),
    "succinct_fast": ("#D55E00", "v"),
    "orzo": ("#6F4E9C", "P"),
    "sdsl": ("#56B4E9", "h"),
    "spider": ("#555555", "X"),
}


def apply_style():
    mpl.rcParams.update(RC_PARAMS)


def line_style(series, index=0):
    """Return plot kwargs; unknown series use a deterministic ordered fallback."""
    styles = tuple(SERIES.values())
    color, marker = SERIES.get(series, styles[index % len(styles)])
    return {"color": color, "marker": marker}


def log_y_axis(ax, percent=False):
    """Use base-two ticks with plain numeric labels on a logarithmic y axis."""
    from matplotlib.ticker import FuncFormatter, LogLocator
    from numpy import format_float_positional

    ax.set_yscale("log", base=2)
    ax.yaxis.set_major_locator(LogLocator(base=2, numticks=100))

    def label(value, position):
        number = (f"{value:.0f}" if value >= 1 else
                  format_float_positional(value, precision=6, fractional=False, trim="-"))
        return number + ("%" if percent else "")

    ax.yaxis.set_major_formatter(FuncFormatter(label))


def build_time_y_axis(ax):
    """Label every time decade, with unlabeled 2–9 subdivisions."""
    from matplotlib.ticker import EngFormatter, FuncFormatter, LogLocator, NullFormatter

    ax.set_yscale("log", base=10)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=100))
    units = EngFormatter(unit="s", sep="")
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: units(value / 1000).replace("µ", "u"))
    )
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(2, 10), numticks=100))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(which="minor", axis="y", alpha=0.12, linewidth=0.4)
