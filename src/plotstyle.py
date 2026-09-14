"""Repo default plot palette: the project's five brand colours as the matplotlib
colour cycle, topped up to ten distinct entries.

Import it for its side effect -- ``import plotstyle`` at the top of a plotting
driver sets the default cycle, so any artist drawn without an explicit colour
(bare ``C0``..``C9``, or successive ``ax.plot`` calls) takes these colours in
order. The five brand colours are C0..C4; C5..C9 are matplotlib defaults that
do not collide, so a plot of up to ten series stays distinct.

    import plotstyle          # applies on import
    plotstyle.PALETTE[0]      # '#3DAE2B', the first brand colour

An explicit ``color=`` on a plot call still wins -- this only sets the default
that unspecified artists fall back to.
"""
import matplotlib as mpl
from cycler import cycler

# The five project colours, in cycle order (C0..C4).
BRAND_COLORS = [
    "#3DAE2B",  # green
    "#00A0DF",  # blue
    "#ED40A9",  # magenta
    "#B14FC5",  # purple
    "#F38B00",  # orange
]


def _default_tail(existing, count):
    """The first ``count`` matplotlib default colours not already in ``existing``.

    Fills C5..C9 from the stock cycle so larger plots keep distinct colours,
    while guaranteeing the brand colours stay exactly the first five.
    """
    default_colors = mpl.rcParamsDefault["axes.prop_cycle"].by_key()["color"]
    already_used = {colour.lower() for colour in existing}
    tail = []
    for colour in default_colors:
        if colour.lower() in already_used:
            continue
        tail.append(colour)
        if len(tail) == count:
            break
    return tail


# Ten-colour cycle: five brand colours followed by five non-clashing defaults.
PALETTE = BRAND_COLORS + _default_tail(BRAND_COLORS, 5)


def apply():
    """Set the repo palette as matplotlib's default colour cycle."""
    mpl.rcParams["axes.prop_cycle"] = cycler(color=PALETTE)


apply()
