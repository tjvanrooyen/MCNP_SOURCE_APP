"""Data-visualisation of the ionising-radiation emission source term.

One overview graph, ``<out_dir>/<stem>_ALL.{svg,pdf,png}``, puts every particle
class on a common energy axis:

* discrete classes as **lollipops** (stem + marker), read on the LEFT axis in
  particles per nuclear transformation (nt);
* the continuous beta spectra as thick **curves**, read on the RIGHT axis in
  particles per nt per MeV.

The photon (DEP) and electron (DEE) classes are drawn from the **binned** source
term: the spectroscopic versions carry tens of thousands of lines that overlap
into an unreadable smear.  Alphas, recoils, neutrons and betas are drawn from
their spectroscopically exact data.

Fission fragments (LFF/HFF, heavy ions at roughly 60-110 MeV) are deliberately
NOT graphed: they are irrelevant from a shielding perspective.  They remain, of
course, in the SDEF decks.

Style
-----
Colours, marker shapes, the line-count scaling of stem width and marker size,
the dashed grid and the black figure frame all follow the reference script
``GraphBank/Plots_RadioNuclideEmissions_SFN_DEP_DEE_BEN_BEP_DEA_DER_DEFF.py``
(see ``GraphBank/SKILLS.md``).

Energy axis
-----------
LINEAR by default.  For a LOGARITHMIC energy axis either change :data:`X_SCALE`
below, or for a single run::

    python MCNP_SOURCE_APP.py --graph-xscale log

The upper energy limit is automatic: a round number just above the
highest-energy emission that is visible, never above 15 MeV.  Fix it with
``--graph-energy-max MEV``.

Nothing graphed is silently absent: every class that carries yield is listed in
the legend with its total yield per nt, and a class with no visible mark (beyond
the energy axis, or below the yield floor) is flagged there as off-scale.

Matplotlib is imported lazily; any failure degrades to a warning so a missing
plotting backend never blocks SDEF generation.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass

from .sourceterm import SourceTerm

# =============================================================================
# ENERGY-AXIS SCALE -- change this ONE line to switch the default:
#       "linear"   (default)        or        "log"
# The command-line option  --graph-xscale {linear,log}  overrides it per run.
# =============================================================================
X_SCALE = "linear"

ENERGY_MIN_LOG_MEV = 1.0e-5   # left edge of a log energy axis
ENERGY_CAP_MEV = 15.0         # the automatic upper energy limit never exceeds this
YIELD_FLOOR = 1.0e-12         # the yield axis never extends below this

FIGSIZE = (16.0, 11.5)        # inches
PNG_DPI = 600                 # SVG and PDF are pure vector; only the PNG is raster

# --- Reference-script styling ------------------------------------------------
CURVE_WIDTH = 6.0             # continuous beta spectra
MARKER_EDGE_WIDTH = 2.0
FRAME_WIDTH = 3.0             # black border round the whole figure
GRID_MAJOR = dict(color="xkcd:charcoal", linestyle="dashed", linewidth=1.10, alpha=0.55)


@dataclass(frozen=True)
class Lollipop:
    """Style of one discrete (stem + marker) class."""

    attr: str        # SourceTerm attribute
    symbol: str      # mathtext symbol shown in the legend
    name: str        # plain-language class name
    marker: str
    face: str
    edge: str
    stem: str


@dataclass(frozen=True)
class Curve:
    """Style of one continuous (beta) class."""

    attr: str
    symbol: str
    name: str
    colour: str


# Draw order = tuple order: later classes sit on top where lines coincide.
# Fission fragments are intentionally absent (see the module docstring).
LOLLIPOPS = (
    Lollipop("dee", r"$Q_{e}(E)$", "Discrete electrons, DEE (binned)",
             "o", "xkcd:neon red", "xkcd:purply blue", "xkcd:violet"),
    Lollipop("photon", r"$Q_{\gamma}(E)$", "Photons, DEP (binned)",
             "D", "xkcd:aqua blue", "xkcd:light royal blue", "xkcd:true blue"),
    Lollipop("recoil", r"$Q_{\mathrm{rec}}(E)$", "Alpha-recoil nuclei, DER",
             "h", "xkcd:grey", "xkcd:black", "xkcd:dark grey"),
    Lollipop("alpha", r"$Q_{\alpha}(E)$", "Alpha particles, DEA",
             "o", "xkcd:rusty red", "xkcd:claret", "xkcd:rust brown"),
    Lollipop("neutron", r"$Q_{n}(E)$", "Neutrons, N",
             "s", "xkcd:spring green", "xkcd:dark green", "xkcd:emerald green"),
)

CURVES = (
    Curve("ben", r"$Q_{\beta^{-}}(E)$", "Beta-minus electrons, BEN", "xkcd:blue/purple"),
    Curve("bep", r"$Q_{\beta^{+}}(E)$", "Beta-plus positrons, BEP", "xkcd:terracotta"),
)

# Legend order: by physical importance, not by draw order.
LEGEND_ORDER = ("photon", "dee", "ben", "bep", "alpha", "recoil", "neutron")


# --- Line-count scaling, as tuned by tests in the reference script ----------
def stem_width(n_lines: int) -> float:
    """Stem line width (pt): thick for a few lines, thinner as they crowd."""
    for limit, lw in ((4000, 1.0), (1500, 1.5), (750, 2.0), (500, 2.5), (250, 3.0), (50, 3.5)):
        if n_lines >= limit:
            return lw
    return 4.0


def marker_size(n_lines: int) -> float:
    """Lollipop marker size (pt): large for a few lines, smaller as they crowd."""
    for limit, ms in ((500, 10), (250, 11), (50, 12), (25, 13), (10, 14)):
        if n_lines >= limit:
            return ms
    return 15


# --- Public entry point -------------------------------------------------------
def plot_source_term(
    st: SourceTerm,
    out_dir: str,
    *,
    stem: str = "IREST",
    title: str = "",
    mixture=None,
    energy_max: float | None = None,
    x_scale: str | None = None,
) -> list[str]:
    """Render ``st`` to ``<out_dir>/<stem>_ALL.{svg,pdf,png}``; return the paths.

    Pass a *binned* source term so the photon and electron classes are tidy.

    ``mixture``    -- optional inventory items (``.nuclide``, ``.raw_name``,
                      ``.activity``) for the composition subtitle.
    ``energy_max`` -- upper energy limit in MeV; ``None`` chooses it
                      automatically.  Raise it (e.g. 120) to show fission fragments.
    ``x_scale``    -- ``"linear"`` or ``"log"``; ``None`` uses :data:`X_SCALE`.
    """
    x_scale = (x_scale or X_SCALE).lower()
    if x_scale not in ("linear", "log"):
        raise ValueError(f"x_scale must be 'linear' or 'log', not {x_scale!r}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:                              # pragma: no cover
        print(f"  ! Plotting skipped (matplotlib unavailable): {exc}")
        return []

    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, f"{stem}_ALL")
    written: list[str] = []
    with plt.rc_context(_house_style()):
        fig = _draw(plt, st, title, mixture, energy_max, x_scale)
        for ext, dpi in (("svg", None), ("pdf", None), ("png", PNG_DPI)):
            path = f"{base}.{ext}"
            fig.savefig(path, dpi=dpi, facecolor="white", edgecolor="black")
            written.append(path)
        plt.close(fig)
    return written


# --- Drawing -----------------------------------------------------------------
def _draw(plt, st, title, mixture, energy_max, x_scale):
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FormatStrFormatter, LogLocator, MultipleLocator

    y_lo, y_hi = _yield_range(st, energy_max)
    e_hi = energy_max or _auto_energy_max(st, y_lo, x_scale)
    if x_scale == "linear":
        x_lo, x_hi = -0.012 * e_hi, e_hi                  # small margin: no marker clipped at 0
    else:
        x_lo, x_hi = ENERGY_MIN_LOG_MEV / 1.2, e_hi

    fig, ax = plt.subplots(figsize=FIGSIZE, layout="constrained")
    fig.patch.set_linewidth(FRAME_WIDTH)
    fig.patch.set_edgecolor("black")

    handles: dict[str, object] = {}
    labels: dict[str, str] = {}

    # Discrete classes: only the lines inside the plotted window are drawn, with
    # stem width and marker size scaled to how many of them there are.
    for s in LOLLIPOPS:
        d = getattr(st, s.attr)
        pts = [(e, y) for e, y in zip(d.energy, d.yield_) if y > 0.0]
        if not pts:
            continue
        marker = dict(marker=s.marker, markerfacecolor=s.face, markeredgecolor=s.edge,
                      markeredgewidth=MARKER_EDGE_WIDTH)
        vis = [(e, y) for e, y in pts if x_lo <= e <= x_hi and y >= y_lo]
        if vis:
            es, ys = zip(*vis)
            cont = ax.stem(es, ys, bottom=y_lo)
            cont.markerline.set(**marker, markersize=marker_size(len(vis)),
                                zorder=6, clip_on=False)
            cont.stemlines.set(color=s.stem, linewidth=stem_width(len(vis)),
                               linestyle="solid", capstyle="butt", zorder=3)
            cont.baseline.set_visible(False)
            handles[s.attr] = cont
        else:                         # emitted, but nothing lands in the window
            handles[s.attr] = Line2D([], [], linestyle="none", markersize=marker_size(1),
                                     **marker)
        es, ys = zip(*pts)
        note = _offscale_note(es, ys, x_lo, x_hi, y_lo)
        labels[s.attr] = _legend_label(s.symbol, s.name, d.total, len(pts), note)

    # Continuous beta spectra: thick curves.
    for c in CURVES:
        b = getattr(st, c.attr)
        pts = [(e, p) for e, p in zip(b.energy, b.pdf) if p > 0.0]
        if not pts:
            continue
        es, ps = zip(*pts)
        (line,) = ax.plot(es, ps, color=c.colour, linewidth=CURVE_WIDTH,
                          solid_capstyle="round", solid_joinstyle="round", zorder=5)
        handles[c.attr] = line
        note = _offscale_note(es, ps, x_lo, x_hi, y_lo)
        labels[c.attr] = _legend_label(c.symbol, c.name, b.total, None, note)

    # Axes, ticks and grid.
    ax.set_xscale(x_scale)
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=40))
    if x_scale == "linear":
        step = _nice_step(x_hi)
        ax.xaxis.set_major_locator(MultipleLocator(step))
        ax.xaxis.set_major_formatter(FormatStrFormatter(f"%.{_decimals(step)}f"))
    else:
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=40))
    ax.minorticks_off()
    ax.grid(True, which="major", axis="both", **GRID_MAJOR)
    ax.set_axisbelow(True)
    ax.tick_params(which="major", length=9, width=1.8, pad=8)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)

    ax.set_xlabel("Energy (MeV)", labelpad=12)
    ax.set_ylabel(r"Discrete: particle yield  $\left(\mathrm{nt}^{-1}\right)$", labelpad=10)

    if any(c.attr in handles for c in CURVES):
        axr = ax.twinx()
        axr.set_yscale("log")
        axr.set_ylim(y_lo, y_hi)
        axr.yaxis.set_major_locator(LogLocator(base=10.0, numticks=40))
        axr.minorticks_off()
        axr.tick_params(which="major", length=9, width=1.8, pad=8)
        for spine in axr.spines.values():
            spine.set_linewidth(2.0)
        axr.set_ylabel(r"Continuous: particle yield  "
                       r"$\left(\mathrm{nt}^{-1}\,\mathrm{MeV}^{-1}\right)$", labelpad=14)

    # Titles: headline, then the mixture composition.
    fig.suptitle(title or r"All-particle emission spectra $Q(E)$",
                 fontsize=30, fontweight="bold")
    ax.set_title(_subtitle(st, mixture), fontsize=17, color="xkcd:charcoal",
                 linespacing=1.6, pad=14)

    # Legend below the axes, so it never hides data.
    order = [k for k in LEGEND_ORDER if k in handles]
    if order:
        fig.legend([handles[k] for k in order], [labels[k] for k in order],
                   loc="outside lower center", ncols=2, frameon=True,
                   fancybox=True, framealpha=1.0, edgecolor="xkcd:charcoal",
                   borderpad=0.8, handlelength=2.6, handleheight=1.6,
                   columnspacing=2.4, labelspacing=0.7)
    return fig


# --- Axis ranges -------------------------------------------------------------
def _discrete_points(st, e_max: float):
    """(energy, yield) of every discrete line at or below ``e_max``."""
    for s in LOLLIPOPS:
        d = getattr(st, s.attr)
        for e, y in zip(d.energy, d.yield_):
            if y > 0.0 and e <= e_max:
                yield e, y


def _yield_range(st, energy_max: float | None) -> tuple[float, float]:
    """Decade-rounded yield limits: no lower than :data:`YIELD_FLOOR`."""
    e_cap = energy_max or ENERGY_CAP_MEV
    ys = [y for _, y in _discrete_points(st, e_cap)]
    peaks = list(ys)
    for c in CURVES:
        peaks.extend(p for p in getattr(st, c.attr).pdf if p > 0.0)
    shown = [y for y in ys if y >= YIELD_FLOOR]
    y_lo = 10.0 ** math.floor(math.log10(min(shown))) if shown else YIELD_FLOOR
    y_lo = max(y_lo, YIELD_FLOOR)
    y_top = max(peaks) if peaks else 1.0
    y_hi = 10.0 ** math.ceil(math.log10(y_top) + 0.15)   # a little headroom
    return y_lo, max(y_hi, 10.0 * y_lo)


def _auto_energy_max(st, y_lo: float, x_scale: str) -> float:
    """A round energy just above the highest visible emission (<= 15 MeV)."""
    tops = [e for e, y in _discrete_points(st, ENERGY_CAP_MEV) if y >= y_lo]
    for c in CURVES:
        b = getattr(st, c.attr)
        tops.extend(e for e, p in zip(b.energy, b.pdf) if p >= y_lo and e <= ENERGY_CAP_MEV)
    e_top = 1.05 * max(tops) if tops else ENERGY_CAP_MEV
    if x_scale == "log":
        k = math.floor(math.log10(e_top))
        for m in (1.0, 2.0, 5.0, 10.0):
            if m * 10.0 ** k >= e_top * (1 - 1e-9):
                return min(m * 10.0 ** k, ENERGY_CAP_MEV)
    step = _nice_step(e_top)
    return min(round(math.ceil(e_top / step - 1e-9) * step, 10), ENERGY_CAP_MEV)


def _nice_step(span: float) -> float:
    """Major tick step giving at most 12 intervals over ``span``."""
    k = math.floor(math.log10(span))
    for m in (0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0, 10.0):
        step = m * 10.0 ** k
        if span / step <= 12.0 + 1e-9:
            return step
    return 10.0 ** (k + 1)


def _decimals(step: float) -> int:
    """Decimal places needed to print multiples of ``step`` exactly."""
    for d in range(0, 7):
        if abs(step * 10 ** d - round(step * 10 ** d)) < 1e-9:
            return d
    return 6


# --- Text --------------------------------------------------------------------
def _legend_label(symbol: str, name: str, total: float, n_lines, note: str | None) -> str:
    if n_lines is None:
        detail = "continuous"
    else:
        detail = f"{n_lines} line" + ("" if n_lines == 1 else "s")
    text = f"{symbol}   {name}   ·   {detail}   ·   {total:.3E} " + r"$\mathrm{nt}^{-1}$"
    return f"{text}\n            → {note}" if note else text


def _offscale_note(energies, values, x_lo, x_hi, y_lo) -> str | None:
    """Say why a class has no visible mark; ``None`` if something is visible."""
    if any(x_lo <= e <= x_hi and v >= y_lo for e, v in zip(energies, values)):
        return None
    if min(energies) > x_hi:
        return (f"off-scale: emitted at {min(energies):.3G}–{max(energies):.3G} MeV, "
                "beyond the energy axis")
    if max(values) < y_lo:
        return f"off-scale: below the {y_lo:.0E} yield floor"
    return "off-scale: outside the plotted window"


def _nuclide_tex(item) -> str:
    """``Co-60`` -> ``$^{60}\\mathrm{Co}$`` (falls back to the raw name)."""
    nuc = getattr(item, "nuclide", None)
    if nuc is None:
        return str(getattr(item, "raw_name", item))
    meta = rf"\mathrm{{{nuc.meta}}}" if nuc.meta else ""
    return rf"$^{{{nuc.a}{meta}}}\mathrm{{{nuc.symbol}}}$"


def _percent(fraction: float) -> str:
    p = 100.0 * fraction
    return f"{p:.3G} %" if p >= 0.01 else f"{p:.1E} %"


def _subtitle(st, mixture, max_shown: int = 7) -> str:
    """The largest activity fractions of the mixture, then the total activity."""
    lines = []
    items = sorted(mixture or [], key=lambda it: -it.activity)
    total = sum(it.activity for it in items)
    if items and total > 0.0:
        shown = items if len(items) <= max_shown + 1 else items[:max_shown]
        parts = [f"{_nuclide_tex(it)} {_percent(it.activity / total)}" for it in shown]
        if len(items) > len(shown):
            parts.append(f"+ {len(items) - len(shown)} more")
        lines.append("Activity fractions:   " + "     ".join(parts))
    lines.append(f"Total activity {st.total_activity:.4E} Bq     ·     "
                 "yields per nuclear transformation (nt) of the mixture     ·     "
                 "DEP and DEE binned, all other classes spectroscopic")
    return "\n".join(lines)


# --- House style -------------------------------------------------------------
def _house_style() -> dict:
    """rcParams for crisp, presentation-grade output (applied in a context)."""
    from matplotlib import font_manager

    have_aptos = any(f.name == "Aptos" for f in font_manager.fontManager.ttflist)
    if have_aptos:                    # the reference script's typeface
        fonts = {"font.family": ["Aptos", "DejaVu Sans"],
                 "mathtext.fontset": "custom", "mathtext.rm": "Aptos",
                 "mathtext.it": "Aptos:italic", "mathtext.bf": "Aptos:bold",
                 "mathtext.sf": "Aptos", "mathtext.fallback": "stixsans"}
    else:
        fonts = {"font.family": ["DejaVu Sans"], "mathtext.fontset": "dejavusans"}
    return {
        **fonts,
        # Glyphs as outlines in the SVG, real embedded fonts in the PDF: both
        # render identically on machines that do not have the font installed.
        "svg.fonttype": "path",
        "pdf.fonttype": 42,
        "font.size": 18,
        "axes.labelsize": 24,
        "xtick.labelsize": 21,
        "ytick.labelsize": 21,
        "legend.fontsize": 16,
        "axes.edgecolor": "black",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "figure.constrained_layout.w_pad": 0.15,
        "figure.constrained_layout.h_pad": 0.15,
    }
