"""Data-visualisation of the ionising-radiation emission source term.

A log-log overview of every particle class on a common energy axis.  By design
the photon (DEP) and electron (DEE) classes are drawn from the **binned** source
term: the spectroscopic versions carry tens of thousands of lines that overlap
into an unreadable smear, whereas the binned versions are clean.  All other
classes -- alphas, neutrons, recoils, beta-minus and beta-plus -- are drawn from
their spectroscopically exact data.

Conventions (per ``GraphBank/SKILLS_GraphBank.md``)
---------------------------------------------------
* Energy axis 1e-5 .. 15 MeV, log scale.  ``energy_max`` is adjustable so the
  fission fragments, which sit near 90 MeV, can be brought into view -- from the
  command line that is ``--graph-energy-max 200``.
* Yield axis lower limit pinned to 1e-12.
* Colours: photon = blue, DEE = violet, BEN = light royal blue, BEP = rouge,
  alpha = brown, neutron = green, recoil = black, fission fragment = royal
  purple.

Readability
-----------
Two habits keep this plot legible with ~1300 lines on it.  The legend sits
*below* the axes rather than on them, so no data is ever hidden behind it.  And
because several classes are routinely off-scale -- neutrons below the 1e-12
yield floor, fission fragments far to the right of 15 MeV -- a caption beneath
the legend names them and gives their yields, so a class that has no visible
mark is still reported rather than silently absent.

Matplotlib is imported lazily; any failure degrades to a warning so a missing
plotting backend never blocks SDEF generation.
"""

from __future__ import annotations

import os

from .sourceterm import SourceTerm

ENERGY_MIN_MEV = 1.0e-5
ENERGY_MAX_MEV = 15.0
YIELD_FLOOR = 1.0e-12

# Colour scheme (GraphBank/SKILLS_GraphBank.md).  The photon and BEN blues are
# deliberately far apart in lightness: DEP is a deep saturated blue, BEN the
# lighter royal/cornflower blue, so the two never read as the same series.
COLOURS = {
    "photon":  "#1034A6",   # blue (Egyptian blue)
    "dee":     "#8F00C7",   # violet
    "ben":     "#6E9AF0",   # light royal blue (cornflower)
    "bep":     "#A6193C",   # rouge
    "alpha":   "#7B3F00",   # brown
    "neutron": "#00A550",   # green
    "recoil":  "#000000",   # black
    "ff":      "#5B2C9B",   # royal purple
}

# Discrete (lollipop) classes: attribute, colour key, legend label.
_DISCRETE = [
    ("photon", "photon", r"$Q_{\gamma}(E)$  DEP (binned)"),
    ("dee", "dee", r"$Q_{e}(E)$  DEE (binned)"),
    ("alpha", "alpha", r"$Q_{\alpha}(E)$  DEA"),
    ("neutron", "neutron", r"$Q_{n}(E)$  N"),
    ("recoil", "recoil", r"$Q_{REC}(E)$  DER"),
]

# Output formats written for every graph: true-vector SVG, true-vector PDF, and
# a high-resolution raster PNG.  Nothing on the axes is ever rasterised, so the
# SVG and PDF carry the full data as resolution-independent vector paths; the
# PNG is the only raster.
PNG_DPI = 600


def _apply_house_style(plt) -> None:
    """Set publication-grade rcParams for crisp, professional output.

    The critical settings for sharp vector files are ``svg.fonttype = 'none'``
    (text stays as selectable text rather than being flattened) and leaving
    every artist un-rasterised so SVG/PDF remain pure vector.
    """
    plt.rcParams.update({
        # Keep text as text in the SVG (crisp at any zoom, still editable);
        # PDF embeds real (Type-3-free) fonts.
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # Type-set look.
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "font.size": 12,
        "axes.titlesize": 15,
        "axes.labelsize": 13,
        "legend.fontsize": 11,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        # Clean, crisp frame.
        "axes.linewidth": 1.0,
        "axes.edgecolor": "#222222",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "lines.antialiased": True,
        "patch.antialiased": True,
    })


def plot_source_term(
    st: SourceTerm,
    out_dir: str,
    *,
    stem: str = "IREST",
    title: str = "",
    energy_max: float = ENERGY_MAX_MEV,
) -> list[str]:
    """Render ``st`` to ``<out_dir>/<stem>_ALL.{svg,pdf,png}``.

    Returns the list of files written (empty if plotting was unavailable).  The
    SVG and PDF are **pure vector** (no artist is rasterised); the PNG is
    rasterised at :data:`PNG_DPI`.  Pass a *binned* source term so the photon and
    electron classes are tidy.  Raise ``energy_max`` (e.g. to 200 MeV) to bring
    the fission fragments into view.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:                              # pragma: no cover
        print(f"  ! Plotting skipped (matplotlib unavailable): {exc}")
        return []

    _apply_house_style(plt)

    os.makedirs(out_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11.0, 7.6))
    offscale: list[str] = []

    def note_offscale(label: str, energies, yields, total: float) -> None:
        """Record a class that carries yield but cannot appear in the window."""
        visible = any(y > YIELD_FLOOR and ENERGY_MIN_MEV <= e <= energy_max
                      for e, y in zip(energies, yields))
        if total > 0.0 and not visible:
            e_lo = min(energies) if energies else 0.0
            e_hi = max(energies) if energies else 0.0
            reason = ("below the yield floor" if e_lo >= ENERGY_MIN_MEV and e_hi <= energy_max
                      else f"at {e_lo:.3G}-{e_hi:.3G} MeV")
            offscale.append(f"{label}:  {total:.3E} /decay, {reason}")

    def lollipops(energies, yields, colour, label):
        pts = [(e, y) for e, y in zip(energies, yields) if y > 0.0]
        if not pts:
            return
        es, ys = zip(*pts)
        # True vector stems + markers (never rasterised) so SVG/PDF stay sharp
        # at any zoom.  The stems are faint: with ~1300 lines they would
        # otherwise merge into a solid curtain and bury the continuous spectra.
        ax.vlines(es, YIELD_FLOOR, ys, colors=colour, linewidth=0.6, alpha=0.30)
        ax.scatter(es, ys, s=14, color=colour, label=label, edgecolors="none", zorder=3)

    # Discrete lollipop classes.
    for attr, ckey, label in _DISCRETE:
        d = getattr(st, attr)
        lollipops(d.energy, d.yield_, COLOURS[ckey], label)
        note_offscale(label.split("  ")[-1], d.energy, d.yield_, d.total)

    # Fission fragments (heavy + light): plotted even though they sit far above
    # the default energy window -- raise energy_max to see them.
    ff_e = list(st.hff.energy) + list(st.lff.energy)
    ff_y = list(st.hff.yield_) + list(st.lff.yield_)
    lollipops(ff_e, ff_y, COLOURS["ff"], r"$Q_{FF}(E)$  FF (LFF + HFF)")
    note_offscale("FF (LFF + HFF)", ff_e, ff_y, st.hff.total + st.lff.total)

    # Continuous beta spectra, drawn last and thick so they read over the stems.
    for attr, ckey, label in (
        ("ben", "ben", r"$Q_{\beta^{-}}(E)$  BEN"),
        ("bep", "bep", r"$Q_{\beta^{+}}(E)$  BEP"),
    ):
        b = getattr(st, attr)
        pts = [(e, y) for e, y in zip(b.energy, b.pdf) if y > 0.0]
        note_offscale(label.split("  ")[-1], b.energy, b.pdf, b.total)
        if not pts:
            continue
        es, ys = zip(*pts)
        ax.plot(es, ys, color=COLOURS[ckey], linewidth=3.2, label=label,
                solid_capstyle="round", zorder=4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(ENERGY_MIN_MEV, energy_max)
    ax.set_ylim(bottom=YIELD_FLOOR)
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Particle emission yield per nuclear transformation")
    ax.set_title(title or "All-particle emission spectra $Q(E)$")
    ax.grid(True, which="major", linewidth=0.8, alpha=0.55)
    ax.grid(True, which="minor", linewidth=0.5, alpha=0.30)
    ax.tick_params(which="both", direction="out", length=4, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    # Legend and off-scale note both live BELOW the axes: with this many lines
    # on the plot there is no in-axes region either could occupy without hiding
    # data, and an annotation that overlaps the spectra defeats its own purpose.
    # Figure-fraction layout, stacked upward: off-scale caption, then legend,
    # then the axes.  Each band is wide enough for the x-axis label to clear the
    # legend below it.
    handles, labels = ax.get_legend_handles_labels()
    legend_y = 0.115 if offscale else 0.030
    bottom = (legend_y + 0.170) if handles else (0.105 if offscale else 0.085)
    fig.subplots_adjust(bottom=bottom, top=0.945, left=0.090, right=0.980)

    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=4, frameon=True,
                   edgecolor="#CCCCCC", fancybox=False, framealpha=1.0,
                   bbox_to_anchor=(0.5, legend_y), borderaxespad=0.0,
                   scatterpoints=1, markerscale=1.6)

    if offscale:
        fig.text(
            0.5, 0.014,
            "Emitted, but outside this window:   " + "   |   ".join(offscale)
            + "\nWiden the energy axis with  --graph-energy-max  to bring them into view.",
            ha="center", va="bottom", fontsize=10, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F7F7F7",
                      edgecolor="#BBBBBB"),
        )

    written: list[str] = []
    base = os.path.join(out_dir, f"{stem}_ALL")
    # SVG/PDF: pure vector, so DPI is irrelevant to their sharpness.
    for ext, dpi in (("svg", None), ("pdf", None), ("png", PNG_DPI)):
        out_path = f"{base}.{ext}"
        fig.savefig(out_path, dpi=dpi)
        written.append(out_path)
    plt.close(fig)
    return written
