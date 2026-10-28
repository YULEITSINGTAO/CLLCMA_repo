"""Small figure helpers used by the stage and inferCNV plot scripts.

Okabe–Ito colors follow Wong, Nature Methods 8, 441 (2011).
`export_figure` writes the requested file types and ignores extra keyword
arguments that the original plotting calls pass for provenance manifests.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}


@contextmanager
def style_context(*_args, **_kwargs):
    yield


def export_figure(fig, path, formats=("pdf", "png"), dpi=300, bbox_inches="tight", **_kwargs):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(out.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches=bbox_inches)
    return out
