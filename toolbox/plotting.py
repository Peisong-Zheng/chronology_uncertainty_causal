"""Panel labels and export helpers; each analysis retains its original style."""

from pathlib import Path
import string

import matplotlib as mpl
import numpy as np


def add_panel_labels(axes, x=-0.10, y=1.04):
    """Label only supplied scientific axes, in their explicit reading order."""
    panels = np.asarray(axes, dtype=object).ravel().tolist()
    for ax in panels:
        for label in list(ax.texts):
            if label.get_gid() == "panel-label":
                label.remove()
    if len(panels) <= 1:
        return []
    if len(panels) > len(string.ascii_lowercase):
        raise ValueError("Specify a panel-label scheme for more than 26 panels.")
    labels = []
    for letter, ax in zip(string.ascii_lowercase, panels):
        label = ax.text(x, y, letter, transform=ax.transAxes, fontsize=14,
                        fontweight="bold", fontstyle="normal", ha="left", va="bottom",
                        color="#111111", clip_on=False, gid="panel-label")
        labels.append(label)
    return labels


def save_figure(fig, output_dir, name, axes=None):
    """Save one intact research figure as a PNG preview and vector PDF."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if axes is not None:
        add_panel_labels(axes)
    # Keep the chosen research canvas; do not crop to an accidental text extent.
    fig.savefig(output_dir / f"{name}.png", dpi=180)
    with mpl.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(output_dir / f"{name}.pdf")
