"""Shared experiment I/O and plotting helpers."""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..metrics import pairwise_orders, asymptotic_fit

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 180,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "lines.marker": "o",
    "lines.markersize": 4,
})


def save_table(frame, filename):
    path = RESULTS / filename
    frame.to_csv(path, index=False, float_format="%.12e")
    return path


def save_figure(fig, filename):
    path = FIGURES / filename
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def add_pairwise_orders(frame, group_columns, metric_columns):
    frame = frame.sort_values(group_columns + ["N"]).copy()
    for metric in metric_columns:
        frame[metric + "_order"] = np.nan
        for _, indices in frame.groupby(group_columns).groups.items():
            ordered = frame.loc[indices].sort_values("N")
            frame.loc[ordered.index, metric + "_order"] = pairwise_orders(
                ordered["N"].to_numpy(), ordered[metric].to_numpy()
            )
    return frame


def fit_summary(frame, group_columns, metric_columns, last=3):
    rows = []
    grouped = [((), frame)] if not group_columns else frame.groupby(group_columns)
    for key, group in grouped:
        if not isinstance(key, tuple):
            key = (key,)
        group = group.sort_values("N")
        for metric in metric_columns:
            fit = asymptotic_fit(group["N"].to_numpy(), group[metric].to_numpy(), last=last)
            row = dict(zip(group_columns, key))
            row.update({"metric": metric, **fit})
            rows.append(row)
    return pd.DataFrame(rows)


def fit_label(group, metric, label, x_column="N", last=3):
    """Compact curve label with fitted order and maximum natural-log deviation."""
    ordered = group.sort_values(x_column)
    fit = asymptotic_fit(ordered[x_column].to_numpy(),
                         ordered[metric].to_numpy(), last=last)
    return (f"{label}: $p={fit['order']:.3f}$, "
            f"$\\delta_{{\\log}}={fit['max_log_deviation']:.2e}$")


def loglog_error_plot(frame, group_column, metric_columns, title, filename,
                      x_column="N", reference_orders=()):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for label, group in frame.groupby(group_column):
        group = group.sort_values(x_column)
        for metric in metric_columns:
            ax.loglog(group[x_column], group[metric], label=f"{label}: {metric}")
    x = np.array(sorted(frame[x_column].unique()), dtype=float)
    if len(x) >= 2:
        anchor = frame[metric_columns[0]].max()
        for order in reference_orders:
            y = anchor * (x / x.max()) ** (-order)
            ax.loglog(x, y, linestyle="--", marker="none", label=f"reference N^-{order}")
    ax.set_xlabel(x_column)
    ax.set_ylabel("error")
    ax.set_title(title)
    ax.legend(fontsize=8)
    return save_figure(fig, filename)
