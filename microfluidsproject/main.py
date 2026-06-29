"""
main.py — Entry point. Builds the full interactive GUI using matplotlib.

Layout (single Figure, multiple Axes):
  ┌────────────────────────┬──────────────┐
  │                        │  Stats axes  │
  │   Network canvas       │  (bar chart) │
  │                        │              │
  │                        ├──────────────┤
  │                        │  Log panel   │
  ├────────────────────────┴──────────────┤
  │      Button row  (Acid | Buffer | Reset)       │
  └─────────────────────────────────────────────────┘
"""

from __future__ import annotations
import sys
import math
import textwrap
from collections import deque

import numpy as np
import matplotlib
matplotlib.use("TkAgg" if "tkinter" in sys.modules else "Qt5Agg")

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.widgets import Button
from matplotlib.animation import FuncAnimation
from matplotlib.patches import FancyBboxPatch

from physics import Network, ACID_PH, BUFFER_PH
from renderer import NetworkRenderer, BG_COLOR, TEXT_C, DIM_TEXT, ACID_C, BUFFER_C, GEL_IDLE, GEL_SWELL


# ─────────────────────────────────────────────────────────────────────────────
# Application state
# ─────────────────────────────────────────────────────────────────────────────

TREE_DEPTH = 3
network    = Network(depth=TREE_DEPTH)

tick_ref   = [0]   # mutable counter shared with renderer

log_lines: deque[str] = deque(maxlen=8)
injection_counts = {"acid": 0, "buffer": 0}
acid_ph   = ACID_PH
buffer_ph = BUFFER_PH


# ─────────────────────────────────────────────────────────────────────────────
# Figure setup
# ─────────────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 9), facecolor=BG_COLOR)
fig.canvas.manager.set_window_title(
    "Chemo-Fluidic Synapse Degradation-Driven Hardware Plasticity"
)

# GridSpec: main canvas | side panel
#           button row at bottom
outer_gs = gridspec.GridSpec(
    2, 1,
    figure=fig,
    height_ratios=[8.5, 1.0],
    hspace=0.02,
    left=0.01, right=0.99, top=0.97, bottom=0.03
)

inner_gs = gridspec.GridSpecFromSubplotSpec(
    2, 2,
    subplot_spec=outer_gs[0],
    width_ratios=[3.0, 1.0],
    height_ratios=[1.6, 1.0],
    wspace=0.02, hspace=0.04
)

ax_net   = fig.add_subplot(inner_gs[:, 0])   # network canvas (both rows, left col)
ax_bar   = fig.add_subplot(inner_gs[0, 1])   # bar chart (top-right)
ax_log   = fig.add_subplot(inner_gs[1, 1])   # log (bottom-right)
ax_btn   = fig.add_subplot(outer_gs[1])      # button strip

for ax in (ax_net, ax_bar, ax_log, ax_btn):
    ax.set_facecolor(BG_COLOR)
    for spine in ax.spines.values():
        spine.set_color("#1a2a3a")

# Network canvas
ax_net.set_aspect("equal")
ax_net.set_xlim(0, 1)
ax_net.set_ylim(0, 1)
ax_net.axis("off")

# Bar chart
ax_bar.set_facecolor("#060c14")
ax_bar.tick_params(colors=DIM_TEXT, labelsize=6)
for s in ax_bar.spines.values():
    s.set_color("#0f1e2e")
ax_bar.set_title("Channel State", color=TEXT_C, fontsize=7, fontfamily="monospace", pad=3)

# Log
ax_log.axis("off")
ax_log.set_facecolor("#060c14")
_log_text = ax_log.text(
    0.02, 0.95, "",
    transform=ax_log.transAxes,
    fontsize=6.5, color="#5a9a5a",
    va="top", fontfamily="monospace",
    linespacing=1.5,
)

# Compute layout (0-1 normalised coords)
network.compute_layout(
    canvas_w=1.0, canvas_h=1.0,
    margin_x=0.09, margin_top=0.08, margin_bottom=0.05
)

renderer = NetworkRenderer(network, ax_net, tick_ref)


# ─────────────────────────────────────────────────────────────────────────────
# Legend
# ─────────────────────────────────────────────────────────────────────────────

legend_elements = [
    mpatches.Patch(color=tuple(ACID_C),    label="Acid  (LTP — erosion)"),
    mpatches.Patch(color=tuple(BUFFER_C),  label="Buffer (LTD — gel swell)"),
    mpatches.Patch(color=tuple(GEL_IDLE),  label="Hydrogel (idle)"),
    mpatches.Patch(color=tuple(GEL_SWELL), label="Hydrogel (swollen)"),
    mpatches.Patch(color="#00c8ff",        label="Eroded channel"),
    mpatches.Patch(color="#1a3a5c",        label="Base channel"),
]
ax_net.legend(
    handles=legend_elements,
    loc="lower left",
    fontsize=5.5,
    framealpha=0.55,
    facecolor="#0a0e1a",
    edgecolor="#1e3a5a",
    labelcolor=TEXT_C,
    handlelength=1.0,
)


# ─────────────────────────────────────────────────────────────────────────────
# Buttons
# ─────────────────────────────────────────────────────────────────────────────

ax_btn.axis("off")

btn_ax_acid   = fig.add_axes([0.04, 0.01, 0.18, 0.07])
btn_ax_buf    = fig.add_axes([0.25, 0.01, 0.18, 0.07])
btn_ax_reset  = fig.add_axes([0.46, 0.01, 0.12, 0.07])
btn_ax_acid2  = fig.add_axes([0.61, 0.01, 0.18, 0.07])  # alt acid pH
btn_ax_buf2   = fig.add_axes([0.82, 0.01, 0.14, 0.07])  # alt buffer pH

for bax in (btn_ax_acid, btn_ax_buf, btn_ax_reset, btn_ax_acid2, btn_ax_buf2):
    bax.set_facecolor(BG_COLOR)

btn_acid  = Button(btn_ax_acid,  "[A] INJECT ACID\n(LTP — Erosion)",
                   color="#2a0808", hovercolor="#661111")
btn_buf   = Button(btn_ax_buf,   "[B] INJECT BUFFER\n(LTD — Gel Swell)",
                   color="#081a10", hovercolor="#114422")
btn_reset = Button(btn_ax_reset, "[R] RESET\nNetwork",
                   color="#0d1520", hovercolor="#1a2a3a")
btn_acid2 = Button(btn_ax_acid2, "[A] STRONG ACID pH=2\n(Max Erosion)",
                   color="#3a0808", hovercolor="#881111")
btn_buf2  = Button(btn_ax_buf2,  "[B] STRONG BUFFER\npH=12",
                   color="#0a1f10", hovercolor="#1a4422")

for btn, fc, ec in [
    (btn_acid,  "#ff4444", "#cc2222"),
    (btn_buf,   "#44ff88", "#22aa55"),
    (btn_reset, "#3a8fd4", "#1e5a8a"),
    (btn_acid2, "#ff7755", "#dd4422"),
    (btn_buf2,  "#55ffaa", "#33cc77"),
]:
    btn.label.set_color(fc)
    btn.label.set_fontsize(7.5)
    btn.label.set_fontfamily("monospace")


# ─────────────────────────────────────────────────────────────────────────────
# Stats bar chart
# ─────────────────────────────────────────────────────────────────────────────

_bar_artists = []

def update_bar_chart():
    global _bar_artists
    for art in _bar_artists:
        try:
            art.remove()
        except Exception:
            pass
    _bar_artists.clear()

    rows = network.summary()
    if not rows:
        return

    n = len(rows)
    y_pos = np.arange(n)
    max_r = max(r["resistance"] for r in rows) or 1
    max_f = max(r["flow_rate"] for r in rows) or 1
    max_w = max(r["width"] for r in rows) or 1

    ax_bar.set_ylim(-0.5, n - 0.5)
    ax_bar.set_xlim(-0.05, 1.05)
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels([r["id"][:5] for r in rows], fontsize=5,
                            color=DIM_TEXT, fontfamily="monospace")
    ax_bar.set_xticks([])

    for i, row in enumerate(rows):
        # Resistance (normalised, blue)
        rw = row["resistance"] / max_r * 0.45
        b1 = ax_bar.barh(i + 0.12, rw, height=0.22, left=0,
                          color="#1e5a9f", alpha=0.85)[0]
        _bar_artists.append(b1)

        # Flow (normalised, cyan)
        fw = row["flow_rate"] / max_f * 0.45
        b2 = ax_bar.barh(i - 0.12, fw, height=0.22, left=0,
                          color="#00a0c8", alpha=0.85)[0]
        _bar_artists.append(b2)

        # Width (normalised, green shading on right side)
        wfrac = (row["width"] - 1.0) / 11.0
        b3 = ax_bar.barh(i, wfrac * 0.45, height=0.05, left=0.50,
                          color=tuple(np.clip([0, wfrac, 1 - wfrac * 0.5], 0, 1)),
                          alpha=0.7)[0]
        _bar_artists.append(b3)

        # Gel block indicator
        gfrac = row["gel_block%"] / 100.0
        if gfrac > 0.01:
            b4 = ax_bar.barh(i, gfrac * 0.45, height=0.05, left=0.50,
                              color=(0.7, 0.3, 1.0), alpha=0.8)[0]
            _bar_artists.append(b4)

    # Labels (static, only set once — we'll just redraw)
    t1 = ax_bar.text(0.225, -0.55, "← Resistance  /  Flow →",
                     fontsize=5, color=DIM_TEXT, ha="center",
                     transform=ax_bar.get_xaxis_transform(),
                     fontfamily="monospace")
    t2 = ax_bar.text(0.725, -0.55, "← Width  |  Gel% →",
                     fontsize=5, color=DIM_TEXT, ha="center",
                     transform=ax_bar.get_xaxis_transform(),
                     fontfamily="monospace")
    _bar_artists += [t1, t2]


# ─────────────────────────────────────────────────────────────────────────────
# Event handlers
# ─────────────────────────────────────────────────────────────────────────────

def _inject(ph: float, label: str):
    network.inject_fluid(ph)
    kind = "ACID" if ph < 7 else "BUFFER"
    injection_counts["acid" if ph < 7 else "buffer"] += 1
    msg = (f"[A:{injection_counts['acid']:03d}|B:{injection_counts['buffer']:03d}] "
           f"{kind} pH={ph:.1f} → {'LTP erosion' if ph < 7 else 'LTD gel-swell'}")
    log_lines.append(msg)
    _log_text.set_text("\n".join(log_lines))

def on_acid(event):
    _inject(ACID_PH, "Acid")

def on_buffer(event):
    _inject(BUFFER_PH, "Buffer")

def on_acid_strong(event):
    _inject(2.0, "Strong Acid")

def on_buffer_strong(event):
    _inject(12.0, "Strong Buffer")

def on_reset(event):
    global network, renderer
    network  = Network(depth=TREE_DEPTH)
    network.compute_layout(
        canvas_w=1.0, canvas_h=1.0,
        margin_x=0.09, margin_top=0.08, margin_bottom=0.05
    )
    renderer = NetworkRenderer(network, ax_net, tick_ref)
    injection_counts["acid"] = injection_counts["buffer"] = 0
    log_lines.clear()
    log_lines.append("Network reset.")
    _log_text.set_text("\n".join(log_lines))

btn_acid.on_clicked(on_acid)
btn_buf.on_clicked(on_buffer)
btn_reset.on_clicked(on_reset)
btn_acid2.on_clicked(on_acid_strong)
btn_buf2.on_clicked(on_buffer_strong)


# ─────────────────────────────────────────────────────────────────────────────
# Title / header
# ─────────────────────────────────────────────────────────────────────────────

fig.text(
    0.5, 0.995,
    "[*]  CHEMO-FLUIDIC SYNAPSE DEGRADATION-DRIVEN HARDWARE PLASTICITY",
    ha="center", va="top",
    fontsize=9, color="#3a8fd4",
    fontfamily="monospace", fontweight="bold",
)
fig.text(
    0.5, 0.978,
    "LTP (Acid → Erosion → ↑ Conductance)   |   LTD (Buffer → Gel Swell → ↑ Resistance)",
    ha="center", va="top",
    fontsize=7, color=DIM_TEXT,
    fontfamily="monospace",
)


# ─────────────────────────────────────────────────────────────────────────────
# Animation loop
# ─────────────────────────────────────────────────────────────────────────────

_frame_counter = [0]
_stats_interval = 8   # update bar chart every N frames

def animate(frame: int):
    tick_ref[0] = frame
    _frame_counter[0] += 1

    network.tick(dt=0.04)
    renderer.draw_frame()

    if _frame_counter[0] % _stats_interval == 0:
        update_bar_chart()

    return []


anim = FuncAnimation(
    fig,
    animate,
    interval=22,       # ~45 fps
    blit=False,
    cache_frame_data=False,
)


# ─────────────────────────────────────────────────────────────────────────────
# Initial solve + draw
# ─────────────────────────────────────────────────────────────────────────────

network.solve_flow()
renderer.draw_frame()
update_bar_chart()
log_lines.append("System ready. Inject fluid to begin.")
_log_text.set_text("\n".join(log_lines))

plt.show()
