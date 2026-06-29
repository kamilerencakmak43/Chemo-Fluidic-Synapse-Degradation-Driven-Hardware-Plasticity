"""
renderer.py — Matplotlib rendering engine for the chemo-fluidic network.

Draws every frame into a matplotlib Figure with:
  • Dark microfluidic-chip aesthetic
  • Channels coloured/widened by erosion state
  • Hydrogel blobs at each branch entrance, scaling with swelling
  • Animated fluid packets (acid = red, buffer = green)
  • Pressure labels on nodes
  • Live legend
"""

from __future__ import annotations
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text

from physics import (
    Network, Channel, Hydrogel, Node,
    BASE_WIDTH, MAX_WIDTH, MIN_WIDTH
)

# ── colour constants (normalised 0-1 RGB tuples) ──────────────────────────────
BG_COLOR    = "#0a0e1a"
GRID_COLOR  = "#111928"

CH_BASE     = np.array([0.05, 0.18, 0.35])  # dark steel blue
CH_ERODED   = np.array([0.00, 0.75, 1.00])  # bright cyan
CH_WALL_A   = 0.35

ACID_C      = np.array([1.00, 0.22, 0.22])
BUFFER_C    = np.array([0.20, 1.00, 0.45])

GEL_IDLE    = np.array([0.25, 0.75, 0.60])
GEL_SWELL   = np.array([0.80, 0.35, 1.00])

NODE_C      = "#1e3a5f"
NODE_EDGE   = "#3a8fd4"
ROOT_C      = "#ffd700"
LEAF_C      = "#1a4a4a"

TEXT_C      = "#8ab8d8"
DIM_TEXT    = "#3a5a70"

# channel lumen → rendered linewidth mapping
LW_BASE = 2.0
LW_MAX  = 14.0


def lerp(a, b, t):
    t = float(np.clip(t, 0, 1))
    return a + (b - a) * t


def ch_lw(ch: Channel) -> float:
    span = MAX_WIDTH - MIN_WIDTH
    frac = (ch.effective_width - MIN_WIDTH) / span
    return LW_BASE + frac * (LW_MAX - LW_BASE)


def erosion_color(ch: Channel):
    frac = (ch.width - BASE_WIDTH) / (MAX_WIDTH - BASE_WIDTH)
    frac = float(np.clip(frac, 0, 1))
    rgb = lerp(CH_BASE, CH_ERODED, frac)
    return tuple(rgb)


def fluid_color(ph: float):
    if ph < 7.0:
        return tuple(ACID_C), tuple(lerp(ACID_C, np.ones(3), 0.3))
    return tuple(BUFFER_C), tuple(lerp(BUFFER_C, np.ones(3), 0.3))


class NetworkRenderer:
    """
    Owns a matplotlib Figure and axes.
    Call `draw_frame()` to fully redraw the network.
    """

    def __init__(self, network: Network, ax: plt.Axes, tick_ref: list[int]):
        self.network   = network
        self.ax        = ax
        self.tick_ref  = tick_ref   # mutable reference so caller can update
        self._artists: list = []

    # ── public ───────────────────────────────────────────────────────────────

    def draw_frame(self):
        ax = self.ax
        for art in self._artists:
            try:
                art.remove()
            except Exception:
                pass
        self._artists.clear()

        self._draw_channels()
        self._draw_fluid_packets()
        self._draw_hydrogels()
        self._draw_nodes()

    # ── internals ────────────────────────────────────────────────────────────

    def _add(self, art):
        self._artists.append(art)
        return art

    def _draw_channels(self):
        ax = self.ax
        net = self.network
        for ch in net.channels.values():
            pn = net.nodes[ch.parent_node_id]
            cn = net.nodes[ch.child_node_id]
            color = erosion_color(ch)
            lw = ch_lw(ch)

            # Shadow / glow
            glow_lw = lw + 5
            glow_c = color + (0.12,)
            ln_glow = ax.plot(
                [pn.x, cn.x], [pn.y, cn.y],
                color=glow_c, lw=glow_lw,
                solid_capstyle="round", zorder=1
            )[0]
            self._add(ln_glow)

            # Wall
            wall_c = CH_BASE * 0.6
            ln_wall = ax.plot(
                [pn.x, cn.x], [pn.y, cn.y],
                color=tuple(wall_c) + (CH_WALL_A,),
                lw=lw + 2, solid_capstyle="round", zorder=2
            )[0]
            self._add(ln_wall)

            # Core
            ln = ax.plot(
                [pn.x, cn.x], [pn.y, cn.y],
                color=color, lw=lw,
                solid_capstyle="round", zorder=3
            )[0]
            self._add(ln)

            # Width annotation (µm) on mid-point
            mx, my = (pn.x + cn.x) / 2, (pn.y + cn.y) / 2
            txt = ax.text(
                mx, my, f"{ch.effective_width:.1f}",
                fontsize=5.5, color=TEXT_C, ha="center", va="center",
                zorder=6, fontfamily="monospace",
            )
            txt.set_path_effects([
                pe.withStroke(linewidth=2, foreground=BG_COLOR)
            ])
            self._add(txt)

    def _draw_fluid_packets(self):
        ax = self.ax
        net = self.network
        tick = self.tick_ref[0]

        for pkt in net.fluid_packets:
            ch = net.channels[pkt.channel_id]
            pn = net.nodes[ch.parent_node_id]
            cn = net.nodes[ch.child_node_id]

            t = pkt.progress
            fx = pn.x + (cn.x - pn.x) * t
            fy = pn.y + (cn.y - pn.y) * t

            fill, bright = fluid_color(pkt.ph)
            rw = ch_lw(ch)

            # Trail
            t0 = max(0.0, t - 0.20)
            tx0 = pn.x + (cn.x - pn.x) * t0
            ty0 = pn.y + (cn.y - pn.y) * t0
            trail_c = fill + (0.25,)
            tr = ax.plot(
                [tx0, fx], [ty0, fy],
                color=trail_c, lw=rw * 0.7,
                solid_capstyle="round", zorder=7
            )[0]
            self._add(tr)

            # Outer halo
            halo_r = rw * 0.012
            halo_c = fill + (0.18,)
            halo = Circle((fx, fy), halo_r * 2.2, color=halo_c, zorder=8)
            ax.add_patch(halo)
            self._add(halo)

            # Core blob
            blob_c = bright + (0.95,)
            blob = Circle((fx, fy), halo_r, color=blob_c, zorder=9)
            ax.add_patch(blob)
            self._add(blob)

    def _draw_hydrogels(self):
        ax = self.ax
        net = self.network
        tick = self.tick_ref[0]

        for ch_id, gel in net.hydrogels.items():
            ch = net.channels[ch_id]
            pn = net.nodes[ch.parent_node_id]
            cn = net.nodes[ch.child_node_id]

            # Place at 15% along the channel
            t = 0.15
            gx = pn.x + (cn.x - pn.x) * t
            gy = pn.y + (cn.y - pn.y) * t

            frac = gel.block_fraction
            color = tuple(lerp(GEL_IDLE, GEL_SWELL, frac))

            # Base radius grows with swelling
            base_r = 0.012
            r = base_r + frac * 0.035

            # Outer pulsing ring when swollen
            if frac > 0.05:
                pulse = math.sin(tick * 0.15) * 0.004 * frac
                ring_c = color + (0.25 * frac,)
                ring = Circle((gx, gy), r + 0.008 + pulse,
                               color=ring_c, fill=False,
                               linewidth=1.0, zorder=4)
                ax.add_patch(ring)
                self._add(ring)

            # Main hydrogel blob
            alpha = 0.55 + 0.40 * frac
            gel_c = color + (alpha,)
            blob = Circle((gx, gy), r, color=gel_c, zorder=5)
            ax.add_patch(blob)
            self._add(blob)

            # Outline
            edge_c = tuple(lerp(GEL_IDLE, np.ones(3), 0.2)) + (0.6,)
            outline = Circle((gx, gy), r, fill=False, edgecolor=edge_c, linewidth=0.7, zorder=5)
            ax.add_patch(outline)
            self._add(outline)

            # Block percentage label
            if frac > 0.1:
                txt = ax.text(
                    gx, gy - r - 0.012, f"{frac*100:.0f}%",
                    fontsize=5, color=TEXT_C, ha="center", va="top",
                    zorder=10, fontfamily="monospace"
                )
                txt.set_path_effects([
                    pe.withStroke(linewidth=1.5, foreground=BG_COLOR)
                ])
                self._add(txt)

    def _draw_nodes(self):
        ax = self.ax
        net = self.network
        r = 0.018

        for node in net.nodes.values():
            x, y = node.x, node.y

            if node.is_root:
                # Golden burst
                for i in range(8):
                    ang = i * math.pi / 4
                    ray_len = r * 2.0
                    rx, ry = x + math.cos(ang) * ray_len, y + math.sin(ang) * ray_len
                    ray = ax.plot(
                        [x, rx], [y, ry], color="#ffd700", lw=1.5,
                        alpha=0.4, zorder=10
                    )[0]
                    self._add(ray)
                c = Circle((x, y), r * 1.3, color=ROOT_C, zorder=11)
                ax.add_patch(c)
                self._add(c)
                t = ax.text(x, y, "IN", fontsize=6, ha="center", va="center",
                            color="#1a1a00", fontweight="bold",
                            fontfamily="monospace", zorder=12)
                self._add(t)

            elif node.is_leaf:
                c = Circle((x, y), r * 0.7, color=LEAF_C, edgecolor="#2a8888", linewidth=0.8, zorder=10)
                ax.add_patch(c)
                self._add(c)
                t = ax.text(x, y + r * 1.2, "OUT", fontsize=4.5,
                            ha="center", va="bottom", color=DIM_TEXT,
                            fontfamily="monospace", zorder=11)
                self._add(t)

            else:
                c = Circle((x, y), r, color=NODE_C, edgecolor=NODE_EDGE, linewidth=0.8, zorder=10)
                ax.add_patch(c)
                self._add(c)
                # Pressure label
                t = ax.text(
                    x, y - r * 1.5,
                    f"{node.pressure:.0f}Pa",
                    fontsize=4.5, ha="center", va="top",
                    color=DIM_TEXT, fontfamily="monospace", zorder=11
                )
                t.set_path_effects([
                    pe.withStroke(linewidth=1, foreground=BG_COLOR)
                ])
                self._add(t)
