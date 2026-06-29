"""
physics.py — Core physics model for the chemo-fluidic network.

Classes
-------
Channel   : A fluidic edge; stores width, resistance, erosion state.
Hydrogel  : Valve placed at a branch entrance; swells under buffer exposure.
Node      : A bifurcation point in the network graph.
Network   : Directed graph of nodes/channels; computes flow via Kirchhoff's laws.
FluidPacket: A discrete parcel of fluid travelling through the network.
"""

from __future__ import annotations
import math
import uuid
from dataclasses import dataclass, field
from typing import Optional
import numpy as np


# ─────────────────────────── constants ────────────────────────────────────────

# Hagen-Poiseuille: R ∝ L / r^4  (we use width as proxy for radius)
BASE_LENGTH         = 1.0        # normalised channel length
BASE_WIDTH          = 4.0        # µm, initial channel half-width
MIN_WIDTH           = 1.0        # µm, fully blocked
MAX_WIDTH           = 12.0       # µm, fully eroded
EROSION_RATE        = 0.35       # width gained per acid injection
SWELL_RATE          = 0.40       # gel expansion per buffer injection
RECOVERY_RATE       = 0.04       # passive recovery per timestep (both erosion & swell)
HYDROGEL_BASE_BLOCK = 0.0        # gel starts fully open
HYDROGEL_MAX_BLOCK  = 3.5        # maximum gel expansion (µm equivalent)
VISCOSITY           = 1.0        # normalised dynamic viscosity
SOURCE_PRESSURE     = 100.0      # Pa applied at root inlet

ACID_PH   = 4.0
BUFFER_PH = 8.5


# ─────────────────────────── data model ───────────────────────────────────────

@dataclass
class Hydrogel:
    """
    Valve at the entrance of a branch channel.

    swelling : float  — current expansion in µm-equivalent (0 = open, MAX = closed)
    """
    channel_id: str
    swelling: float = 0.0

    # ── derived ──────────────────────────────────────────────────────────────
    @property
    def block_fraction(self) -> float:
        """0.0 = fully open, 1.0 = fully occluded."""
        return min(self.swelling / HYDROGEL_MAX_BLOCK, 1.0)

    @property
    def effective_width_reduction(self) -> float:
        return self.swelling

    # ── dynamics ─────────────────────────────────────────────────────────────
    def expose_to_buffer(self, ph: float) -> None:
        """Buffer (pH > 7) causes swelling."""
        strength = (ph - 7.0) / 7.0          # 0..1
        self.swelling = min(self.swelling + SWELL_RATE * strength * 2,
                            HYDROGEL_MAX_BLOCK)

    def passive_recovery(self) -> None:
        """Hydrogel slowly shrinks back when not exposed."""
        self.swelling = max(0.0, self.swelling - RECOVERY_RATE * 0.5)


@dataclass
class Channel:
    """
    A directed fluidic channel (edge in the network graph).

    width       : current lumen half-width (µm)
    hydrogel    : optional Hydrogel valve at the *entrance* of this channel
    flow_rate   : most recent computed volumetric flow (normalised)
    fluid_ph    : pH of fluid currently in channel (None = empty)
    """
    id: str
    parent_node_id: str
    child_node_id: str
    length: float = BASE_LENGTH
    width: float = BASE_WIDTH
    hydrogel: Optional[Hydrogel] = None
    flow_rate: float = 0.0
    fluid_ph: Optional[float] = None
    fluid_progress: float = 0.0   # 0..1 animation position

    # ── resistance (Hagen-Poiseuille) ────────────────────────────────────────
    @property
    def effective_width(self) -> float:
        gel_reduction = self.hydrogel.effective_width_reduction if self.hydrogel else 0.0
        return max(MIN_WIDTH, self.width - gel_reduction)

    @property
    def resistance(self) -> float:
        r = self.effective_width
        return (8 * VISCOSITY * self.length) / (math.pi * r ** 4)

    # ── acid erosion (LTP) ───────────────────────────────────────────────────
    def expose_to_acid(self, ph: float) -> None:
        strength = max(0, (7.0 - ph) / 7.0)   # stronger at lower pH
        self.width = min(self.width + EROSION_RATE * strength * 2, MAX_WIDTH)

    # ── passive recovery ─────────────────────────────────────────────────────
    def passive_recovery(self) -> None:
        self.width = max(BASE_WIDTH, self.width - RECOVERY_RATE * 0.3)
        if self.hydrogel:
            self.hydrogel.passive_recovery()


@dataclass
class Node:
    """A bifurcation (or root / leaf) point in the network."""
    id: str
    depth: int = 0
    pressure: float = 0.0
    is_root: bool = False
    is_leaf: bool = False
    # pixel positions assigned by layout engine
    x: float = 0.0
    y: float = 0.0


@dataclass
class FluidPacket:
    """A discrete parcel of fluid travelling through a specific channel."""
    channel_id: str
    ph: float
    progress: float = 0.0   # 0 = just entered, 1 = exited


# ─────────────────────────── network ──────────────────────────────────────────

class Network:
    """
    Directed root-like (binary-tree) fluidic network.

    Parameters
    ----------
    depth : int — number of bifurcation levels (depth=3 → 8 leaves)
    """

    def __init__(self, depth: int = 3):
        self.depth = depth
        self.nodes:    dict[str, Node]    = {}
        self.channels: dict[str, Channel] = {}
        self.hydrogels: dict[str, Hydrogel] = {}

        # adjacency: node_id → list[channel_id] (outgoing)
        self.children:  dict[str, list[str]] = {}
        # reverse: channel_id → parent node_id
        self.channel_parent: dict[str, str] = {}

        self.fluid_packets: list[FluidPacket] = []
        self.root_id: str = ""
        self.leaf_ids: list[str] = []

        self._build_tree(depth)

    # ── construction ─────────────────────────────────────────────────────────

    def _make_node(self, depth: int, is_root=False, is_leaf=False) -> Node:
        nid = str(uuid.uuid4())[:8]
        n = Node(id=nid, depth=depth, is_root=is_root, is_leaf=is_leaf)
        self.nodes[nid] = n
        self.children[nid] = []
        return n

    def _make_channel(self, parent: Node, child: Node,
                      has_hydrogel: bool = False) -> Channel:
        cid = str(uuid.uuid4())[:8]
        gel = None
        if has_hydrogel:
            gel = Hydrogel(channel_id=cid)
            self.hydrogels[cid] = gel
        ch = Channel(id=cid, parent_node_id=parent.id,
                     child_node_id=child.id, hydrogel=gel)
        self.channels[cid] = ch
        self.children[parent.id].append(cid)
        self.channel_parent[cid] = parent.id
        return ch

    def _build_tree(self, depth: int) -> None:
        """Recursively build a binary-tree fluidic network."""
        root = self._make_node(0, is_root=True)
        self.root_id = root.id

        def recurse(parent_node: Node, current_depth: int):
            if current_depth >= depth:
                parent_node.is_leaf = True
                self.leaf_ids.append(parent_node.id)
                return
            # Two child branches — BOTH get a hydrogel at entrance
            for _ in range(2):
                child = self._make_node(current_depth + 1,
                                        is_leaf=(current_depth + 1 >= depth))
                ch = self._make_channel(parent_node, child,
                                        has_hydrogel=True)
                if child.is_leaf:
                    self.leaf_ids.append(child.id)
                else:
                    recurse(child, current_depth + 1)

        recurse(root, 0)

    # ── layout (assign x/y pixel coords) ─────────────────────────────────────

    def compute_layout(self, canvas_w: float, canvas_h: float,
                       margin_x: float = 60, margin_top: float = 40,
                       margin_bottom: float = 40) -> None:
        """
        Assign (x, y) to every node for rendering.
        Root at top-centre; leaves spread across the bottom.
        """
        usable_w = canvas_w - 2 * margin_x
        usable_h = canvas_h - margin_top - margin_bottom
        level_height = usable_h / self.depth if self.depth else usable_h

        # BFS to assign positions
        from collections import deque

        # First pass: count leaves per subtree to space x-positions
        leaf_count: dict[str, int] = {}

        def count_leaves(node_id: str) -> int:
            ch_ids = self.children[node_id]
            if not ch_ids:
                leaf_count[node_id] = 1
                return 1
            total = 0
            for cid in ch_ids:
                ch = self.channels[cid]
                total += count_leaves(ch.child_node_id)
            leaf_count[node_id] = total
            return total

        count_leaves(self.root_id)
        total_leaves = leaf_count[self.root_id]

        # Second pass: assign x by leaf-weighted position
        x_cursor: list[float] = [margin_x]

        def assign_pos(node_id: str, depth: int):
            node = self.nodes[node_id]
            ch_ids = self.children[node_id]
            node.y = margin_top + depth * level_height

            if not ch_ids:
                node.x = x_cursor[0] + (usable_w / total_leaves) * 0.5
                x_cursor[0] += usable_w / total_leaves
                return

            # Position node at midpoint of its children
            left_x_before = x_cursor[0]
            for cid in ch_ids:
                assign_pos(self.channels[cid].child_node_id, depth + 1)
            right_x_after = x_cursor[0]
            node.x = (left_x_before + right_x_after) / 2

        assign_pos(self.root_id, 0)

    # ── flow solver (Kirchhoff / nodal analysis) ──────────────────────────────

    def solve_flow(self, inlet_pressure: float = SOURCE_PRESSURE) -> None:
        """
        Solve for nodal pressures using conductance matrix (Kirchhoff).
        Assigns self.channels[*].flow_rate.
        """
        node_list = list(self.nodes.keys())
        n = len(node_list)
        idx = {nid: i for i, nid in enumerate(node_list)}

        G = np.zeros((n, n))
        b = np.zeros(n)

        for ch in self.channels.values():
            p = idx[ch.parent_node_id]
            c = idx[ch.child_node_id]
            g = 1.0 / ch.resistance
            G[p][p] += g
            G[c][c] += g
            G[p][c] -= g
            G[c][p] -= g

        # Boundary: root node = inlet pressure
        root_i = idx[self.root_id]
        G[root_i, :] = 0
        G[root_i, root_i] = 1
        b[root_i] = inlet_pressure

        # Boundary: leaf nodes = 0 (atmospheric)
        for lid in self.leaf_ids:
            li = idx[lid]
            G[li, :] = 0
            G[li, li] = 1
            b[li] = 0.0

        pressures = np.linalg.solve(G, b)

        for nid, node in self.nodes.items():
            node.pressure = pressures[idx[nid]]

        for ch in self.channels.values():
            dP = (self.nodes[ch.parent_node_id].pressure -
                  self.nodes[ch.child_node_id].pressure)
            ch.flow_rate = max(0.0, dP / ch.resistance)

    # ── injection ─────────────────────────────────────────────────────────────

    def inject_fluid(self, ph: float) -> None:
        """
        Inject a fluid pulse from the root.
        Spawns FluidPackets in every channel, weighted by flow_rate.
        Applies physics immediately (simplified: instant propagation per step).
        """
        self.solve_flow()
        self.fluid_packets.clear()

        # Propagate pH through the tree and apply physics
        self._propagate(self.root_id, ph)

    def _propagate(self, node_id: str, ph: float) -> None:
        ch_ids = self.children[node_id]
        if not ch_ids:
            return
        for cid in ch_ids:
            ch = self.channels[cid]
            # Spawn visual packet
            pkt = FluidPacket(channel_id=cid, ph=ph)
            self.fluid_packets.append(pkt)
            ch.fluid_ph = ph
            ch.fluid_progress = 0.0
            # Apply physics
            if ph < 7.0:           # Acid → LTP (erosion)
                ch.expose_to_acid(ph)
            else:                  # Buffer → LTD (hydrogel swell)
                if ch.hydrogel:
                    ch.hydrogel.expose_to_buffer(ph)
            # Recurse (fluid flows if not fully blocked)
            if ch.hydrogel is None or ch.hydrogel.block_fraction < 0.95:
                self._propagate(ch.child_node_id, ph)

    # ── animation step ────────────────────────────────────────────────────────

    def tick(self, dt: float = 0.05) -> None:
        """Advance fluid packet animations and apply passive recovery."""
        for pkt in self.fluid_packets:
            ch = self.channels[pkt.channel_id]
            speed = max(0.1, ch.flow_rate / (SOURCE_PRESSURE * 0.1)) * dt * 2
            pkt.progress = min(1.0, pkt.progress + speed)
            ch.fluid_progress = pkt.progress

        # Remove completed packets
        self.fluid_packets = [p for p in self.fluid_packets if p.progress < 1.0]

        # Passive recovery
        for ch in self.channels.values():
            ch.passive_recovery()
            if ch.fluid_progress >= 1.0:
                ch.fluid_ph = None

    # ── diagnostics ──────────────────────────────────────────────────────────

    def summary(self) -> list[dict]:
        rows = []
        for ch in self.channels.values():
            rows.append({
                "id":          ch.id,
                "width":       round(ch.width, 2),
                "eff_width":   round(ch.effective_width, 2),
                "resistance":  round(ch.resistance, 4),
                "flow_rate":   round(ch.flow_rate, 4),
                "gel_swell":   round(ch.hydrogel.swelling, 2) if ch.hydrogel else 0,
                "gel_block%":  round((ch.hydrogel.block_fraction * 100), 1) if ch.hydrogel else 0,
            })
        return rows
