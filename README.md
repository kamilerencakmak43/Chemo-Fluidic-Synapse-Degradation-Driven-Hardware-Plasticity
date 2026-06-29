# Chemo-Fluidic Synaptic Core Simulation

## Overview
This repository contains a physics-based Python simulation of a **Chemo-Fluidic Synaptic Core**. This independent research project explores unconventional computing and neuromorphic hardware by utilizing microfluidic physical limitations as features. 

Instead of traditional silicon-based matrix multiplications, this system models how chemical flow, structural wear (channel degradation), and hydrogel swelling can be used to mimic biological synaptic plasticity (learning and forgetting) directly at the hardware level.

## Core Mechanisms
The simulation visually and mathematically represents a microfluidic channel acting as a synapse. The synaptic weight (W) is determined by the fluidic resistance, which dynamically changes based on chemical interactions:

*   **LTP (Long-Term Potentiation / Learning):** Triggered by an acidic injection. The acid degrades the channel walls, widening the effective radius. This decreases fluidic resistance and increases the synaptic weight permanently.
*   **LTD (Long-Term Depression / Forgetting):** Triggered by a pH buffer injection. The buffer causes the integrated hydrogel layers to swell, narrowing the channel. This increases fluidic resistance and decreases the synaptic weight.
*   **Idle / Neutral Flow:** Represents the resting state where the hydrogel slowly stabilizes without permanent structural changes.

## Mathematical Foundation
The weight adjustment is modeled around the principles of fluid dynamics, heavily inspired by the Hagen-Poiseuille equation. Since fluidic resistance ($R$) is inversely proportional to the fourth power of the radius ($r^4$), even micro-scale physical changes in the channel via degradation or hydrogel swelling result in significant and controllable shifts in synaptic weight.

## Features
*   **Real-time Physical Simulation:** A 2D visual interface built with Tkinter that animates fluid particles, channel degradation, and hydrogel expansion in real-time.
*   **Dynamic Weight Calculation:** Live updating of flow speed and channel width based on user injections.

## How to Run
Requirements
pip install matplotlib numpy
Optional for a nicer GUI window backend:

pip install PyQt5   # or PyQt6 / tk
Run
python main.py
Controls
Button	Action
[A] INJECT ACID (LTP)	Injects pH=4 acid → channels erode (widen), resistance drops
[B] INJECT BUFFER (LTD)	Injects pH=8.5 buffer → hydrogels swell, flow restricted
[A] STRONG ACID pH=2	Maximum erosion pulse
[B] STRONG BUFFER pH=12	Maximum hydrogel swelling
[R] RESET	Rebuild fresh network
Architecture
physics.py   — Network, Channel, Hydrogel, Node, FluidPacket dataclasses + Kirchhoff solver
renderer.py  — Matplotlib draw engine (channels, fluid packets, hydrogel blobs, nodes)
main.py      — Figure layout, animation loop, button handlers, bar chart, event log
Physics Summary
Resistance follows Hagen-Poiseuille: R ∝ L / r⁴
LTP (acid): erosion widens effective radius → lower R → higher flow
LTD (buffer): hydrogel swells → reduces effective width → higher R → lower flow
Kirchhoff nodal analysis solves pressure at every bifurcation each frame
Passive recovery: channels and hydrogels slowly revert toward baseline over time

Author
Kamil Eren Çakmak

Independent Researcher | Unconventional Computing & Microfluidics
