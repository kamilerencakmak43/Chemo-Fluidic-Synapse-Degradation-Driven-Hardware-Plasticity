# Chemo-Fluidic Synapse Degradation-Driven Hardware Plasticity

## Requirements
```
pip install matplotlib numpy
```
Optional for a nicer GUI window backend:
```
pip install PyQt5   # or PyQt6 / tk
```

## Run
```
python main.py
```

## Controls
| Button | Action |
|--------|--------|
| [A] INJECT ACID (LTP) | Injects pH=4 acid → channels erode (widen), resistance drops |
| [B] INJECT BUFFER (LTD) | Injects pH=8.5 buffer → hydrogels swell, flow restricted |
| [A] STRONG ACID pH=2 | Maximum erosion pulse |
| [B] STRONG BUFFER pH=12 | Maximum hydrogel swelling |
| [R] RESET | Rebuild fresh network |

## Architecture
```
physics.py   — Network, Channel, Hydrogel, Node, FluidPacket dataclasses + Kirchhoff solver
renderer.py  — Matplotlib draw engine (channels, fluid packets, hydrogel blobs, nodes)
main.py      — Figure layout, animation loop, button handlers, bar chart, event log
```

## Physics Summary
- **Resistance** follows Hagen-Poiseuille: R ∝ L / r⁴
- **LTP (acid)**: erosion widens effective radius → lower R → higher flow
- **LTD (buffer)**: hydrogel swells → reduces effective width → higher R → lower flow
- **Kirchhoff nodal analysis** solves pressure at every bifurcation each frame
- **Passive recovery**: channels and hydrogels slowly revert toward baseline over time
