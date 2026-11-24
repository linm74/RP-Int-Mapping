#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Set, Tuple
import matplotlib.patheffects as pe

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# =========================
# USER SETTINGS
# =========================

NETWORK = "g3"                    # "all", "g2", "g3", or "g4"
LAYERS = [2, 3, 4, 5, 9999]       # counts per layer: inner->outer. 9999 = "all remaining"

INNER_RADIUS = 1.4
RING_SPACING_BASE = 1.7
RING_SPACING_GROWTH = 0.75

NEIGHBOR_EXTRA_RING = True        # put non-layer neighbors on an extra outer ring

# Output folder name
OUTPUT_DIR_NAME = "layered_network_outputs"

# =========================


EDGE_A_COL = "#node1"
EDGE_B_COL = "node2"


def load_ranked_genes(rank_csv: Path) -> List[str]:
    df = pd.read_csv(rank_csv)
    if "Gene" not in df.columns:
        raise SystemExit(f"{rank_csv.name}: missing 'Gene' column. Found: {list(df.columns)}")
    return df["Gene"].astype(str).tolist()


def split_into_layers(ranked_genes: List[str], layer_sizes: List[int]) -> List[List[str]]:
    layers: List[List[str]] = []
    idx = 0
    n = len(ranked_genes)

    for size in layer_sizes:
        if idx >= n:
            break
        if size >= 9999:
            layers.append(ranked_genes[idx:])
            idx = n
            break
        layers.append(ranked_genes[idx: idx + size])
        idx += size

    return layers


def build_graph_from_tsv(tsv_path: Path) -> nx.Graph:
    df = pd.read_csv(tsv_path, sep="\t")
    for col in (EDGE_A_COL, EDGE_B_COL):
        if col not in df.columns:
            raise SystemExit(f"{tsv_path.name}: missing '{col}'. Found: {list(df.columns)}")

    G = nx.Graph()
    for a, b in zip(df[EDGE_A_COL].astype(str), df[EDGE_B_COL].astype(str)):
        a = a.strip()
        b = b.strip()
        if not a or not b or a == b:
            continue
        G.add_edge(a, b)
    return G


def get_one_step_neighbors(G: nx.Graph, seeds: Set[str]) -> Set[str]:
    nbrs: Set[str] = set()
    for s in seeds:
        if s in G:
            nbrs.update(G.neighbors(s))
    return nbrs


def circular_positions(nodes: List[str], radius: float, angle_offset: float = 0.0) -> Dict[str, Tuple[float, float]]:
    pos: Dict[str, Tuple[float, float]] = {}
    if not nodes:
        return pos
    for i, n in enumerate(nodes):
        theta = angle_offset + (2.0 * math.pi * i / len(nodes))
        pos[n] = (radius * math.cos(theta), radius * math.sin(theta))
    return pos


def layered_radial_layout(
    layers: List[List[str]],
    extra_ring_nodes: List[str],
    inner_radius: float,
    ring_spacing_base: float,
    ring_spacing_growth: float,
) -> Dict[str, Tuple[float, float]]:
    """
    Places important inner layers close together and pushes outer layers farther away.
    Each ring gets progressively farther from the previous one.
    """
    pos: Dict[str, Tuple[float, float]] = {}
    angle_offset = 0.0

    radius = inner_radius

    for li, layer_nodes in enumerate(layers):
        pos.update(circular_positions(layer_nodes, radius, angle_offset=angle_offset))
        angle_offset += 0.25

        # progressively larger gaps as we move outward
        radius += ring_spacing_base + (li ** 2) * ring_spacing_growth

    if extra_ring_nodes:
        radius += ring_spacing_base + len(layers) * ring_spacing_growth
        pos.update(circular_positions(extra_ring_nodes, radius, angle_offset=angle_offset))

    return pos


def layer_color_map(num_layers: int):
    """
    Inner layers are darker/more important.
    Outer layers become lighter and less opaque.
    """
    base_colors = [
        (0.15, 0.32, 0.70, 1.00),
        (0.30, 0.50, 0.82, 0.90),
        (0.52, 0.68, 0.90, 0.75),
        (0.72, 0.82, 0.94, 0.60),
        (0.86, 0.91, 0.97, 0.45),
    ]

    if num_layers <= len(base_colors):
        return base_colors[:num_layers]

    return base_colors + [base_colors[-1]] * (num_layers - len(base_colors))


def main():
    prefix = f"interactions_{NETWORK}"
    rank_csv = Path(f"{prefix}_degree_rankings.csv")
    tsv_path = Path(f"{prefix}.tsv")

    if not rank_csv.exists():
        raise SystemExit(f"Missing ranking CSV: {rank_csv.name} (run your degree ranking script first)")
    if not tsv_path.exists():
        raise SystemExit(f"Missing TSV: {tsv_path.name}")

    ranked_genes = load_ranked_genes(rank_csv)
    layers = split_into_layers(ranked_genes, LAYERS)

    # Genes included in layers
    layer_genes: Set[str] = set(g for layer in layers for g in layer)

    G = build_graph_from_tsv(tsv_path)

    # Expand to show connections around layer genes
    neighbors = get_one_step_neighbors(G, layer_genes)
    extra_nodes = sorted(list(neighbors - layer_genes))

    # Subnetwork
    keep_nodes = set(layer_genes) | neighbors
    sub = G.subgraph(keep_nodes).copy()

    # Output folder
    output_dir = Path(OUTPUT_DIR_NAME)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output file name includes layer definition
    layer_tag = "_".join("rest" if x >= 9999 else str(x) for x in LAYERS)
    output_svg = output_dir / f"new_{prefix}_degree_layers_{layer_tag}.svg"

    # Layout positions
    extra_ring_nodes = extra_nodes if NEIGHBOR_EXTRA_RING else []
    pos = layered_radial_layout(
        layers=[sorted(layer) for layer in layers],
        extra_ring_nodes=extra_ring_nodes,
        inner_radius=INNER_RADIUS,
        ring_spacing_base=RING_SPACING_BASE,
        ring_spacing_growth=RING_SPACING_GROWTH,
    )

    # If not using extra neighbor ring, place any missing nodes on an additional last ring
    if not NEIGHBOR_EXTRA_RING:
        missing = [n for n in sub.nodes() if n not in pos]
        if missing:
            r_last = INNER_RADIUS + (len(layers) - 1) * RING_SPACING_BASE
            pos.update(circular_positions(sorted(missing), r_last + RING_SPACING_BASE, angle_offset=0.0))

    # Colors: each layer gets progressively lighter
    layer_colors = layer_color_map(len(layers))

    gene_to_layer: Dict[str, int] = {}
    for li, layer in enumerate(layers):
        for g in layer:
            gene_to_layer[g] = li

    neighbor_color = "#E6EEF5"
    node_colors = []
    for n in sub.nodes():
        if n in gene_to_layer:
            node_colors.append(layer_colors[gene_to_layer[n]])
        else:
            node_colors.append(neighbor_color)

    # Node sizes: inner layers larger, outer layers smaller
    node_sizes: Dict[str, int] = {}
    for li, layer in enumerate(layers):
        for g in layer:
            node_sizes[g] = max(450, 1350 - li * 220)

    # Neighbor-only nodes small
    for n in sub.nodes():
        if n not in node_sizes:
            node_sizes[n] = 280

    sizes = [node_sizes[n] for n in sub.nodes()]

    # Draw SVG
    plt.figure(figsize=(15, 12))

    nx.draw_networkx_edges(sub, pos, alpha=0.35, width=1.1, edge_color="gray")

    nx.draw_networkx_nodes(
        sub,
        pos,
        node_size=sizes,
        node_color=node_colors,
        edgecolors="black",
        linewidths=0.8,
    )

    # Labels centered on nodes, but readable
    for node, (x, y) in pos.items():
        txt = plt.text(
            x,
            y,
            node,
            fontsize=10 if node in layer_genes else 9,
            color="black",
            ha="center",
            va="center",
            fontweight="normal",
            zorder=10,
        )

        txt.set_path_effects([
            pe.Stroke(linewidth=3.0, foreground="white"),
            pe.Normal(),
        ])

    layer_desc = ", ".join(str(x) for x in LAYERS)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_svg, format="svg")
    plt.close()

    print(f"Built subnetwork: nodes={sub.number_of_nodes()} edges={sub.number_of_edges()}")
    print(f"Layer genes total: {len(layer_genes)}")
    print(f"Neighbor-only nodes: {len(set(sub.nodes()) - layer_genes)}")
    print(f"Saved SVG: {output_svg}")


if __name__ == "__main__":
    main()