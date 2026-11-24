#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import json
import math
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

# ----------------------- DEFAULT GROUP FILES ---------------------------------
GROUP_FILES = {
    1: Path("group1_genes.txt"),
    2: Path("group2_genes.txt"),
    3: Path("group3_genes.txt"),
    4: Path("group4_genes.txt"),
}

# ----------------------- Paths / canvas --------------------------------------
INFO_DIR = Path("info_files")
SVG_DIR = Path("svg_files")
CACHE_FILE = INFO_DIR / "string_cache.json"

OUT_LEGEND = INFO_DIR / "rp_overlap_category_vertical_legend.txt"

CANVAS_W = 2200
CANVAS_H = 3200

NODE_R = 11
LABEL_FONT = 13

DEGREE_LAYER_RANKING_FILE = Path("interactions_g3_degree_rankings.csv")
DEGREE_LAYERS = [2, 3, 4, 5, 9999]

DEGREE_LAYER_COLORS = [
    "#1f4fa8",  # darkest, most central
    "#3f73c7",
    "#6f9fe0",
    "#a9c7ee",
    "#d8e7f7",  # lightest, least central
]
DEGREE_UNRANKED_G3_COLOR = "#eeeeee"

GROUP_COLORS = {
    1: "#e74c3c",
    2: "#3498db",
    3: "#27ae60",
    4: "#8e44ad",
}

GROUP_TITLES = {
    1: "Group 1 - RPE",
    2: "Group 2 - Outer Segment",
    3: "Group 3 - Connecting Cilium",
    4: "Group 4 - Nucleus",
}

VISUAL_TIER_INFO = {}

# ----------------------- STRING / bridge config ------------------------------
STRING_SPECIES = 9606
MIN_SCORE = 0.40
BRIDGE_UNASSIGNED_COLOR = "#b39ddb"
INTERMEDIATE_COLOR = "#00897b"
BRIDGE_EDGE_COLOR = "#555555"
BRIDGE_EDGE_OPACITY = 0.40
BRIDGE_EDGE_WIDTH = 1.2

# Highlight-mode colors
HIGHLIGHT_COLOR = "#ff0000"      # red
ONE_STEP_COLOR = "#4dd0e1"       # light cyan
TWO_STEP_COLOR = "#ffd54f"       # light yellow/gold
NON_HIGHLIGHT_COLOR = "#bbbbbb"  # gray

# ----------------------- Helpers ---------------------------------------------
def ensure_dirs():
    INFO_DIR.mkdir(parents=True, exist_ok=True)
    SVG_DIR.mkdir(parents=True, exist_ok=True)


def load_gene_file(path: Path):
    genes, seen = [], set()
    for raw in path.read_text().splitlines():
        g = raw.strip().upper()
        if g and g not in seen:
            seen.add(g)
            genes.append(g)
    return genes


def load_gene_list(path: Path):
    genes, seen = [], set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        g = raw.strip().replace(",", "").replace("'", "").replace('"', "")
        up = g.upper()
        if up and up not in seen:
            seen.add(up)
            genes.append(up)
    return genes


def load_edges(path: Path, allowed):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(header)}

    edges = []
    seen = set()

    for row in lines[1:]:
        parts = row.split("\t")
        a = parts[idx["#node1"]].upper()
        b = parts[idx["node2"]].upper()

        if a in allowed and b in allowed and a != b:
            key = tuple(sorted((a, b)))
            if key not in seen:
                seen.add(key)
                score = float(parts[idx["combined_score"]])
                edges.append((a, b, score))
    return edges


def compute_circle_radius(n):
    if n <= 1:
        return 60
    est = max(95, int((n * 34) / (2 * math.pi)))
    return min(340, est)


def circular_positions(genes, cx, cy, r):
    pos = {}
    genes = sorted(genes)
    n = len(genes)

    if n == 1:
        pos[genes[0]] = (cx, cy)
        return pos

    for i, g in enumerate(genes):
        theta = 2 * math.pi * i / n - math.pi / 2
        pos[g] = (cx + r * math.cos(theta), cy + r * math.sin(theta))
    return pos


def get_highlight_center():
    return CANVAS_W * 0.50, CANVAS_H * 0.52


def apply_highlight_visual_layout(pos, groups, center_pos, radii, highlight_gene, highlight_gid,
                                  highlight_neighbors, highlight_two_step):
    global VISUAL_TIER_INFO
    VISUAL_TIER_INFO = {}

    if not highlight_gene or highlight_gid not in (1, 2, 3, 4):
        return pos

    def tier_radius(n):
        if n <= 1:
            return 0
        needed = int((n * (2 * NODE_R + LABEL_FONT + 14)) / (2 * math.pi))
        return max(55, min(220, needed))

    def place_tiers(gid, tier_order, top_limit=None, bottom_limit=None):
        cx, cy = center_pos[gid]
        base_r = radii[gid]

        active_tiers = [(name, genes) for name, genes in tier_order if genes]
        if not active_tiers:
            return []

        tier_rs = [tier_radius(len(genes)) for _, genes in active_tiers]

        label_clearance = NODE_R + LABEL_FONT + 30
        gaps = []
        for i in range(len(active_tiers) - 1):
            gaps.append(tier_rs[i] + tier_rs[i + 1] + label_clearance)

        total_span = sum(gaps)

        # In highlight/cross mode, side arms should extend sideways,
        # not downward. Top/bottom arms still stack vertically.
        use_horizontal_tiers = False
        if highlight_gid in (1, 2, 3, 4):
            hx, hy = get_highlight_center()
            dx_from_center = cx - hx
            dy_from_center = cy - hy
            use_horizontal_tiers = abs(dx_from_center) > abs(dy_from_center)

        if use_horizontal_tiers:
            hx, hy = get_highlight_center()
            direction = 1 if cx > hx else -1

            start_x = cx - direction * (total_span / 2)
            tier_xs = []
            current_x = start_x

            for i in range(len(active_tiers)):
                if i == 0:
                    tier_xs.append(current_x)
                else:
                    current_x += direction * gaps[i - 1]
                    tier_xs.append(current_x)

            tier_info = []

            for (tier_name, genes), tier_x, r_tier in zip(active_tiers, tier_xs, tier_rs):
                if len(genes) == 1:
                    pos[genes[0]] = (tier_x, cy)
                else:
                    pos.update(circular_positions(genes, tier_x, cy, r_tier))

                for g in genes:
                    VISUAL_TIER_INFO[g] = {
                        "gid": gid,
                        "tier": tier_name,
                        "cx": tier_x,
                        "cy": cy,
                        "r": r_tier,
                    }

                tier_info.append({
                    "name": tier_name,
                    "genes": genes,
                    "cx": tier_x,
                    "cy": cy,
                    "r": r_tier,
                    "top": cy - r_tier - LABEL_FONT - NODE_R,
                    "bottom": cy + r_tier + LABEL_FONT + NODE_R,
                })

            return tier_info

        # Original vertical tier stacking for top/down arms
        tier_ys = []
        start_y = cy - total_span / 2
        current_y = start_y

        for i in range(len(active_tiers)):
            if i == 0:
                tier_ys.append(current_y)
            else:
                current_y += gaps[i - 1]
                tier_ys.append(current_y)

        if top_limit is None:
            title_y = cy - base_r - 50
            top_limit = title_y + NODE_R + LABEL_FONT + 65

        if bottom_limit is None:
            bottom_limit = CANVAS_H - 80

        top_edge = tier_ys[0] - tier_rs[0] - LABEL_FONT - NODE_R
        if top_edge < top_limit:
            shift_down = top_limit - top_edge
            tier_ys = [y + shift_down for y in tier_ys]

        bottom_edge = tier_ys[-1] + tier_rs[-1] + LABEL_FONT + NODE_R
        if bottom_edge > bottom_limit:
            shift_up = bottom_edge - bottom_limit
            tier_ys = [y - shift_up for y in tier_ys]

        tier_info = []

        for (tier_name, genes), tier_y, r_tier in zip(active_tiers, tier_ys, tier_rs):
            if len(genes) == 1:
                pos[genes[0]] = (cx, tier_y)
            else:
                pos.update(circular_positions(genes, cx, tier_y, r_tier))

            for g in genes:
                VISUAL_TIER_INFO[g] = {
                    "gid": gid,
                    "tier": tier_name,
                    "cx": cx,
                    "cy": tier_y,
                    "r": r_tier,
                }

            tier_info.append({
                "name": tier_name,
                "genes": genes,
                "cx": cx,
                "cy": tier_y,
                "r": r_tier,
                "top": tier_y - r_tier - LABEL_FONT - NODE_R,
                "bottom": tier_y + r_tier + LABEL_FONT + NODE_R,
            })

        return tier_info

    cx, cy = center_pos[highlight_gid]
    r = radii[highlight_gid]
    highlight_center_x, highlight_center_y = get_highlight_center()

    top_group_genes = [g for g in groups[highlight_gid] if g != highlight_gene]

    top_other = sorted([
        g for g in top_group_genes
        if g not in highlight_neighbors and g not in highlight_two_step
    ])
    top_two_step = sorted([g for g in top_group_genes if g in highlight_two_step])
    top_one_step = sorted([g for g in top_group_genes if g in highlight_neighbors])

    top_tiers = [
        ("other", top_other),
        ("two", top_two_step),
        ("one", top_one_step),
    ]

    rough_highlight_y = highlight_center_y

    top_limit = max(40, cy - r - 10)
    bottom_limit = highlight_center_y - NODE_R - LABEL_FONT - 90

    top_tier_info = place_tiers(
        highlight_gid,
        top_tiers,
        top_limit=top_limit,
        bottom_limit=bottom_limit
    )

    orange_tier = None
    for info in top_tier_info:
        if info["name"] == "one":
            orange_tier = info
            break

    if orange_tier is not None:
        highlight_clearance = NODE_R + LABEL_FONT + 70
        highlight_y = orange_tier["bottom"] + highlight_clearance
    else:
        highlight_y = rough_highlight_y

    # Highlighted gene goes in the middle of the cross.
    pos[highlight_gene] = (highlight_center_x, highlight_center_y)

    for gid in (1, 2, 3, 4):
        if gid == highlight_gid:
            continue

        one_step = sorted([g for g in groups[gid] if g in highlight_neighbors])
        two_step = sorted([g for g in groups[gid] if g in highlight_two_step])
        other = sorted([
            g for g in groups[gid]
            if g not in highlight_neighbors and g not in highlight_two_step
        ])

        bottom_tiers = [
            ("one", one_step),
            ("two", two_step),
            ("other", other),
        ]

        place_tiers(gid, bottom_tiers)

    return pos


def get_layout_centers(highlight_gid=None):
    # Highlight mode: cross layout
    # highlighted gene sits in the middle, with four groups around it.
    if highlight_gid in (1, 2, 3, 4):
        hx, hy = get_highlight_center()

        arm_dist = CANVAS_H * 0.20

        centers = {
            highlight_gid: (hx, hy - arm_dist),
        }

        remaining = [gid for gid in (1, 2, 3, 4) if gid != highlight_gid]

        # Put Group 4 downward instead of Group 3
        if 4 in remaining:
            down_gid = 4
            remaining.remove(4)
        else:
            down_gid = remaining.pop(0)

        centers[down_gid] = (hx, hy + arm_dist * 0.60)

        if len(remaining) >= 1:
            centers[remaining[0]] = (hx - arm_dist, hy)

        if len(remaining) >= 2:
            centers[remaining[1]] = (hx + arm_dist, hy)

        return centers

    raise ValueError("Center gene must belong to one of the four groups")


# ----------------------- STRING helpers --------------------------------------
def load_cache():
    if CACHE_FILE.exists():
        try:
            return json.loads(CACHE_FILE.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def save_cache(cache):
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp = CACHE_FILE.with_suffix(".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, sort_keys=True)
    tmp.replace(CACHE_FILE)


def resolve_string_id(symbol: str):
    r = requests.get(
        "https://string-db.org/api/json/get_string_ids",
        params={"identifiers": symbol, "species": STRING_SPECIES, "limit": 1, "caller_identity": "rp_mapping"},
        timeout=30,
    )
    if r.status_code != 200:
        return None, None
    data = r.json()
    if not data:
        return None, None
    preferred = (data[0].get("preferredName") or symbol).upper()
    string_id = data[0].get("stringId")
    return preferred, string_id


def fetch_network(string_id: str, limit=500):
    r = requests.get(
        "https://string-db.org/api/json/network",
        params={"identifiers": string_id, "species": STRING_SPECIES, "limit": limit, "caller_identity": "rp_mapping"},
        timeout=60,
    )
    if r.status_code != 200:
        return []
    try:
        return r.json()
    except Exception:
        return []


def get_neighbors(symbol: str, cache: dict):
    key = f"NEI::{symbol.upper()}"
    if key in cache:
        return cache[key]
    pref, sid = resolve_string_id(symbol)
    if not sid:
        cache[key] = {"preferred": None, "edges": []}
        return cache[key]
    edges = fetch_network(sid, limit=500)
    cache[key] = {"preferred": pref, "edges": edges}
    return cache[key]


def fetch_many_neighbors(symbols, cache, max_workers=8):
    to_fetch = []
    for s in symbols:
        if not s:
            continue
        key = f"NEI::{s.upper()}"
        if key not in cache:
            to_fetch.append(s)

    if not to_fetch:
        return

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(get_neighbors, s, cache): s for s in to_fetch}
        for _ in as_completed(futures):
            pass


def find_intermediate_bridges(unassigned_genes, groups, cache, min_score):
    group_1_4 = set().union(groups[1], groups[2], groups[3], groups[4])

    fetch_many_neighbors(unassigned_genes, cache)

    candidate_m = set()
    for u in unassigned_genes:
        info_u = get_neighbors(u, cache)
        pref_u = info_u["preferred"]
        if not pref_u:
            continue

        for e in info_u["edges"]:
            a = (e.get("preferredName_A") or "").upper()
            b = (e.get("preferredName_B") or "").upper()
            score = float(e.get("score") or 0.0)
            if score < min_score:
                continue

            if a == pref_u and b:
                candidate_m.add(b)
            elif b == pref_u and a:
                candidate_m.add(a)

    fetch_many_neighbors(candidate_m, cache)

    best_paths = []

    for u in sorted(unassigned_genes):
        info_u = get_neighbors(u, cache)
        pref_u = info_u["preferred"]
        if not pref_u:
            continue

        best = None
        for e in info_u["edges"]:
            a = (e.get("preferredName_A") or "").upper()
            b = (e.get("preferredName_B") or "").upper()
            score_um = float(e.get("score") or 0.0)
            if score_um < min_score:
                continue

            if a == pref_u and b:
                m_name = b
            elif b == pref_u and a:
                m_name = a
            else:
                continue

            info_m = get_neighbors(m_name, cache)
            pref_m = info_m["preferred"]
            if not pref_m:
                continue
            m = pref_m

            if m in group_1_4:
                continue

            for e2 in info_m["edges"]:
                a2 = (e2.get("preferredName_A") or "").upper()
                b2 = (e2.get("preferredName_B") or "").upper()
                score_mg = float(e2.get("score") or 0.0)
                if score_mg < min_score:
                    continue

                if a2 == m and b2 in group_1_4:
                    ggene = b2
                elif b2 == m and a2 in group_1_4:
                    ggene = a2
                else:
                    continue

                candidate = (
                    score_um + score_mg,
                    min(score_um, score_mg),
                    max(score_um, score_mg),
                    u, m, ggene, score_um, score_mg
                )
                if best is None or candidate > best:
                    best = candidate

        if best is not None:
            _, _, _, u, m, ggene, score_um, score_mg = best
            best_paths.append((u, m, ggene, score_um, score_mg))

    return best_paths


def move_bridge_anchor_genes_to_sides(bridging_edges, positions, groups, centers, radii, gene_to_group):
    """
    Visual-only adjustment:
    genes that have intermediate bridges are swapped toward outward positions
    within their same displayed circle/tier.

    For normal/top/bottom circles, anchors are moved to left/right side positions.
    For left/right arms in highlight mode, anchors are moved to top/bottom positions
    so their intermediates can extend vertically instead of sideways.

    The displaced gene fills the bridge-anchor gene's old slot, so no two group
    genes stack.
    """
    bridge_anchors = []
    seen = set()
    for _, _, ggene, _, _ in bridging_edges:
        if ggene not in seen and ggene in positions:
            seen.add(ggene)
            bridge_anchors.append(ggene)

    bridge_anchor_set = set(bridge_anchors)
    used_target_genes = set()
    side_arm_anchor_counts = defaultdict(int)

    def display_circle_for_gene(gene, gid):
        if gene in VISUAL_TIER_INFO:
            t = VISUAL_TIER_INFO[gene]
            return t["cx"], t["cy"], max(t["r"], 1), t.get("tier")
        cx, cy = centers[gid]
        return cx, cy, radii[gid], None

    def is_left_or_right_arm(gene):
        if gene not in VISUAL_TIER_INFO:
            return False

        hx, hy = get_highlight_center()
        t = VISUAL_TIER_INFO[gene]
        tier_cx = t["cx"]
        tier_cy = t["cy"]

        return abs(tier_cx - hx) > abs(tier_cy - hy)

    for anchor in bridge_anchors:
        gid = gene_to_group.get(anchor)
        if gid not in groups or anchor not in positions:
            continue

        cx, cy, r, tier_name = display_circle_for_gene(anchor, gid)
        side_arm = is_left_or_right_arm(anchor)

        if tier_name is not None:
            candidates = [
                g for g in groups[gid]
                if g in positions
                and g != anchor
                and g not in bridge_anchor_set
                and g not in used_target_genes
                and g in VISUAL_TIER_INFO
                and VISUAL_TIER_INFO[g].get("tier") == tier_name
            ]
        else:
            candidates = [
                g for g in groups[gid]
                if g in positions
                and g != anchor
                and g not in bridge_anchor_set
                and g not in used_target_genes
            ]

        if not candidates:
            continue

        if side_arm:
            # For left/right arms, move anchor genes to top/bottom positions.
            idx = side_arm_anchor_counts[(gid, tier_name)]
            side_arm_anchor_counts[(gid, tier_name)] += 1
            want_bottom = (idx % 2 == 0)

            def target_score(g):
                x, y = positions[g]
                dx = x - cx
                dy = y - cy

                if want_bottom:
                    vertical_pref = -dy
                else:
                    vertical_pref = dy

                horizontal_penalty = abs(dx) * 0.35
                edge_penalty = abs(math.hypot(dx, dy) - r) * 0.25

                return vertical_pref + horizontal_penalty + edge_penalty

        else:
            # For top/bottom arms, keep the side-position behavior.
            idx = side_arm_anchor_counts[(gid, tier_name)]
            side_arm_anchor_counts[(gid, tier_name)] += 1
            want_right = (idx % 2 == 0)

            def target_score(g):
                x, y = positions[g]
                dx = x - cx
                dy = y - cy

                if want_right:
                    side_pref = -dx
                else:
                    side_pref = dx

                vertical_penalty = abs(dy) * 0.35
                edge_penalty = abs(math.hypot(dx, dy) - r) * 0.25

                return side_pref + vertical_penalty + edge_penalty

        target = min(candidates, key=target_score)

        old_anchor_pos = positions[anchor]
        old_target_pos = positions[target]

        positions[anchor] = old_target_pos
        positions[target] = old_anchor_pos

        used_target_genes.add(target)

    return positions

def place_bridged_nodes(bridging_edges, positions, centers, radii, gene_to_group):
    extra_positions = {}
    anchor_counts = defaultdict(int)
    intermediate_counts = defaultdict(int)
    placed_bridge_segments = []

    anchor_to_bridges = defaultdict(list)
    for bridge in bridging_edges:
        u, m, ggene, score_um, score_mg = bridge
        anchor_to_bridges[ggene].append(bridge)
        intermediate_counts[m] += 1

    outside_gap = NODE_R + LABEL_FONT + 26
    unassigned_gap = NODE_R + LABEL_FONT + 34
    step_perp = NODE_R + LABEL_FONT + 18
    step_out = NODE_R + LABEL_FONT + 16

    occupied = set()

    def reserve_existing_positions():
        for x, y in positions.values():
            occupied.add((int(round(x)), int(round(y))))

    def too_close(x, y, min_dist=None):
        if min_dist is None:
            min_dist = (2 * NODE_R) + LABEL_FONT + 10

        for ox, oy in occupied:
            if (x - ox) ** 2 + (y - oy) ** 2 < min_dist ** 2:
                return True
        return False

    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

    def segments_intersect(a, b, c, d):
        shared_dist = 20
        for p in (a, b):
            for q in (c, d):
                if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 < shared_dist ** 2:
                    return False

        return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)

    def crosses_existing_segment(x1, y1, x2, y2):
        new_seg = ((x1, y1), (x2, y2))
        for old_a, old_b in placed_bridge_segments:
            if segments_intersect(new_seg[0], new_seg[1], old_a, old_b):
                return True
        return False

    def find_free_spot(x, y):
        x = int(round(x))
        y = int(round(y))

        if not too_close(x, y):
            occupied.add((x, y))
            return x, y

        search_step = NODE_R + LABEL_FONT + 6
        for radius in [search_step * i for i in range(1, 8)]:
            for dx, dy in [
                ( radius, 0), (-radius, 0), (0,  radius), (0, -radius),
                ( radius,  radius), ( radius, -radius),
                (-radius,  radius), (-radius, -radius),
                ( radius,  radius // 2), ( radius, -radius // 2),
                (-radius,  radius // 2), (-radius, -radius // 2),
                ( radius // 2,  radius), (-radius // 2,  radius),
                ( radius // 2, -radius), (-radius // 2, -radius),
            ]:
                nx, ny = x + dx, y + dy
                if not too_close(nx, ny):
                    occupied.add((nx, ny))
                    return nx, ny

        occupied.add((x, y))
        return x, y

    reserve_existing_positions()

    bridges_sorted = sorted(
        bridging_edges,
        key=lambda t: (gene_to_group[t[2]], t[2], t[0], t[1])
    )

    for u, m, ggene, score_um, score_mg in bridges_sorted:
        gid = gene_to_group[ggene]
        gx, gy = positions[ggene]

        if ggene in VISUAL_TIER_INFO:
            tier = VISUAL_TIER_INFO[ggene]
            cx = tier["cx"]
            cy = tier["cy"]
        else:
            cx, cy = centers[gid]

        # For left/right arms, bridges extend vertically from the anchor.
        # For top/bottom arms, bridges extend sideways.
        if ggene in VISUAL_TIER_INFO:
            hx, hy = get_highlight_center()
            tier = VISUAL_TIER_INFO[ggene]
            tier_cx = tier["cx"]
            tier_cy = tier["cy"]

            side_arm = abs(tier_cx - hx) > abs(tier_cy - hy)

            if side_arm:
                if gy >= tier_cy:
                    out_x, out_y = 0.0, 1.0
                else:
                    out_x, out_y = 0.0, -1.0

                perp_x, perp_y = 1.0, 0.0
            else:
                if gx >= cx:
                    out_x, out_y = 1.0, 0.0
                else:
                    out_x, out_y = -1.0, 0.0

                perp_x, perp_y = 0.0, 1.0
        else:
            if gx >= cx:
                out_x, out_y = 1.0, 0.0
            else:
                out_x, out_y = -1.0, 0.0

            perp_x, perp_y = 0.0, 1.0

        idx = anchor_counts[ggene]
        anchor_counts[ggene] += 1

        n_for_anchor = len(anchor_to_bridges[ggene])
        centered_idx = idx - (n_for_anchor - 1) / 2.0

        if n_for_anchor == 1:
            perp_offset = 0
            outward_extra = 0
        else:
            perp_offset = centered_idx * step_perp
            outward_extra = abs(centered_idx) * step_out

        mx = gx + out_x * (outside_gap + outward_extra) + perp_x * perp_offset
        my = gy + out_y * (outside_gap + outward_extra) + perp_y * perp_offset

        ux_pos = gx + out_x * (outside_gap + unassigned_gap + outward_extra) + perp_x * perp_offset
        uy_pos = gy + out_y * (outside_gap + unassigned_gap + outward_extra) + perp_y * perp_offset

        if m not in extra_positions:
            mx, my = find_free_spot(mx, my)
            extra_positions[m] = (mx, my)
        else:
            mx, my = extra_positions[m]

        if u not in extra_positions:
            if intermediate_counts[m] > 1:
                used_for_m = sum(
                    1 for existing_u, existing_m, _, _, _ in bridges_sorted
                    if existing_m == m and existing_u < u
                )
                fan_idx = used_for_m - (intermediate_counts[m] - 1) / 2.0
                ux_pos = mx + out_x * unassigned_gap + perp_x * fan_idx * step_perp
                uy_pos = my + out_y * unassigned_gap + perp_y * fan_idx * step_perp
            else:
                ux_pos = mx + out_x * unassigned_gap
                uy_pos = my + out_y * unassigned_gap

            if crosses_existing_segment(mx, my, ux_pos, uy_pos):
                for mult in (1, -1, 2, -2, 3, -3):
                    test_x = ux_pos + perp_x * mult * step_perp
                    test_y = uy_pos + perp_y * mult * step_perp

                    if not crosses_existing_segment(mx, my, test_x, test_y):
                        ux_pos, uy_pos = test_x, test_y
                        break

            ux_pos, uy_pos = find_free_spot(ux_pos, uy_pos)
            extra_positions[u] = (ux_pos, uy_pos)

        if u in extra_positions and m in extra_positions:
            placed_bridge_segments.append(
                (extra_positions[m], extra_positions[u])
            )

    return extra_positions

def compute_label_offsets(positions, label_genes):
    """Returns {gene: (dx, dy)}. Labels may cross lines, never circles or other labels."""
    node_boxes = {g: (x - NODE_R, y - NODE_R, x + NODE_R, y + NODE_R)
                  for g, (x, y) in positions.items()}
    placed_label_boxes = []
    label_offsets = {}
    unresolved = 0

    def text_box(x, y_baseline, gene, pad=2):
        half_w = max(18, len(gene) * LABEL_FONT * 0.34) + pad
        return (x - half_w, y_baseline - LABEL_FONT - pad,
                x + half_w, y_baseline + 3 + pad)

    def overlaps(a, b):
        return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])

    def candidates(gene):
        half_w = max(18, len(gene) * LABEL_FONT * 0.34)
        above = -(NODE_R + 5)
        below = NODE_R + LABEL_FONT + 2
        side = NODE_R + 4 + half_w
        c = []
        for dy in (above, below):                       # over / under the circle,
            for dx in (0, 10, -10, 20, -20, 32, -32, 46, -46):   # nudged sideways
                c.append((dx, dy))
        c += [(side, 4), (-side, 4)]                    # beside the circle
        for dy in (above - 16, below + 16):             # second row
            for dx in (0, 12, -12, 28, -28):
                c.append((dx, dy))
        c += [(side, -10), (-side, -10), (side, 18), (-side, 18)]
        return c

    def crowding(gene):
        x, y = positions[gene]
        return sum(1 for k, (a, b) in positions.items()
                   if k != gene and (a - x) ** 2 + (b - y) ** 2 < 70 ** 2)

    # Most crowded genes choose first, since they have the fewest options.
    ordered = sorted(label_genes,
                     key=lambda g: (-crowding(g), positions[g][1], positions[g][0], g))

    for gene in ordered:
        x, y = positions[gene]
        best = None
        for dx, dy in candidates(gene):
            box = text_box(x + dx, y + dy, gene)
            n = (sum(overlaps(box, b) for b in placed_label_boxes)
                 + sum(overlaps(box, b) for g, b in node_boxes.items()
                       if g != gene))   # own circle excluded: every slot already clears it
            if best is None or n < best[0]:
                best = (n, dx, dy, box)
            if n == 0:          # first clean slot wins (closest to the node)
                break
        if best[0] > 0:
            unresolved += 1
        placed_label_boxes.append(best[3])
        label_offsets[gene] = (best[1], best[2])

    if unresolved:
        print(f"[WARN] {unresolved} labels could not avoid overlapping a circle or label.")
    return label_offsets


# ----------------------- SVG --------------------------------------------------
def svg_edge(x1, y1, x2, y2, color, opacity, width):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-opacity="{opacity}" stroke-width="{width}"/>\n'


def svg_circle(x, y, fill):
    return f'<circle cx="{x}" cy="{y}" r="{NODE_R}" fill="{fill}" stroke="#222" stroke-width="0.8"/>\n'


def _label_xy(x, y, offset):
    if offset is None:
        return x, y - (NODE_R + 5)
    if isinstance(offset, (tuple, list)):
        return x + offset[0], y + offset[1]
    return x, y + offset          # old-style single dy


def svg_text(x, y, gene, offset=None):
    lx, ly = _label_xy(x, y, offset)
    return f'<text x="{int(lx)}" y="{int(ly)}" font-size="{LABEL_FONT}" text-anchor="middle">{gene}</text>\n'

def load_group3_degree_layer_colors():
    """
    Loads interactions_g3_degree_rankings.csv and assigns colors by degree layer.
    Layer 1 gets darkest color, outer/rest layers get lighter colors.
    """
    if not DEGREE_LAYER_RANKING_FILE.exists():
        print(f"[WARN] Missing {DEGREE_LAYER_RANKING_FILE}; using normal Group 3 color.")
        return {}

    lines = DEGREE_LAYER_RANKING_FILE.read_text(encoding="utf-8").splitlines()
    if not lines:
        return {}

    header = lines[0].split(",")
    if "Gene" not in header:
        print(f"[WARN] {DEGREE_LAYER_RANKING_FILE} missing Gene column.")
        return {}

    gene_idx = header.index("Gene")

    ranked_genes = []
    for row in lines[1:]:
        parts = row.split(",")
        if len(parts) > gene_idx:
            ranked_genes.append(parts[gene_idx].strip().upper())

    gene_to_color = {}
    idx = 0

    for layer_i, size in enumerate(DEGREE_LAYERS):
        color = DEGREE_LAYER_COLORS[min(layer_i, len(DEGREE_LAYER_COLORS) - 1)]

        if size >= 9999:
            layer_genes = ranked_genes[idx:]
        else:
            layer_genes = ranked_genes[idx: idx + size]

        for gene in layer_genes:
            gene_to_color[gene] = color

        idx += size
        if idx >= len(ranked_genes):
            break

    return gene_to_color


def svg_text_halo(x, y, gene, offset=None):
    lx, ly = _label_xy(x, y, offset)
    return (
        f'<text x="{int(lx)}" y="{int(ly)}" font-size="{LABEL_FONT}" text-anchor="middle" '
        f'font-family="Arial" stroke="white" stroke-width="3" '
        f'stroke-linejoin="round">{gene}</text>\n'
        f'<text x="{int(lx)}" y="{int(ly)}" font-size="{LABEL_FONT}" text-anchor="middle" '
        f'font-family="Arial" fill="black">{gene}</text>\n'
    )


def build_group3_isolated_svg(pos, groups, edges, bridging_edges, intermediates, bridged_unassigned, label_offsets):
    """
    Builds a separate SVG containing only Group 3 genes and the intermediate
    bridge structures attached to Group 3 genes.

    Positions/connections are copied from the full intermediate diagram.
    Group 3 gene colors come from the degree-layer ranking diagram.
    """
    group3_genes = set(groups[3])

    group3_bridges = [
        bridge for bridge in bridging_edges
        if bridge[2] in group3_genes
    ]

    group3_intermediates = {
        m for _, m, _, _, _ in group3_bridges
        if m in pos
    }

    group3_unassigned = {
        u for u, _, _, _, _ in group3_bridges
        if u in pos
    }

    keep_nodes = set(group3_genes) | group3_intermediates | group3_unassigned
    keep_nodes = {g for g in keep_nodes if g in pos}

    if not keep_nodes:
        print("[WARN] No Group 3 nodes found for isolated SVG.")
        return

    xs = [pos[g][0] for g in keep_nodes]
    ys = [pos[g][1] for g in keep_nodes]

    padding = 180
    min_x = min(xs) - padding
    max_x = max(xs) + padding
    min_y = min(ys) - padding
    max_y = max(ys) + padding

    out_w = int(max_x - min_x)
    out_h = int(max_y - min_y)

    def shift_xy(gene):
        x, y = pos[gene]
        return int(x - min_x), int(y - min_y)

    degree_colors = load_group3_degree_layer_colors()

    out_svg = SVG_DIR / "group3_isolated_intermediates_degree_colors.svg"

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{out_w}" height="{out_h}">\n')

    # Title
    lines.append(
        f'<text x="{out_w // 2}" y="45" text-anchor="middle" '
        f'font-family="Arial" font-size="26" fill="#111">'
        f'Group 3 - Connecting Cilium</text>\n'
    )

    # Original Group 3 internal edges
    for a, b, score in edges:
        if a in group3_genes and b in group3_genes and a in pos and b in pos:
            x1, y1 = shift_xy(a)
            x2, y2 = shift_xy(b)
            lines.append(svg_edge(x1, y1, x2, y2, "#777777", 0.22, 1.0))

    # Group 3 bridge edges
    for u, m, ggene, score_um, score_mg in group3_bridges:
        if ggene in group3_genes and m in pos and ggene in pos:
            x1, y1 = shift_xy(m)
            x2, y2 = shift_xy(ggene)
            lines.append(svg_edge(x1, y1, x2, y2, BRIDGE_EDGE_COLOR, BRIDGE_EDGE_OPACITY, BRIDGE_EDGE_WIDTH))

        if u in pos and m in pos:
            x1, y1 = shift_xy(u)
            x2, y2 = shift_xy(m)
            lines.append(svg_edge(x1, y1, x2, y2, BRIDGE_EDGE_COLOR, BRIDGE_EDGE_OPACITY, BRIDGE_EDGE_WIDTH))

    # Group 3 nodes, colored by degree-layer diagram
    for gene in sorted(group3_genes):
        if gene in pos:
            x, y = shift_xy(gene)
            fill = degree_colors.get(gene, DEGREE_UNRANKED_G3_COLOR)
            lines.append(svg_circle(x, y, fill))

    # Intermediate nodes
    for gene in sorted(group3_intermediates):
        x, y = shift_xy(gene)
        lines.append(svg_circle(x, y, INTERMEDIATE_COLOR))

    # Bridged unassigned nodes
    for gene in sorted(group3_unassigned):
        x, y = shift_xy(gene)
        lines.append(svg_circle(x, y, BRIDGE_UNASSIGNED_COLOR))

    # Labels last
    for gene in sorted(keep_nodes):
        x, y = shift_xy(gene)
        lines.append(svg_text_halo(x, y, gene, label_offsets.get(gene)))

    lines.append("</svg>")
    out_svg.write_text("".join(lines), encoding="utf-8")

    print(f"[OK] Wrote isolated Group 3 SVG: {out_svg}")


def build_group3_isolated_rho_colors_svg(
    pos, groups, edges, bridging_edges,
    intermediates, bridged_unassigned,
    label_offsets, highlight_gene="RHO"
):
    """
    Generate an isolated Group 3 SVG with the same layout,
    connections, and RHO-based colors as the full diagram.
    """
    group3_genes = set(groups[3])

    group3_bridges = [
        bridge for bridge in bridging_edges
        if bridge[2] in group3_genes
    ]

    group3_intermediates = {
        m for _, m, _, _, _ in group3_bridges if m in pos
    }
    group3_unassigned = {
        u for u, _, _, _, _ in group3_bridges if u in pos
    }

    keep_nodes = (
        group3_genes | group3_intermediates | group3_unassigned
    )
    keep_nodes = {g for g in keep_nodes if g in pos}

    if not keep_nodes:
        return

    xs = [pos[g][0] for g in keep_nodes]
    ys = [pos[g][1] for g in keep_nodes]

    padding = 180
    min_x = min(xs) - padding
    max_x = max(xs) + padding
    min_y = min(ys) - padding
    max_y = max(ys) + padding

    out_w = int(max_x - min_x)
    out_h = int(max_y - min_y)

    def shift_xy(gene):
        x, y = pos[gene]
        return int(x - min_x), int(y - min_y)

    # Reconstruct the full diagram's RHO connectivity colors.
    adjacency = defaultdict(set)
    for a, b, score in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    highlight_gene = highlight_gene.upper()
    one_step = set(adjacency.get(highlight_gene, set()))

    two_step = set()
    for neighbor in one_step:
        two_step.update(adjacency.get(neighbor, set()))

    two_step.discard(highlight_gene)
    two_step -= one_step

    out_svg = (
        SVG_DIR /
        "group3_isolated_intermediates_RHO_colors.svg"
    )

    lines = []
    lines.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{out_w}" height="{out_h}">\n'
    )

    lines.append(
        f'<text x="{out_w // 2}" y="45" '
        f'text-anchor="middle" font-family="Arial" '
        f'font-size="26" fill="#111">'
        f'Group 3 - Connecting Cilium</text>\n'
    )

    # Same internal Group 3 edges as the existing close-up.
    for a, b, score in edges:
        if a in group3_genes and b in group3_genes:
            if a in pos and b in pos:
                x1, y1 = shift_xy(a)
                x2, y2 = shift_xy(b)
                lines.append(
                    svg_edge(
                        x1, y1, x2, y2,
                        "#777777", 0.22, 1.0
                    )
                )

    # Same intermediate bridging edges.
    for u, m, ggene, score_um, score_mg in group3_bridges:
        if m in pos and ggene in pos:
            x1, y1 = shift_xy(m)
            x2, y2 = shift_xy(ggene)
            lines.append(
                svg_edge(
                    x1, y1, x2, y2,
                    BRIDGE_EDGE_COLOR,
                    BRIDGE_EDGE_OPACITY,
                    BRIDGE_EDGE_WIDTH
                )
            )

        if u in pos and m in pos:
            x1, y1 = shift_xy(u)
            x2, y2 = shift_xy(m)
            lines.append(
                svg_edge(
                    x1, y1, x2, y2,
                    BRIDGE_EDGE_COLOR,
                    BRIDGE_EDGE_OPACITY,
                    BRIDGE_EDGE_WIDTH
                )
            )

    # Group 3 colors identical to the full diagram.
    for gene in sorted(group3_genes):
        if gene not in pos:
            continue

        x, y = shift_xy(gene)

        if gene == highlight_gene:
            fill = HIGHLIGHT_COLOR
        elif gene in one_step:
            fill = ONE_STEP_COLOR
        elif gene in two_step:
            fill = TWO_STEP_COLOR
        else:
            fill = NON_HIGHLIGHT_COLOR

        lines.append(svg_circle(x, y, fill))

    # Intermediate and unassigned colors remain unchanged.
    for gene in sorted(group3_intermediates):
        x, y = shift_xy(gene)
        lines.append(svg_circle(x, y, INTERMEDIATE_COLOR))

    for gene in sorted(group3_unassigned):
        x, y = shift_xy(gene)
        lines.append(svg_circle(x, y, BRIDGE_UNASSIGNED_COLOR))

    # Same label placement.
    for gene in sorted(keep_nodes):
        x, y = shift_xy(gene)
        lines.append(
            svg_text_halo(
                x, y, gene, label_offsets.get(gene)
            )
        )

    lines.append("</svg>")
    out_svg.write_text("".join(lines), encoding="utf-8")

    print(f"[OK] Wrote isolated Group 3 RHO-colored SVG: {out_svg}")


# ----------------------- Main -------------------------------------------------
def build_diagram(edge_file, genes_file, highlight_gene="RHO", min_score=MIN_SCORE):
    groups = {}
    gene_to_group = {}
    all_genes = set()

    for gid, path in GROUP_FILES.items():
        genes = load_gene_file(path)
        groups[gid] = genes
        for g in genes:
            gene_to_group[g] = gid
            all_genes.add(g)

    edges = load_edges(edge_file, all_genes)

    highlight_gene = highlight_gene.upper()
    if highlight_gene not in gene_to_group:
        raise ValueError(f"Center gene {highlight_gene} is not in group1_genes.txt through group4_genes.txt")

    adjacency = defaultdict(set)
    for a, b, _ in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    highlight_neighbors = set(adjacency.get(highlight_gene, set()))
    highlight_two_step = set()
    for n1 in highlight_neighbors:
        highlight_two_step.update(adjacency.get(n1, set()))
    highlight_two_step.discard(highlight_gene)
    highlight_two_step -= highlight_neighbors
    highlight_gid = gene_to_group[highlight_gene]

    center_pos = get_layout_centers(highlight_gid)
    radii = {gid: compute_circle_radius(len(groups[gid])) for gid in (1, 2, 3, 4)}

    pos = {}
    for gid in (1, 2, 3, 4):
        cx, cy = center_pos[gid]
        pos.update(circular_positions(groups[gid], cx, cy, radii[gid]))

    pos = apply_highlight_visual_layout(
        pos,
        groups,
        center_pos,
        radii,
        highlight_gene,
        highlight_gid,
        highlight_neighbors,
        highlight_two_step
    )

    # ---------- bridge logic added on top of old-paper base ----------
    user_genes = set(load_gene_list(Path(genes_file)))
    group_genes = set(all_genes)
    unassigned_genes = sorted(user_genes - group_genes)

    cache = load_cache()
    n_cached_before = len(cache)
    bridging_edges = find_intermediate_bridges(
        unassigned_genes,
        {gid: set(groups[gid]) for gid in groups},
        cache,
        min_score
    )
    n_new = len(cache) - n_cached_before
    if n_new > 0:
        print(f"[INFO] Fetched {n_new} new genes from STRING; updating cache.")
        save_cache(cache)
    else:
        print("[INFO] All genes found in the existing cache; cache file left untouched.")

    bridged_unassigned = {u for u, _, _, _, _ in bridging_edges}
    intermediates = {m for _, m, _, _, _ in bridging_edges if m not in group_genes and m not in unassigned_genes}

    pos = move_bridge_anchor_genes_to_sides(
        bridging_edges,
        pos,
        groups,
        center_pos,
        radii,
        gene_to_group
    )

    extra_positions = place_bridged_nodes(bridging_edges, pos, center_pos, radii, gene_to_group)
    pos.update(extra_positions)

    out_txt = INFO_DIR / "unassigned_with_intermediates.txt"
    if bridging_edges:
        with out_txt.open("w", encoding="utf-8") as f:
            header = (
                f"{'U':<17}"
                f"{'M':<22}"
                f"{'G':<22}"
                f"{'Score_unassigned→intermediate':<28}"
                f"{'Score_intermediate→group':<28}\n"
            )
            f.write(header)
            f.write("-" * len(header) + "\n")
            for u, m, ggene, score_um, score_mg in sorted(bridging_edges):
                f.write(
                    f"{u:<17}"
                    f"{m:<22}"
                    f"{ggene:<22}"
                    f"{score_um:<28.3f}"
                    f"{score_mg:<28.3f}\n"
                )
    else:
        out_txt.write_text("No valid intermediate bridges found.\n", encoding="utf-8")

    out_svg = SVG_DIR / f"rp_with_intermediates_new_highlight_{highlight_gene}.svg"

    label_genes = set()
    for gid in (1, 2, 3, 4):
        label_genes.update(groups[gid])
    label_genes.update(intermediates)
    label_genes.update(bridged_unassigned)
    label_genes = {g for g in label_genes if g in pos}
    label_offsets = compute_label_offsets(pos, label_genes)

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{CANVAS_W}" height="{CANVAS_H}">\n')

    # group labels
    for gid in (1, 2, 3, 4):
        cx, cy = center_pos[gid]
        r = radii[gid]

        hx, hy = get_highlight_center()
        dx = cx - hx
        dy = cy - hy

        # Find visible extent of all tier circles belonging to this group.
        group_tier_infos = [
            info for info in VISUAL_TIER_INFO.values()
            if info["gid"] == gid
        ]

        if group_tier_infos:
            top_extent = min(info["cy"] - info["r"] for info in group_tier_infos)
            bottom_extent = max(info["cy"] + info["r"] for info in group_tier_infos)
            left_extent = min(info["cx"] - info["r"] for info in group_tier_infos)
            right_extent = max(info["cx"] + info["r"] for info in group_tier_infos)
        else:
            top_extent = cy - r
            bottom_extent = cy + r
            left_extent = cx - r
            right_extent = cx + r

        label_padding = 70

        if gid in (1, 3) and abs(dx) > abs(dy):
            arm_genes = set(groups[gid])
            for u, m, ggene, _, _ in bridging_edges:
                if gene_to_group.get(ggene) == gid:
                    arm_genes.update((u, m))
            arm_genes = [g for g in arm_genes if g in pos]

            arm_left = min(pos[g][0] - NODE_R for g in arm_genes)
            arm_right = max(pos[g][0] + NODE_R for g in arm_genes)
            arm_top = min(
                min(pos[g][1] - NODE_R,
                    pos[g][1] + label_offsets.get(g, (0, -(NODE_R + 5)))[1] - LABEL_FONT)
                for g in arm_genes
            )
            label_x = (arm_left + arm_right) / 2
            label_y = arm_top - 45  
            # Keep even the longer Connecting Cilium title inside the SVG.
            estimated_half_width = len(GROUP_TITLES[gid]) * 24 * 0.31
            label_x = max(estimated_half_width + 20,
                          min(CANVAS_W - estimated_half_width - 20, label_x))
            anchor = "middle"

        elif abs(dy) >= abs(dx) and dy < 0:
            # top arm
            label_x = cx
            label_y = top_extent - label_padding
            anchor = "middle"

        elif abs(dy) >= abs(dx) and dy > 0:
            # bottom arm
            label_x = cx
            label_y = bottom_extent + label_padding + 20
            anchor = "middle"

        elif dx < 0:
            # left arm
            label_x = left_extent - label_padding
            label_y = cy
            anchor = "end"

        else:
            # right arm
            label_x = right_extent + label_padding
            label_y = cy
            anchor = "start"

        lines.append(
            f'<text x="{int(label_x)}" y="{int(label_y)}" text-anchor="{anchor}" '
            f'font-family="Arial" font-size="24" fill="#111">{GROUP_TITLES[gid]}</text>\n'
        )

    # original edges
    for a, b, score in edges:
        x1, y1 = pos[a]
        x2, y2 = pos[b]

        if a == highlight_gene or b == highlight_gene:
            color = HIGHLIGHT_COLOR
            opacity = 0.62
            width = 2.0
        elif (a in highlight_neighbors and b in highlight_two_step) or \
             (b in highlight_neighbors and a in highlight_two_step):
            color = ONE_STEP_COLOR
            opacity = 0.34
            width = 1.4
        else:
            color = "#777777"
            opacity = 0.24
            width = 0.8

        lines.append(svg_edge(int(x1), int(y1), int(x2), int(y2), color, opacity, width))

    # bridge edges
    for u, m, ggene, score_um, score_mg in bridging_edges:
        if u in pos and m in pos:
            x1, y1 = pos[u]
            x2, y2 = pos[m]
            lines.append(svg_edge(int(x1), int(y1), int(x2), int(y2), BRIDGE_EDGE_COLOR, BRIDGE_EDGE_OPACITY, BRIDGE_EDGE_WIDTH))
        if m in pos and ggene in pos:
            x1, y1 = pos[m]
            x2, y2 = pos[ggene]
            lines.append(svg_edge(int(x1), int(y1), int(x2), int(y2), BRIDGE_EDGE_COLOR, BRIDGE_EDGE_OPACITY, BRIDGE_EDGE_WIDTH))

    # ---------- draw ALL circles first ----------
    for gid in (1, 2, 3, 4):
        for gene in groups[gid]:
            x, y = pos[gene]

            if gene == highlight_gene:
                fill = HIGHLIGHT_COLOR
            elif gene in highlight_neighbors:
                fill = ONE_STEP_COLOR
            elif gene in highlight_two_step:
                fill = TWO_STEP_COLOR
            else:
                fill = NON_HIGHLIGHT_COLOR

            lines.append(svg_circle(int(x), int(y), fill))

    for gene in sorted(intermediates):
        if gene in pos:
            x, y = pos[gene]
            lines.append(svg_circle(int(x), int(y), INTERMEDIATE_COLOR))

    for gene in sorted(bridged_unassigned):
        if gene in pos:
            x, y = pos[gene]
            lines.append(svg_circle(int(x), int(y), BRIDGE_UNASSIGNED_COLOR))

    # ---------- then draw ALL text last ----------
    for gid in (1, 2, 3, 4):
        for gene in groups[gid]:
            x, y = pos[gene]
            lines.append(svg_text(int(x), int(y), gene, label_offsets.get(gene)))

    for gene in sorted(intermediates):
        if gene in pos:
            x, y = pos[gene]
            lines.append(svg_text(int(x), int(y), gene, label_offsets.get(gene)))

    for gene in sorted(bridged_unassigned):
        if gene in pos:
            x, y = pos[gene]
            lines.append(svg_text(int(x), int(y), gene, label_offsets.get(gene)))

    lines.append("</svg>")
    out_svg.write_text("".join(lines))
    build_group3_isolated_svg(
        pos,
        groups,
        edges,
        bridging_edges,
        intermediates,
        bridged_unassigned,
        label_offsets,
    )
    build_group3_isolated_rho_colors_svg(
        pos,
        groups,
        edges,
        bridging_edges,
        intermediates,
        bridged_unassigned,
        label_offsets,
        highlight_gene=highlight_gene,
    )


    print(f"[OK] Wrote: {out_svg}")
    print(f"[INFO] Wrote: {out_txt}")


# ----------------------- CLI --------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("genes_file", help="Text file with one gene symbol per line.")
    p.add_argument("--edges", required=True)
    p.add_argument("--highlight-gene", default="RHO", help="Center gene (default: RHO).")
    p.add_argument("--min-score", type=float, default=MIN_SCORE)
    args = p.parse_args()

    ensure_dirs()

    build_diagram(
        Path(args.edges),
        Path(args.genes_file),
        highlight_gene=args.highlight_gene,
        min_score=args.min_score,
    )


if __name__ == "__main__":
    main()