#!/usr/bin/env python3
"""Export STRING RP gene interactions as minimal, figure-compatible TSVs.

Usage:
    python3 -m pip install requests
    python3 generate_string_interactions_full.py --min-score 0.40
"""

import argparse
import csv
import io
import json
import sys
import time
from pathlib import Path

import requests

API_ROOT = "https://string-db.org/api"
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_COLUMNS = ["#node1", "node2", "combined_score"]

SEED_GROUPS = {
    1: {"RPE65", "LRAT", "MERTK", "RBP3", "RGR"},
    2: {"RHO", "PRPH2", "SAG", "CNGA1", "CNGB1", "FSCN2", "ROM1", "IMPG2", "PROM1"},
    3: {"RPGR", "RP1", "CLRN1", "MAK", "USH2A"},
    4: {"NR2E3", "CRX", "ZNF513", "PRPF31", "PRPF8", "PRPF3", "PRPF4", "PRPF6", "RP9", "SNRNP200", "DHX38"},
}


def read_genes(path):
    """Read unique gene symbols, one per line."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing gene list: {path}")
    return list(dict.fromkeys(
        symbol for line in path.read_text(encoding="utf-8-sig").splitlines()
        if (symbol := line.strip().strip("'\" ,").upper())
    ))


def request_string(session, endpoint, params, response_type):
    """Request STRING data, retrying transient network errors."""
    url = f"{API_ROOT}/{response_type}/{endpoint}"
    for attempt in range(4):
        try:
            response = session.post(url, data=params, timeout=180)
            response.raise_for_status()
            return response.json() if response_type == "json" else response.text
        except (requests.RequestException, ValueError) as exc:
            if attempt == 3:
                raise RuntimeError(f"STRING request failed ({endpoint}): {exc}") from exc
            time.sleep(2 ** attempt)


def resolve_ids(session, genes, species, caller):
    """Resolve gene symbols to STRING IDs for consistent output naming."""
    mapping = {}
    for start in range(0, len(genes), 100):
        batch = genes[start:start + 100]
        records = request_string(session, "get_string_ids", {
            "identifiers": "\r".join(batch),
            "species": species,
            "limit": 1,
            "caller_identity": caller,
        }, "json")
        for record in records:
            symbol = str(record.get("queryItem", "")).upper()
            protein_id = record.get("stringId")
            if symbol in batch and protein_id and symbol not in mapping:
                mapping[symbol] = protein_id
    return mapping


def fetch_network(session, protein_ids, species, min_score, caller):
    """Retrieve functional associations without additional STRING nodes."""
    response = request_string(session, "network", {
        "identifiers": "\r".join(protein_ids),
        "species": species,
        "required_score": round(min_score * 1000),
        "network_type": "functional",
        "add_nodes": 0,
        "caller_identity": caller,
    }, "tsv")
    reader = csv.DictReader(io.StringIO(response), delimiter="\t")
    required = {"stringId_A", "stringId_B", "score"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise RuntimeError(f"Unexpected STRING columns: {reader.fieldnames}")
    rows = list(reader)
    if any(None in row for row in rows):
        raise RuntimeError("Malformed STRING TSV response")
    return rows


def convert_row(row, id_to_symbol):
    """Convert a STRING association to the three figure-required fields."""
    gene_a = id_to_symbol.get(row["stringId_A"])
    gene_b = id_to_symbol.get(row["stringId_B"])
    if not gene_a or not gene_b or gene_a == gene_b:
        return None
    gene_a, gene_b = sorted((gene_a, gene_b))
    return {"#node1": gene_a, "node2": gene_b, "combined_score": row["score"]}


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path} ({len(rows)} interactions)")


def write_degree_rankings(path, rows):
    """Rank by degree, then summed STRING confidence, then gene name."""
    degrees = {}
    score_sums = {}
    for row in rows:
        a, b = row["#node1"], row["node2"]
        score = float(row["combined_score"])
        for gene in (a, b):
            degrees[gene] = degrees.get(gene, 0) + 1
            score_sums[gene] = score_sums.get(gene, 0.0) + score

    ranked = sorted(degrees, key=lambda gene: (-degrees[gene], -score_sums[gene], gene))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Gene", "Degree", "Rank"])
        for rank, gene in enumerate(ranked, start=1):
            writer.writerow([gene, degrees[gene], rank])
    print(f"Wrote {path} ({len(ranked)} ranked genes)")


def fetch_grouping_neighbors(session, gene, species, caller, cache):
    """Match rp_groups.py: query each gene's top 500 STRING neighbors."""
    key = f"NEI::{gene}"
    if key in cache:
        return cache[key]

    records = request_string(session, "get_string_ids", {
        "identifiers": gene, "species": species, "limit": 1,
        "caller_identity": caller,
    }, "json")
    if not records or not records[0].get("stringId"):
        result = {"preferred": None, "edges": []}
    else:
        record = records[0]
        preferred = (record.get("preferredName") or gene).upper()
        edges = request_string(session, "network", {
            "identifiers": record["stringId"], "species": species,
            "limit": 500, "caller_identity": caller,
        }, "json")
        result = {"preferred": preferred, "edges": edges}
    cache[key] = result
    return result


def build_grouping_graph(session, genes, species, min_score, caller, cache):
    """Reproduce rp_groups.py's directed adjacency and maximum-score updates."""
    universe = set(genes)
    graph = {}
    for gene in universe:
        info = fetch_grouping_neighbors(session, gene, species, caller, cache)
        preferred = info["preferred"]
        if not preferred:
            continue
        for edge in info["edges"]:
            a = (edge.get("preferredName_A") or "").upper()
            b = (edge.get("preferredName_B") or "").upper()
            score = float(edge.get("score") or 0.0)
            if score < min_score:
                continue
            other = b if a == preferred else a if b == preferred else None
            if other in universe:
                graph.setdefault(preferred, {})[other] = max(
                    graph.get(preferred, {}).get(other, 0.0), score
                )
    return graph


def assign_groups(genes, graph):
    """Replicate rp_groups.py: alphabetical passes, forward AND reverse votes."""
    groups = {number: set(seeds) for number, seeds in SEED_GROUPS.items()}
    assigned = set().union(*groups.values())
    remaining = [gene for gene in genes if gene not in assigned]

    changed = True
    while changed and remaining:
        changed = False
        still_unplaced = []
        for gene in remaining:
            counts = {number: 0 for number in groups}
            sums = {number: 0.0 for number in groups}
            for other, score in graph.get(gene, {}).items():
                for number, members in groups.items():
                    if other in members:
                        counts[number] += 1
                        sums[number] += score
            for source, neighbors in graph.items():
                if gene in neighbors:
                    for number, members in groups.items():
                        if source in members:
                            counts[number] += 1
                            sums[number] += neighbors[gene]
            best = max((counts[number], sums[number], number) for number in groups)
            if best[0] > 0:
                groups[best[2]].add(gene)
                changed = True
            else:
                still_unplaced.append(gene)
        remaining = still_unplaced
    return groups, set(remaining)


def write_gene_list(path, genes):
    path.write_text("".join(f"{gene}\n" for gene in sorted(genes)), encoding="utf-8")
    print(f"Wrote {path} ({len(genes)} genes)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, default=SCRIPT_DIR)
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR / "string_outputs")
    parser.add_argument("--min-score", type=float, default=0.40)
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--caller", default="rp_gene_network_export")
    args = parser.parse_args()
    if not 0 <= args.min_score <= 1:
        parser.error("--min-score must be between 0 and 1")

    user_genes = read_genes(args.input_dir / "my_genes.txt")
    # The original rp_groups.py added every seed, then sorted the full universe.
    genes = sorted(set(user_genes).union(*SEED_GROUPS.values()))
    print(f"Reading {len(user_genes)} unique input genes; {len(genes)} including seeds")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with requests.Session() as session:
        symbol_to_id = resolve_ids(session, genes, args.species, args.caller)
        unresolved = sorted(set(genes) - symbol_to_id.keys())
        if unresolved:
            print(f"WARNING: Unresolved genes: {', '.join(unresolved)}", file=sys.stderr)
        if len(set(symbol_to_id.values())) != len(symbol_to_id):
            raise RuntimeError("Multiple gene symbols resolve to the same STRING ID")
        if not symbol_to_id:
            raise RuntimeError("No gene symbols resolved in STRING")
        native_rows = fetch_network(session, list(symbol_to_id.values()), args.species, args.min_score, args.caller)

        # Grouping intentionally uses the original per-gene STRING requests,
        # independent of the network-wide request used for TSV/CSV exports.
        cache_path = args.output_dir / "grouping_string_cache.json"
        cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
        grouping_graph = build_grouping_graph(
            session, genes, args.species, args.min_score, "rp_mapping", cache
        )

    id_to_symbol = {protein_id: symbol for symbol, protein_id in symbol_to_id.items()}
    edges = {}
    for native in native_rows:
        converted = convert_row(native, id_to_symbol)
        if converted is not None:
            pair = (converted["#node1"], converted["node2"])
            edges.setdefault(pair, converted)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_rows = [edges[pair] for pair in sorted(edges)]
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")
    groups, unassigned = assign_groups(genes, grouping_graph)
    for group_number, members in groups.items():
        write_gene_list(args.output_dir / f"group{group_number}_genes.txt", members)
    write_gene_list(args.output_dir / "unassigned_genes.txt", unassigned)
    write_tsv(args.output_dir / "interactions_all.tsv", all_rows)
    write_degree_rankings(args.output_dir / "interactions_all_degree_rankings.csv", all_rows)
    for group_number, members in groups.items():
        group_rows = [edges[pair] for pair in sorted(edges) if pair[0] in members and pair[1] in members]
        write_tsv(args.output_dir / f"interactions_g{group_number}.tsv", group_rows)
        write_degree_rankings(args.output_dir / f"interactions_g{group_number}_degree_rankings.csv", group_rows)


if __name__ == "__main__":
    main()
