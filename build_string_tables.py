#!/usr/bin/env python3
"""Export STRING RP gene interactions as minimal, figure-compatible TSVs.


Usage:
    python3 -m pip install requests
    python3 generate_string_interactions_raw.py --min-score 0.40
"""

import argparse
import csv
import io
import sys
import time
from pathlib import Path

import requests

API_ROOT = "https://string-db.org/api"
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_COLUMNS = ["#node1", "node2", "combined_score"]


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
    """Rank genes by degree, breaking ties by total STRING confidence score."""
    degrees = {}
    score_sums = {}

    for row in rows:
        a, b = row["#node1"], row["node2"]
        score = float(row["combined_score"])

        for gene in (a, b):
            degrees[gene] = degrees.get(gene, 0) + 1
            score_sums[gene] = score_sums.get(gene, 0.0) + score

    ranked = sorted(
        degrees,
        key=lambda gene: (-degrees[gene], -score_sums[gene], gene)
    )

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Gene", "Degree", "Rank"])

        for rank, gene in enumerate(ranked, start=1):
            writer.writerow([gene, degrees[gene], rank])

    print(f"Wrote {path} ({len(ranked)} ranked genes)")



def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, default=SCRIPT_DIR)
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR)
    parser.add_argument("--min-score", type=float, default=0.40)
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--caller", default="rp_gene_network_export")
    args = parser.parse_args()
    if not 0 <= args.min_score <= 1:
        parser.error("--min-score must be between 0 and 1")

    master = set(read_genes(args.input_dir / "my_genes.txt"))
    groups = {i: set(read_genes(args.input_dir / f"group{i}_genes.txt")) for i in range(1, 5)}
    genes = sorted(master | set().union(*groups.values()))
    print(f"Reading {len(genes)} unique gene symbols")

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

    id_to_symbol = {protein_id: symbol for symbol, protein_id in symbol_to_id.items()}
    edges = {}
    for native in native_rows:
        converted = convert_row(native, id_to_symbol)
        if converted is not None:
            pair = (converted["#node1"], converted["node2"])
            edges.setdefault(pair, converted)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_rows = [edges[pair] for pair in sorted(edges)]
    write_tsv(args.output_dir / "interactions_all.tsv", all_rows)
    write_degree_rankings(args.output_dir / "interactions_all_degree_rankings.csv", all_rows)
    for group_number, members in groups.items():
        group_rows = [edges[pair] for pair in sorted(edges) if pair[0] in members and pair[1] in members]
        write_tsv(args.output_dir / f"interactions_g{group_number}.tsv", group_rows)
        write_degree_rankings(args.output_dir / f"interactions_g{group_number}_degree_rankings.csv", group_rows)


if __name__ == "__main__":
    main()
