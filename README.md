# RP-Int-Mapping

Code and input data for mapping protein–protein interactions among retinitis pigmentosa (RP) genes using the [STRING](https://string-db.org/) database. The scripts build the interaction tables, rank genes by network degree, and draw the network figures used in the manuscript.

## Contents

| Path | Purpose |
|---|---|
| `build_string_tables.py` | Queries STRING for all genes and writes the interaction tables (`interactions_*.tsv`) and degree-ranking files (`*_degree_rankings.csv`). |
| `plot_group_network.py` | Draws the main group network figure (centered on one gene, default RHO), the unassigned-gene → intermediate → assigned-gene bridges, and two close-ups of Group 3. |
| `plot_degree_layers.py` | Draws the degree-ranked layered network (default: Group 3). |
| `optional_regenerate_groups.py` | Optional. Re-derives the four gene groups from scratch (see [Reproducibility notes](#reproducibility-notes)). |
| `original_input_files/` | The original input files used to produce the paper figures (listed below). |
| `paper_figures/` | The figures as they appear in the manuscript, for comparison with figures you regenerate. |

### Files in `original_input_files/`

- `my_genes.txt`: master list of RP genes, one gene symbol per line.
- `group1_genes.txt` … `group4_genes.txt`: the four gene groups used in the paper (Group 1 – RPE, Group 2 – Outer Segment, Group 3 – Connecting Cilium, Group 4 – Nucleus).
- `interactions_all.tsv`, `interactions_g1.tsv` … `interactions_g4.tsv`: STRING interactions among all genes / within each group (columns: `#node1`, `node2`, `combined_score`).
- `interactions_*_degree_rankings.csv`: genes ranked by degree (columns: `Gene`, `Degree`, `Rank`).
- `info_files/string_cache.json`: saved STRING responses used to find intermediates. [Included here / download link — edit.]

## Requirements

Tested with Python 3.11.

```
python -m pip install -r requirements.txt
```

(`requests`, `pandas`, `networkx`, `matplotlib`)

## Reproducing the figures

**Setup.** The scripts read and write files relative to the folder you run them from, so first copy the input files next to the scripts (run from the repository folder):

```
# Windows (PowerShell)
Copy-Item original_input_files\* . -Recurse

# macOS / Linux
cp -R original_input_files/* .
```

This leaves `original_input_files/` untouched. The STRING cache must end up at `info_files/string_cache.json`. Generated outputs (`svg_files/`, `info_files/`, `layered_network_outputs/`) are created in this working folder. Run all commands below from it.

**1. Interaction tables (optional).** The tables are included in `original_input_files/`. To regenerate them from STRING, which needs internet access:

```
python build_string_tables.py --min-score 0.40
```

This reads `my_genes.txt` and `group1–4_genes.txt` and writes `interactions_all.tsv`, `interactions_g1–g4.tsv` and the degree-ranking CSVs. It overwrites the copies in your working folder, not the originals.

**2. Main network figure.**

```
python plot_group_network.py my_genes.txt --edges interactions_all.tsv --min-score 0.40
```

Optional: `--highlight-gene GENE` sets the central gene (default `RHO`; it must belong to one of the four groups).

Needs: `my_genes.txt`, `group1–4_genes.txt`, `interactions_all.tsv`, `interactions_g3_degree_rankings.csv`.

Writes:
- `svg_files/rp_with_intermediates_new_highlight_RHO.svg` (main figure)
- `svg_files/group3_isolated_intermediates_degree_colors.svg` and `svg_files/group3_isolated_intermediates_RHO_colors.svg` (Group 3 close-ups)
- `info_files/unassigned_with_intermediates.txt` (table of unassigned genes, their intermediates, and the STRING scores)

**3. Degree-layered network figure.**

```
python plot_degree_layers.py
```

Needs: `interactions_g3.tsv` and `interactions_g3_degree_rankings.csv`. Writes `layered_network_outputs/new_interactions_g3_degree_layers_2_3_4_5_rest.svg`. The network (`NETWORK`) and layer sizes (`LAYERS`) are set at the top of the script.

**4. Compare with the paper.** Open the regenerated SVGs next to the figures in `paper_figures/`. Exact agreement requires the same inputs and the same STRING cache (see below).

## Methods notes

- **Interactions:** STRING functional associations for human (species 9606), combined score ≥ 0.40.
- **Degree ranking:** genes are ranked by number of interactions, ties broken by the sum of their STRING scores.
- **Intermediates:** an *unassigned* gene is a gene in `my_genes.txt` that is not in any of the four groups. For each one, the script looks for a two-step path (unassigned → intermediate → assigned gene) through a STRING neighbor that is not itself in a group. It keeps the path with the highest summed score; both links must be ≥ 0.40. STRING scores are undirected.

## Reproducibility notes

- **STRING changes over time.** The intermediates come from live STRING queries, so results can differ between releases. To reproduce the paper's figure exactly, use the provided `info_files/string_cache.json` (copied into the working folder as in Setup): the script reads it and does not contact STRING for genes it already contains. The script only rewrites the cache when it has to fetch new genes. STRING version used for the paper: [version / date — add].
- **Gene groups.** The four group files used in the paper are in `original_input_files/` and are the intended inputs. `optional_regenerate_groups.py` re-derives the groups from STRING with the original grouping logic, but the assignment is order-dependent and depends on the STRING release, so its output can differ from the published groups. Use it to explore, not to reproduce the figures.
- **Layout.** Node positions and label placement are deterministic for the same inputs.
