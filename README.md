# KG Biolink Sanity Check

This project checks node category assignments in a Biolink knowledge graph. It
reads nodes from a newline-delimited JSON (JSONL) file and uses the installed
`biolink-model` schema to identify combinations of categories that may conflict
within the Biolink class hierarchy.

## Requirements

- Python 3.12 or later
- The project dependencies: `biolink-model` and `pyyaml`

Install the project and its dependencies from the repository root:

```powershell
python -m pip install -e .
```

## Run the category-conflict check

Place the graph's node data at `data/nodes.jsonl`, then run:

```powershell
python main.py
```

The command calls `check_category_conflicts` and writes flagged nodes to
`data/category_conflicts.jsonl`. It also prints the number of input lines,
the number of nodes with conflicts, and counts grouped by conflicting category
combination.

Each input line should be a JSON object. The check reads its `category` field
as a list of Biolink category CURIEs (or as a single category string). It
removes categories that are ancestors of another category in the same node,
then examines the remaining, most-specific categories pairwise. A pair is
reported when its classes share a nearest `is_a` ancestor below the generic
Biolink roots `NamedThing` and `Entity`. The `Gene` + `Protein` combination is
explicitly allowed and is not reported.

The Biolink hierarchy is loaded from the installed `biolink-model` package;
the version is pinned in `pyproject.toml`.

## Optional required-category filter

`category_checks.py` also contains `filter_and_export_nodes`, which looks for
nodes whose `category` field contains both:

- `biolink:MolecularMixture`
- `biolink:SmallMolecule`

It exports every node containing both categories, including nodes with other
categories, to `problematic_biolink_nodes.jsonl` by default. It separately
counts matches whose raw `category` list has exactly two entries, and nodes
with no categories. This function is not called by `main.py`. To run it with
its default paths:

```powershell
python category_checks.py
```

The filter's “exactly two” count is based on the raw node category-list length.
It is therefore not necessarily the same as a graph-metadata count for a
category combination: node records may include additional ancestor or related
Biolink categories alongside the two categories of interest.