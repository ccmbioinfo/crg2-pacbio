#!/usr/bin/env python3
"""Add indirect HPO matches to the setup-generated patient HPO file."""

import logging
import re
from pathlib import Path

import pandas as pd

from hpo_indirect_matcher import get_indirect_hpo_gene_matches, load_references


HPO_ID_PATTERN = re.compile(r"HP:\d{7}")


def get_patient_hpo_ids(hpo_file):
    """Read patient terms from the comment and legacy multi-term rows."""
    hpo_ids = []
    with open(hpo_file, encoding="ISO-8859-1") as handle:
        for line in handle:
            if line.lstrip().startswith("#") and "patient_hpo_ids" in line.lower():
                hpo_ids.extend(HPO_ID_PATTERN.findall(line))

    hpo = pd.read_csv(hpo_file, sep="\t", comment="#", dtype=str)
    hpo.columns = hpo.columns.str.strip()
    hpo = hpo.rename(columns={"Gene symbol": "Gene Symbol"})
    for value in hpo["HPO IDs"].dropna():
        hpo_ids.extend(HPO_ID_PATTERN.findall(value))

    return hpo, list(dict.fromkeys(hpo_ids))


def format_score(score):
    score = f"{float(score):.3f}".rstrip("0").rstrip(".")
    return score if "." in score else f"{score}.0"


def join_unique(values):
    return "; ".join(dict.fromkeys(str(value) for value in values))


def prepare_hpo_matches(hpo_file, hpo_dir, output_file):
    """Combine exact and indirect matches in the original one-row-per-gene format."""
    hpo, hpo_ids = get_patient_hpo_ids(hpo_file)
    hpo_mapping = pd.read_csv(f"{hpo_dir}/genes_to_phenotype.txt", sep="\t")
    hpo_data = load_references(f"{hpo_dir}/hp.json", hpo_mapping)
    hpo_ids = [hpo_data["alternative_ids"].get(hpo_id, hpo_id) for hpo_id in hpo_ids]

    exact_hpo = hpo_mapping[hpo_mapping["hpo_id"].isin(hpo_ids)][
        ["gene_symbol", "hpo_id", "hpo_name"]
    ].copy()
    exact_hpo["HPO Match Score"] = 1.0

    indirect_hpo = get_indirect_hpo_gene_matches(hpo_ids, hpo_data)
    hpo_matches = pd.concat([exact_hpo, indirect_hpo], ignore_index=True)
    hpo_matches = hpo_matches.sort_values("HPO Match Score", ascending=False)
    hpo_matches = hpo_matches.drop_duplicates(["gene_symbol", "hpo_id"])

    # Reuse setup-resolved IDs for direct genes and the standard HGNC map for
    # genes introduced by indirect matching.
    hgnc = pd.read_csv(f"{hpo_dir}/HGNC_ensembl_map.csv", dtype=str)
    gene_ids = dict(zip(hgnc["hgnc_symbol"], hgnc["ensembl_gene_id"]))
    gene_ids.update(dict(zip(hpo["Gene Symbol"], hpo["Gene ID"])))
    hpo_matches["Gene ID"] = hpo_matches["gene_symbol"].map(gene_ids)

    hpo_matches["Features"] = hpo_matches.apply(
        lambda row: f'{row["hpo_name"]} ({format_score(row["HPO Match Score"])})',
        axis=1,
    )
    hpo_matches["HPO Match Score"] = hpo_matches["HPO Match Score"].map(format_score)

    hpo_matches = hpo_matches.groupby(
        ["gene_symbol", "Gene ID"], as_index=False, dropna=False
    ).agg(
        **{
            "Number of occurrences": ("hpo_id", "nunique"),
            "Features": ("Features", join_unique),
            "HPO IDs": ("hpo_id", join_unique),
            "HPO Match Score": (
                "HPO Match Score",
                lambda scores: "; ".join(scores),
            ),
        }
    )
    hpo_matches = hpo_matches.rename(columns={"gene_symbol": "Gene Symbol"})

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    hpo_matches.to_csv(output_file, sep="\t", index=False)
    logging.info("Wrote HPO matches for %d genes to %s", len(hpo_matches), output_file)


if "snakemake" in globals():
    Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=str(snakemake.log[0]), level=logging.INFO)
    prepare_hpo_matches(
        snakemake.input.hpo,
        snakemake.params.hpo_dir,
        snakemake.output.hpo,
    )
