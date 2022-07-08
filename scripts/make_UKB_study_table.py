"""
Processes UK Biobank's manifest to extract all studies and their metadata in the OTG format.
"""

import argparse
import logging
from collections import OrderedDict
import re

import pandas as pd


def main(input_path: str, output_path: str) -> None:

    logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

    # Only keep required cols
    to_keep = OrderedDict(
        [("code", "study_id"), ('trait', 'trait_raw'), ("n_total", "n_total"), ("n_cases", "n_cases")]
    )

    # Load manifest
    manifest = pd.read_csv(input_path, sep="\t", header=0, dtype=object).filter(items=to_keep).rename(columns=to_keep)

    logging.info(f"{input_path} has been loaded. Formatting...")

    #
    # Add other columns -------------------------------------------------------
    #

    # Process traits to a nicer format
    manifest["trait_reported"] = manifest['trait_raw'].apply(make_trait_reported_string)

    # Vector to get Neale or SAIGE studies
    is_neale = manifest["study_id"].str.startswith("NEALE2_")
    is_saige = manifest["study_id"].str.startswith("SAIGE_")
    assert (is_neale | is_saige).all()

    # Make columns
    manifest.loc[:, "pmid"] = ""
    manifest.loc[is_neale, "pub_date"] = "2018-08-01"
    manifest.loc[is_saige, "pub_date"] = "2018-10-24"
    manifest.loc[:, "pub_journal"] = ""
    manifest.loc[:, "pub_title"] = ""
    manifest.loc[is_neale, "pub_author"] = "UKB Neale v2"
    manifest.loc[is_saige, "pub_author"] = "UKB SAIGE"
    manifest.loc[:, "n_initial"] = manifest["n_total"].apply(to_int_safe)
    manifest.loc[:, "n_cases"] = manifest["n_cases"].apply(to_int_safe)
    manifest.loc[:, "n_replication"] = 0
    manifest.loc[:, "ancestry_initial"] = "European=" + manifest["n_initial"].astype(str)
    manifest.loc[:, "ancestry_replication"] = ""

    # Ouput required columns
    cols = [
        "study_id",
        "pmid",
        "pub_date",
        "pub_journal",
        "pub_title",
        "pub_author",
        "trait_reported",
        "ancestry_initial",
        "ancestry_replication",
        "n_initial",
        "n_cases",
        "n_replication",
    ]

    manifest = manifest.filter(items=cols).drop_duplicates()

    # Assert trait_reported is not empty
    assert manifest["trait_reported"].notna().all(), "There are studies where the trait is not defined."

    # Write
    manifest.to_json(args.output, orient="records", lines=True)
    logging.info(f"{len(manifest)} studies have been saved in {output_path}. Exiting.")


def to_int_safe(i):
    try:
        return int(float(i))
    except ValueError:
        return None


def make_trait_reported_string(s_raw):
    '''Takes the raw trait name and outputs transformed name'''

    # Replace any double spaces with single
    s_raw = re.sub(r' +', r' ', s_raw)

    # Assert no "|" in trait name
    assert "|" not in s_raw, f"Reported trait ({s_raw}), contains invalid character."

    # Split prefix
    parts = s_raw.split(': ', 1)

    # Move prefix to end if exists
    if len(parts) == 2:
        trait = " | ".join([parts[1], parts[0]])
    else:
        trait = s_raw

    # Capitalise the first letter
    trait = trait.capitalize()

    return trait


def parse_args():
    """Load command line args"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        metavar="<str>",
        help=("TSV file with UK Biobank (Neale V2 and SAIGE) manifest"),
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        metavar="<str>",
        help=("JSON file with all Saige and Neale studies and their metadata."),
        type=str,
        required=True,
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # Parse args
    args = parse_args()

    main(
        input_path=args.input,
        output_path=args.output,
    )
