from typing import List, Optional

import pandas as pd
from pronto import Ontology
import requests
from retry import retry

# A dict of all therapeutic areas and their ids ranked by order of relevance
THERAPEUTIC_AREAS = {
    'therapeutic_area': [
        'cell proliferation disorder',
        'infectious disease',
        'pregnancy or perinatal disease',
        'animal disease',
        'disease of visual system',
        'cardiovascular disease',
        'pancreas disease',
        'gastrointestinal disease',
        'reproductive system or breast disease',
        'integumentary system disease',
        'endocrine system disease',
        'respiratory or thoracic disease',
        'urinary system disease',
        'musculoskeletal or connective tissue disease',
        'disease of ear',
        'immune system disease',
        'hematologic disease',
        'nervous system disease',
        'psychiatric disorder',
        'nutritional or metabolic disease',
        'genetic, familial or congenital disease',
        'injury, poisoning or other complication',
        'phenotype',
        'measurement',
        'biological process',
    ],
    'id': [
        'MONDO_0045024',
        'EFO_0005741',
        'OTAR_0000014',
        'EFO_0005932',
        'MONDO_0024458',
        'EFO_0000319',
        'EFO_0009605',
        'EFO_0010282',
        'OTAR_0000017',
        'EFO_0010285',
        'EFO_0001379',
        'OTAR_0000010',
        'EFO_0009690',
        'OTAR_0000006',
        'MONDO_0021205',
        'EFO_0000540',
        'EFO_0005803',
        'EFO_0000618',
        'MONDO_0002025',
        'MONDO_0024297',
        'OTAR_0000018',
        'OTAR_0000009',
        'EFO_0000651',
        'EFO_0001444',
        'GO_0008150',
    ],
}
SORTED_TAS_DF = pd.DataFrame(data=THERAPEUTIC_AREAS)

EFO_RELEASE_API_TEMPLATE = 'https://api.github.com/repos/EBISPOT/efo/releases/{}'
OWL_FILENAME = 'efo_otar_slim.owl'


@retry(tries=5, delay=3, backoff=1.5, jitter=(1, 3))
def fetch_otar_owl_from_github(efo_release):
    """Queries the GitHub API to fetch the latest EFO OTAR SLIM OWL URL."""

    if efo_release == 'latest':
        url = EFO_RELEASE_API_TEMPLATE.format(efo_release)
    else:
        url = EFO_RELEASE_API_TEMPLATE.format(f'tags/{efo_release}')
    response = requests.get(url)
    response.raise_for_status()  # In case of HTTP errors, this will be caught by the @retry decorator.

    otar_slim_assets = [asset for asset in response.json()['assets'] if asset['name'] == OWL_FILENAME]
    if len(otar_slim_assets) == 0:
        raise AssertionError(f'EFO release {efo_release!r} on GitHub does not contain the file {OWL_FILENAME!r}.')
    if len(otar_slim_assets) > 1:
        raise AssertionError(f'EFO release {efo_release!r} contains multiple files named {OWL_FILENAME!r}.')

    return otar_slim_assets[0]['browser_download_url']


def extract_therapeutic_areas_from_owl() -> pd.DataFrame:
    """
    A dataframe with all the EFO IDs and their therapeutic areas are parsed from the EFO OTAR SLIM OWL file.
    """

    owl_url = fetch_otar_owl_from_github('latest')
    efo_terms = Ontology(owl_url, timeout=10).terms()
    owl_parsed = []

    for term in efo_terms:
        # The TAs are extracted by iterating through the ancestors of a term and looking up if it's in THERAPEUTIC_AREAS
        therapeutic_areas = []
        for ancestor in term.superclasses():
            ancestor_id = normalise_ontology_identifier(ancestor.id)
            if ancestor_id in THERAPEUTIC_AREAS['id']:
                therapeutic_areas.append(ancestor_id)

        efo_id = normalise_ontology_identifier(term.id)
        owl_parsed.append((efo_id, therapeutic_areas))

    return pd.DataFrame(owl_parsed, columns=['efo_id', 'therapeutic_areas'])


def get_prioritised_therapeutic_area(
    therapeutic_areas: List,
) -> str:
    """
    SORTED_TAS_DF is a df where the therapeutic areas are arranged in order of relevance.
    The more relevant TA is extracted by selecting which one has the minimal index.
    """
    try:
        if len(therapeutic_areas) > 0:
            min_index = float('inf')
            for ta in therapeutic_areas:
                idx = SORTED_TAS_DF.index[SORTED_TAS_DF['id'] == ta]
                if idx < min_index:
                    min_index = idx
            ta = SORTED_TAS_DF.iloc[min_index]["therapeutic_area"].values[0]
            return ta
    except TypeError:
        return "Uncategorised"
    except Exception as e:
        raise e


def normalise_ontology_identifier(identifier: str) -> Optional[str]:
    """
    Normalise ontology identifier representation in order to make direct string-to-string comparison possible.
    Ex:
    'http://www.orpha.net/ORDO/Orphanet_178506' --> 'Orphanet_178506'
    'BTO:0000305' --> 'BTO_0000305'
    """

    return identifier.split('/')[-1].replace(':', '_')
