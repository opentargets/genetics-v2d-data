import pandas as pd
import numpy as np
import json

assert(pd.__version__ >= '0.24')

rule make_gwas_cat_studies_table:
    ''' Make GWAS Catalog table. It needs to use the intermediate top loci table
        as we must create new study IDs based on the subphenotype in the field
        "P-VALUE (TEXT)". This isn't available in the GWAS Catalog study file.
    '''
    input:
        toploci=tmpdir + '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv',
        ancestries=HTTPRemoteProvider().remote('https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry', keep_local=KEEP_LOCAL)
    output:
        main=tmpdir + '/{version}/gwas-catalog_study_table.tsv',
        lut=tmpdir + '/{version}/gwas-catalog_study_id_lut.tsv'
    shell:
        'python scripts/make_gwas_cat_study_table.py '
        '--in_toploci {input.toploci} '
        '--in_ancestries {input.ancestries} '
        '--outf {output.main} '
        '--out_id_lut {output.lut}'

rule process_efo_curation:
    ''' Downloads EFO curation, applies filters, merges
    Finally, maps EFOs to mapped trait names using ? API.
    '''
    input:
        icd10=GSRemoteProvider().remote(config['neale_efo_icd10'], keep_local=KEEP_LOCAL),
        selfrep=GSRemoteProvider().remote(config['neale_efo_self'], keep_local=KEEP_LOCAL)
    output:
        tmpdir + '/{version}/nealeUKB_efo_curation.tsv'
    shell:
        'python scripts/process_nealeUKB_efo_curations.py '
        '--in_icd10 "{input.icd10}" '
        '--in_self "{input.selfrep}" '
        '--outf {output}'

rule make_nealeUKB_studies_table:
    ''' Makes study table for Neale et al UKB summary statistics
    '''
    input:
        manifest=GSRemoteProvider().remote(config['neale_manifest'], keep_local=KEEP_LOCAL),
        categories=GSRemoteProvider().remote(config['neale_categories'], keep_local=KEEP_LOCAL),
        efos=tmpdir + '/{version}/nealeUKB_efo_curation.tsv'
    output:
        tmpdir + '/{version}/nealeUKB_study_table.tsv'
    shell:
        'python scripts/make_nealeUKB_study_table.py '
        '--in_manifest {input.manifest} '
        '--in_efos {input.efos} '
        '--in_categories {input.categories} '
        '--outf {output}'

rule merge_study_tables:
    ''' Merges the GWAS Catalog and Neale UK Biobank study tables together
    '''
    input:
        gwas=tmpdir + '/{version}/gwas-catalog_study_table.tsv',
        neale=tmpdir + '/{version}/nealeUKB_study_table.tsv',
        top_loci='output/{version}/toploci.parquet'
    output:
        'output/{version}/studies.parquet'
    shell:
        'python scripts/merge_study_tables.py '
        '--in_gwascat {input.gwas} '
        '--in_neale {input.neale} '
        '--in_toploci {input.top_loci} '
        '--output {output}'

rule study_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/studies.parquet'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/studies.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
