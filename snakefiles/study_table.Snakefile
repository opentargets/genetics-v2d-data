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

# rule process_efo_curation:
#     ''' Downloads EFO curation, applies filters, merges
#     Finally, maps EFOs to mapped trait names using ? API.
#     '''
#     input:
#         icd10=GSRemoteProvider().remote(config['neale_efo_icd10'], keep_local=KEEP_LOCAL),
#         selfrep=GSRemoteProvider().remote(config['neale_efo_self'], keep_local=KEEP_LOCAL)
#     output:
#         tmpdir + '/{version}/nealeUKB_efo_curation.tsv'
#     shell:
#         'python scripts/process_nealeUKB_efo_curations.py '
#         '--in_icd10 "{input.icd10}" '
#         '--in_self "{input.selfrep}" '
#         '--outf {output}'

rule make_UKB_studies_table:
    ''' Makes study table for UKB summary statistics (Neale v2 and SAIGE)
    '''
    input:
        manifest=GSRemoteProvider().remote(config['ukb_manifest'], keep_local=KEEP_LOCAL),
        # efos=GSRemoteProvider().remote(config['ukb_efo_curation'], keep_local=KEEP_LOCAL)
    output:
        study_table=tmpdir + '/{version}/UKB_study_table.tsv',
        prefix_counts=tmpdir + '/{version}/UKB_prefix_counts.tsv'
    shell:
        'python scripts/make_UKB_study_table.py '
        '--in_manifest {input.manifest} '
        # '--in_efos {input.efos} '
        '--prefix_counts {output.prefix_counts} '
        '--outf {output.study_table}'

rule merge_study_tables:
    ''' Merges the GWAS Catalog and Neale UK Biobank study tables together
    '''
    input:
        gwas=tmpdir + '/{version}/gwas-catalog_study_table.tsv',
        ukb=tmpdir + '/{version}/UKB_study_table.tsv'
    output:
        tmpdir + '/{version}/merged_study_table.tsv'
    shell:
        'python scripts/merge_study_tables.py '
        '--in_gwascat {input.gwas} '
        '--in_ukb {input.ukb} '
        '--output {output}'

rule list_studies_with_sumstats:
    ''' Makes a list of files with sumstats
    '''
    params:
        url=config['gcs_sumstat_pattern']
    output:
        tmpdir + '/{version}/studies_with_sumstats.tsv'
    shell:
        'gsutil -m ls -d "{params.url}" > {output}'

rule study_table_to_parquet:
    ''' Converts study table to final parquet
    '''
    input:
        study_table = tmpdir + '/{version}/merged_study_table.tsv',
        sumstat_studies = tmpdir + '/{version}/studies_with_sumstats.tsv',
        top_loci='output/{version}/toploci.parquet'
    output:
        'output/{version}/studies.parquet'
    shell:
        'python scripts/study_table_to_parquet.py '
        '--in_study_table {input.study_table} '
        '--in_toploci {input.top_loci} '
        '--sumstat_studies {input.sumstat_studies} '
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
