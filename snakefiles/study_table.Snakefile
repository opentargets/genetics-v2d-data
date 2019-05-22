import pandas as pd
import numpy as np
import json

assert(pd.__version__ >= '0.24')

rule make_gwas_cat_studies_table:
    ''' Make GWAS Catalog table.

        We need to join GWAS Catalog's study table with the top loci
        (association) table in order to get subphenotype in the
        "P-VALUE (TEXT)" field.

        We can't use the top loci table by itself as this would only include
        studies with >0 associations (may not be valid for phewas studies).
    '''
    input:
        gwas_study = HTTPRemoteProvider().remote(
            'https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative', keep_local=KEEP_LOCAL),
        toploci=tmpdir + '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv',
        ancestries=HTTPRemoteProvider().remote(
            'https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry', keep_local=KEEP_LOCAL)
    output:
        main=tmpdir + '/{version}/gwas-catalog_study_table.json',
        lut=tmpdir + '/{version}/gwas-catalog_study_id_lut.tsv'
    shell:
        'python scripts/make_gwas_cat_study_table.py '
        '--in_gwascat_study {input.gwas_study} '
        '--in_toploci {input.toploci} '
        '--in_ancestries {input.ancestries} '
        '--outf {output.main} '
        '--out_id_lut {output.lut}'

rule make_UKB_studies_table:
    ''' Makes study table for UKB summary statistics (Neale v2 and SAIGE)
    '''
    input:
        manifest=GSRemoteProvider().remote(config['ukb_manifest'], keep_local=KEEP_LOCAL),
        efos=GSRemoteProvider().remote(config['ukb_efo_curation'], keep_local=KEEP_LOCAL)
    output:
        study_table = tmpdir + '/{version}/UKB_study_table.json',
        prefix_counts=tmpdir + '/{version}/UKB_prefix_counts.tsv'
    shell:
        'python scripts/make_UKB_study_table.py '
        '--in_manifest {input.manifest} '
        '--in_efos {input.efos} '
        '--prefix_counts {output.prefix_counts} '
        '--outf {output.study_table}'

rule merge_study_tables:
    ''' Merges the GWAS Catalog and Neale UK Biobank study tables together
    '''
    input:
        gwas = rules.make_gwas_cat_studies_table.output.main,
        ukb = rules.make_UKB_studies_table.output.study_table
    output:
        tmpdir + '/{version}/merged_study_table.json'
    shell:
        'python scripts/merge_study_tables.py '
        '--in_gwascat {input.gwas} '
        '--in_ukb {input.ukb} '
        '--output {output}'

rule get_efo_categories:
    ''' Uses OLS API to get "therapeutic area" / category for each EFO
    '''
    input:
        study_table = rules.merge_study_tables.output,
        therapeutic_areas = config['efo_therapeutic_areas']
    output:
        tmpdir + '/{version}/efo_categories.json'
    shell:
        'python scripts/get_therapeutic_areas.py '
        '--in_study {input.study_table} '
        '--in_ta {input.therapeutic_areas} '
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
        study_table = rules.merge_study_tables.output,
        sumstat_studies = rules.list_studies_with_sumstats.output,
        efo_anno = rules.get_efo_categories.output,
        top_loci = rules.merge_gwascat_and_sumstat_toploci.output,
        therapeutic_areas = config['efo_therapeutic_areas']
    output:
        'output/{version}/studies.parquet'
    params:
        unknown_label = config['therapeutic_area_unknown_label']
    shell:
        'python scripts/study_table_to_parquet.py '
        '--in_study_table {input.study_table} '
        '--in_toploci {input.top_loci} '
        '--sumstat_studies {input.sumstat_studies} '
        '--efo_categories {input.efo_anno} '
        '--in_ta {input.therapeutic_areas} '
        '--unknown_label {params.unknown_label} '
        '--output {output}'

rule study_to_GCS:
    ''' Copy to GCS
    '''
    input:
        rules.study_table_to_parquet.output
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/studies.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
