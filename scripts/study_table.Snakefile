import pandas as pd
import numpy as np
import json

rule make_gwas_cat_studies_table:
    ''' Make GWAS Catalog table
    '''
    input:
        studies=HTTPRemoteProvider().remote('https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative', keep_local=KEEP_LOCAL),
        ancestries=HTTPRemoteProvider().remote('https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry', keep_local=KEEP_LOCAL)
    output:
        tmpdir + '/{version}/gwas-catalog_study_table.tsv'
    shell:
        'python scripts/make_gwas_cat_study_table.py '
        '--in_studies {input.studies} '
        '--in_ancestries {input.ancestries} '
        '--outf {output}'

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
        top_loci='output/{version}/toploci.tsv'
    output:
        'output/{version}/studies.tsv'
    run:

        # Load
        # casting everything as object since, pandas does not support NaN for int
        # see: http://pandas.pydata.org/pandas-docs/stable/gotchas.html#support-for-integer-na
        gwas = pd.read_csv(input['gwas'], sep='\t', header=0, dtype=object)
        neale = pd.read_csv(input['neale'], sep='\t', header=0, dtype=object)
        # Merge gwas cat and neale
        merged = pd.concat([gwas, neale], sort=False)
        # Load the number of associated loci per study
        num_loci = ( pd.read_csv(input['top_loci'], sep='\t', header=0)
                       .groupby('study_id')
                       .size().reset_index(name='counts')
                       .rename(columns={'counts': 'num_assoc_loci'}) )
        merged = pd.merge(merged, num_loci, how='left', on='study_id')
        merged['num_assoc_loci'] = merged['num_assoc_loci'].fillna(value=0).astype(int)
        # Save
        merged.to_csv(output[0], sep='\t', index=None)

rule make_json_study_table:
    ''' Transform into json
    '''
    input:
        'output/{version}/studies.tsv'
    output:
        'output/{version}/studies.json'
    shell:
        'python scripts/study_tsv2json.py {input} {output} '

rule study_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/studies.tsv'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/studies.tsv'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'

rule studyjson_to_GCS_json:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/studies.json'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/studies.json'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
