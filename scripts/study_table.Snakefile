import pandas as pd

# rule make_gwas_cat_studies_table:
#     ''' Use GWAS Catalog API to get all studies
#     '''
#     output:
#         tmpdir + '/{version}/gwas-catalog_study_table.tsv'
#     shell:
#         'python scripts/make_gwas_cat_study_table.py '
#         '--outf {output}'

rule make_gwas_cat_studies_table:
    ''' Make GWAS Catalog table
    '''
    input:
        studies=HTTPRemoteProvider().remote('https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative', keep_local=True),
        ancestries=HTTPRemoteProvider().remote('https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry', keep_local=True)
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
        icd10=GSRemoteProvider().remote(config['neale_efo_icd10'], keep_local=True),
        selfrep=GSRemoteProvider().remote(config['neale_efo_self'], keep_local=True)
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
        manifest=GSRemoteProvider().remote(config['neale_manifest'], keep_local=True),
        efos=tmpdir + '/{version}/nealeUKB_efo_curation.tsv'
    output:
        tmpdir + '/{version}/nealeUKB_study_table.tsv'
    shell:
        'python scripts/make_nealeUKB_study_table.py '
        '--in_manifest {input.manifest} '
        '--in_efos {input.efos} '
        '--outf {output}'

rule merge_study_tables:
    ''' Merges the GWAS Catalog and Neale UK Biobank study tables together
    '''
    input:
        gwas=tmpdir + '/{version}/gwas-catalog_study_table.tsv',
        neale=tmpdir + '/{version}/nealeUKB_study_table.tsv'
    output:
        'output/{version}/studies.tsv'
    run:
        # Load
        gwas = pd.read_csv(input['gwas'], sep='\t', header=0)
        neale = pd.read_csv(input['neale'], sep='\t', header=0)
        # Merge
        merged = pd.concat([gwas, neale], sort=False)
        # Save
        merged.to_csv(output[0], sep='\t', index=None)

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
