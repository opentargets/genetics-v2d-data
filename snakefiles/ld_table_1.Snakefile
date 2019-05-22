#
# Rules to make input manifest for LD calculation ------------------------------
#


rule make_ld_input_queries:
    ''' Make query variant and superpopulation proportions table. This table
        will be the input for making our improved look-up table
    '''
    input:
        loci = rules.merge_gwascat_and_sumstat_toploci.output,
        study = rules.study_table_to_parquet.output,
        pop_map=config['gwascat_2_superpop']
    output:
        'output/{version}/ld_analysis_input.tsv.gz'
    shell:
        'python scripts/create_ld_input_table.py '
        '--in_loci {input.loci} '
        '--in_study {input.study} '
        '--in_popmap {input.pop_map} '
        '--outf {output}'

rule ld_input_to_GCS:
    ''' Copy to GCS
    '''
    input:
        rules.make_ld_input_queries.output
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/extras/ld_analysis_input.tsv.gz'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
