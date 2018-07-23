rule calculate_overlaps:
    ''' Calcs overlap between trait associated loci
    '''
    input:
        top_loci='output/ot_genetics_toploci_table.{version}.tsv',
        ld='output/ot_genetics_ld_table.{version}.tsv.gz',
        finemap='output/ot_genetics_finemapping_table.{version}.tsv.gz'
    output:
        'output/ot_genetics_locus_overlap_table.{version}.tsv.gz'
    shell:
        'python scripts/calculate_locus_set_overlaps.py '
        '--top_loci {input.top_loci} '
        '--ld {input.ld} '
        '--finemap {input.finemap} '
        '--outf {output}'



# rule ld_to_GCS:
#     ''' Copy to GCS
#     '''
#     input:
#         'output/ot_genetics_ld_table.{version}.tsv.gz'
#     output:
#         GSRemoteProvider().remote(
#             '{gs_dir}/{{version}}/ot_genetics_ld_table.{{version}}.tsv.gz'.format(gs_dir=config['gs_dir'])
#             )
#     shell:
#         'cp {input} {output}'
