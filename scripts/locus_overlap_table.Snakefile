rule calculate_overlaps:
    ''' Calcs overlap between trait associated loci
    '''
    input:
        top_loci='output/{version}/toploci.tsv',
        ld='output/{version}/ld.tsv.gz',
        finemap='output/{version}/finemapping.tsv.gz'
    output:
        'output/{version}/locus_overlap.tsv.gz'
    shell:
        'pypy3 scripts/calculate_locus_set_overlaps.py '
        '--top_loci {input.top_loci} '
        '--ld {input.ld} '
        '--finemap {input.finemap} '
        '--outf {output}'

rule overlap_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/locus_overlap.tsv.gz'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/locus_overlap.tsv.gz'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
