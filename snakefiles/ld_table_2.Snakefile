
#
# Get 1000 Genomes haplotypes from GCS -----------------------------------------
#

hap1000G_pops = ['EUR', 'EAS', 'AMR', 'AFR', 'SAS']
hap1000G_chroms = list(range(1, 23)) + ['X']
# hap1000G_chroms = ['22']

rule get_1000G_from_GCS:
    ''' Copy 1000G plink files from GCS
    '''
    output:
        expand(tmpdir + '/{version}/ld/1000Genomep3/{pop}/{pop}.{chrom}.1000Gp3.20130502.{ext}',
               pop=hap1000G_pops,
               chrom=hap1000G_chroms,
               ext=['bed', 'bim', 'fam'],
               version=config['version'])
    params:
        url_1000G = config['url_1000G'],
        outdir=tmpdir + '/{version}/ld/1000Genomep3'.format(version=config['version'])
    shell:
        'gsutil -m rsync -r {params.url_1000G} {params.outdir}'

#
# Calculate LD using plink -----------------------------------------------------
#

rule write_variant_list:
    ''' Write list variants to a file. Variants must be in the same format as
        in the bim files chrom:pos:ref:alt
    '''
    input:
        in_manifest
    output:
        tmpdir + '/{version}/ld/variant_list.txt'.format(version=config['version'])
    params:
        varlist = varid_list
    run:
        print(params['varlist'])
        with open(output[0], 'w') as out_h:
            for var in params['varlist']:
                out_h.write(
                    var.replace('_', ':') + '\n'
                )

    # shell:
    #     "zcat < {input} | tail -n +2 | cut -f 2 | sed 's/_/:/g' | sort | uniq > {output}"

rule calculate_r_using_plink:
    ''' Uses plink to calculate LD for an input list of variant IDs
    '''
    input:
        bfiles = rules.get_1000G_from_GCS.output,
        varfile = rules.write_variant_list.output
    output:
        expand(tmpdir + '/' + str(config['version']) + '/ld/ld_each_variant/{varid}.ld.tsv.gz',
               varid=varid_list)
    params:
        bfile_pref=tmpdir + '/{version}/ld/1000Genomep3/POPULATION/POPULATION.CHROM.1000Gp3.20130502'.format(version=config['version']),
        pops=hap1000G_pops,
        ld_window=config['ld_window'],
        min_r2=config['min_r2'],
        outdir = tmpdir + '/{version}/ld/ld_each_variant'.format(version=config['version'])
    threads: 300 # This is the max threads and will be scaled down by --cores argument
    shell:
        'python scripts/calc_ld_1000G.v2.py '
        '--varfile {input.varfile} '
        '--bfile {params.bfile_pref} '
        '--pops {params.pops} '
        '--ld_window {params.ld_window} '
        '--min_r2 {params.min_r2} '
        '--max_cores {threads} '
        '--outdir {params.outdir} '
        '--delete_temp'

rule concat_ld_scores:
    ''' Concat LD caluclated using plink to a single table
    '''
    input:
        rules.calculate_r_using_plink.output
    output:
        tmpdir + '/{version}/ld/top_loci_variants.ld.gz'
    params:
        in_pattern = tmpdir + '/' + str(config['version']) + '/ld/ld_each_variant/\*.ld.tsv.gz'
    shell:
        'python scripts/merge_ld_outputs.py '
        '--inpattern {params.in_pattern} '
        '--output {output}'

rule calc_study_specific_weighted_r2:
    ''' Merge LD to manifest and calculate overall weighted R2
    '''
    input:
        ld = rules.concat_ld_scores.output,
        manifest=in_manifest
    output:
        tmpdir + '/{version}/ld/study_weighted.ld.gz'
    shell:
        'python scripts/calc_study_weighted_ld.py '
        '--in_ld {input.ld} '
        '--in_manifest {input.manifest} '
        '--out {output}'

rule weight_studies_to_final:
    ''' Make finalised output of LD table
    '''
    input:
        ld = rules.calc_study_specific_weighted_r2.output,
        manifest = in_manifest
    output:
        'output/{version}/ld.parquet'
    params:
        min_r2=config['min_r2']
    shell:
        'python scripts/format_ld_table.py '
        '--inf {input.ld} '
        '--in_manifest {input.manifest} '
        '--outf {output} '
        '--min_r2 {params.min_r2}'

rule ld_to_GCS:
    ''' Copy to GCS
    '''
    input:
        rules.weight_studies_to_final.output
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/ld.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
