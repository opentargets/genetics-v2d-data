
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

rule process_ld:
    ''' Process the LD table, including:
        - study/population specific weighting
        - PICS LD fine mapping metho
        - output parquet
    '''
    input:
        ld_files=rules.calculate_r_using_plink.output,
        manifest=in_manifest,
        toploci='output/{version}/toploci.parquet'
    output:
        directory('output/{version}/ld.parquet')
    params:
        in_ld_pattern = tmpdir + '/{version}/ld/ld_each_variant/*.ld.tsv.gz'.format(version=config['version']),
        min_r2 = config['min_r2']
    shell:
        'python scripts/process_ld.py '
        '--in_ld_pattern {params.in_ld_pattern} '
        '--in_manifest {input.manifest}'
        '--in_top_loci {input.toploci}'
        '--min_r2 {params.min_r2}'
        '--out {output}'

