
#
# Get 1000 Genomes haplotypes from GCS -----------------------------------------
#

hap1000G_pops = ['EUR', 'EAS', 'AMR', 'AFR', 'SAS']
hap1000G_chroms = list(range(1, 23)) + ['X']

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
        outdir=tmpdir + '/{version}/ld/1000Genomep3'.format(version=config['version'])
    shell:
        'gsutil -m rsync -r gs://genetics-portal-input/1000Genomes_phase3/plink_format {params.outdir}'

#
# Calculate LD using plink -----------------------------------------------------
#

def make_input_bfiles(wildcards):
    chrom = wildcards['varid'].split('_')[0]
    bfiles = []
    for pop in hap1000G_pops:
        for ext in ['bed', 'bim', 'fam']:
            bfiles.append(
                tmpdir + '/{version}/ld/1000Genomep3/{pop}/{pop}.{chrom}.1000Gp3.20130502.{ext}'.format(
                    version=config['version'],
                    pop=pop,
                    chrom=chrom,
                    ext=ext
                )
            )
    return bfiles

rule calculate_r_using_plink:
    ''' Uses plink to calculate LD
    '''
    input: make_input_bfiles
    output:
        tmpdir + '/{version}/ld/plink_r_calc/{varid}/{varid}.index_var.ld.gz'
    params:
        bfile_pref=lambda wildcards: tmpdir + '/{version}/ld/1000Genomep3/POPULATION/POPULATION.CHROM.1000Gp3.20130502'.format(version=wildcards['version']),
        pops=hap1000G_pops,
        ld_window=config['ld_window'],
        min_r2=config['min_r2']
    shell:
        'python scripts/calc_ld_1000G.py '
        '--varid {wildcards.varid} '
        '--bfile {params.bfile_pref} '
        '--pops {params.pops} '
        '--ld_window {params.ld_window} '
        '--min_r2 {params.min_r2} '
        '--outf {output} '

rule concat_ld_scores:
    ''' Concat LD caluclated using plink to a single table
    '''
    input:
        expand(tmpdir + '/' + str(config['version']) + '/ld/plink_r_calc/{varid}/{varid}.index_var.ld.gz',
               varid=varid_list)
    output:
        tmpdir + '/{version}/ld/top_loci_variants.ld.gz'
    params:
        in_pattern=tmpdir + '/' + str(config['version']) + '/ld/plink_r_calc/\*/\*.index_var.ld.gz'
    shell:
        'python scripts/merge_ld_outputs.py '
        '--inpattern {params.in_pattern} '
        '--output {output}'

rule calc_study_specific_weighted_r2:
    ''' Merge LD to manifest and calculate overall weighted R2
    '''
    input:
        ld=tmpdir + '/{version}/ld/top_loci_variants.ld.gz',
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
        tmpdir + '/{version}/ld/study_weighted.ld.gz'
    output:
        'output/{version}/ld.parquet'
    params:
        min_r2=config['min_r2']
    shell:
        'python scripts/format_ld_table.py '
        '--inf {input} '
        '--outf {output} '
        '--min_r2 {params.min_r2}'

rule ld_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/ld.parquet'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/ld.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'