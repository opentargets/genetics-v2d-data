import os

rule download_credible_set_directory:
    ''' Input is a directory of json parts. Need to do some trickery to download from GCS.
    '''
    input:
        GSRemoteProvider().remote(config['credsets'], keep_local=KEEP_LOCAL)
    output:
        temp(directory('tmp/finemapping/{version}/credset')) # ~ 2GB
        # directory('tmp/finemapping/{version}/credset')
    run:
        # Create input name 
        in_path = os.path.dirname(input[0])
        if not in_path.startswith('gs://'):
            in_path = 'gs://' + in_path
        # Make output directory
        os.makedirs(output[0], exist_ok=True)
        # Download using 
        shell(
            'gsutil -m rsync -r {src} {dest}'.format(
                src=in_path,
                dest=output[0]
            )
        )

rule convert_finemapping_to_standard:
    ''' Extract required fields from credible set file
    '''
    input:
        rules.download_credible_set_directory.output
    output:
        directory('output/{version}/finemapping.parquet')
    shell:
        'python scripts/format_finemapping_table.py '
        '--inf {input} '
        '--outf {output}'

rule finemap_to_GCS:
    ''' Copy to GCS
    '''
    input:
        rules.convert_finemapping_to_standard.output
    output:
        GSRemoteProvider().remote(
            directory(
                '{gs_dir}/{{version}}/finemapping.parquet'.format(gs_dir=config['gs_dir'])
            )
        )
    params:
        outgs = config['gs_dir'] + '/{version}/finemapping.parquet'
    shell:
        'gsutil -m rsync -r {input} {params}'
