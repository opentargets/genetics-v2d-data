import pandas as pd

rule get_gwas_cat_assoc:
    ''' Download GWAS catalog association file
    '''
    input:
        FTPRemoteProvider().remote('ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv')
    output:
        tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv'
    shell:
        'cp {input} {output}'

rule get_ensembl_variation_grch37:
    ''' Download all Ensembl variation data
    '''
    input:
        FTPRemoteProvider().remote('ftp://ftp.ensembl.org/pub/grch37/update/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz')
    output:
        tmpdir + '/Ensembl.Homo_sapiens.grch37.vcf.gz'
    shell:
        'cp {input} {output}'

rule get_HRC_variation_grch37:
    ''' Download all HRC variation data
    '''
    input:
        FTPRemoteProvider().remote('ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz')
    output:
        tmpdir + '/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz'
    shell:
        'cp {input} {output}'

rule extract_gwas_rsids_from_HRC:
    ''' Makes set of GWAS Catalog rsids and chrom:pos strings. Then reads
        these from the VCF file. Takes ~10 mins.
    '''
    input:
        gwascat= tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv',
        vcf= tmpdir + '/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz'
    output:
        tmpdir + '/{version}/HRC.r1-1.GRCh37.gwasCat_only.vcf.gz'
    shell:
        'pypy3 scripts/extract_from_vcf.py '
        '--gwas {input.gwascat} '
        '--vcf {input.vcf} '
        '--out {output}'

rule extract_gwas_rsids_from_Ensembl:
    ''' Makes set of GWAS Catalog rsids and chrom:pos strings. Then reads
        these from the VCF file. Takes ~10 mins.
    '''
    input:
        gwascat= tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv',
        vcf= tmpdir + '/Ensembl.Homo_sapiens.grch37.vcf.gz'
    output:
        tmpdir + '/{version}/Ensembl.Homo_sapiens.grch37.gwasCat_only.vcf.gz'
    shell:
        'pypy3 scripts/extract_from_vcf.py '
        '--gwas {input.gwascat} '
        '--vcf {input.vcf} '
        '--out {output}'

rule annotate_gwas_cat_with_variant_ids:
    ''' Annotates rows in the gwas catalog assoc file with variant IDs from
        a VCF
    '''
    input:
        gwascat= tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv',
        vcf_hrc= tmpdir + '/{version}/HRC.r1-1.GRCh37.gwasCat_only.vcf.gz',
        vcf_ensembl= tmpdir + '/{version}/Ensembl.Homo_sapiens.grch37.gwasCat_only.vcf.gz'
    output:
        tmpdir + '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv'
    shell:
        'python scripts/annotate_gwascat_varaintids.py '
        '--gwas {input.gwascat} '
        '--vcf_hrc {input.vcf_hrc} '
        '--vcf_ensembl {input.vcf_ensembl} '
        '--out {output}'

rule convert_gwas_catalog_to_standard:
    ''' Outputs the GWAS Catalog association data into a standard format
    '''
    input:
        tmpdir + '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv'
    output:
        out_assoc = tmpdir + '/{version}/gwas-catalog-associations_ot-format.tsv',
        log = 'logs/{version}/gwas-cat-assocs.log'
    shell:
        'python scripts/format_gwas_assoc.py '
        '--inf {input} '
        '--outf {output.out_assoc} '
        '--log {output.log}'

rule convert_nealeUKB_to_standard:
    ''' Converts the credible set results into a table of top loci in standard
        format.
    '''
    input:
        GSRemoteProvider().remote(config['credset'], keep_local=KEEP_LOCAL) # DEBUG
    output:
        tmpdir + '/{version}/nealeUKB-associations_ot-format.tsv'
    shell:
        'python scripts/format_nealeUKB_assoc.py '
        '--inf {input} '
        '--outf {output}'

rule merge_gwascat_and_nealeUKB_toploci:
    ''' Merges association files from gwas_cat and neale UKB
    '''
    input:
        gwas = tmpdir + '/{version}/gwas-catalog-associations_ot-format.tsv',
        neale = tmpdir + '/{version}/nealeUKB-associations_ot-format.tsv'
    output:
        'output/{version}/toploci.tsv'
    run:
        # Load
        gwas = pd.read_csv(input['gwas'], sep='\t', header=0)
        neale = pd.read_csv(input['neale'], sep='\t', header=0)
        # Merge
        merged = pd.concat([gwas, neale], sort=False)
        # Save
        merged.to_csv(output[0], sep='\t', index=None)

rule toploci_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/toploci.tsv'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/toploci.tsv'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
