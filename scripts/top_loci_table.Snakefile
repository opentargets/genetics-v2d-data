import pandas as pd

rule get_gwas_cat_assoc:
    ''' Download GWAS catalog association file
    '''
    input:
        FTPRemoteProvider().remote('ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv')
    output:
        tmpdir + '/gwas-catalog-associations_ontology-annotated.{version}.tsv'
    shell:
        'cp {input} {output}'

# rule get_ensembl_variation_grch37:
#     ''' Download all Ensembl variation data
#     '''
#     input:
#         FTPRemoteProvider().remote('ftp://ftp.ensembl.org/pub/grch37/update/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz')
#     output:
#         tmpdir + '/Homo_sapiens.grch37.vcf.gz'
#     shell:
#         'cp {input} {output}'

rule get_HRC_variation_grch37:
    ''' Download all HRC variation data
    '''
    input:
        FTPRemoteProvider().remote('ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz')
    output:
        tmpdir + '/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz'
    shell:
        'cp {input} {output}'

rule extract_gwas_rsids:
    ''' Makes set of GWAS Catalog rsids and chrom:pos strings. Then reads
        these from the VCF file. Takes ~10 mins.
    '''
    input:
        gwascat= tmpdir + '/gwas-catalog-associations_ontology-annotated.{version}.tsv',
        vcf= tmpdir + '/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz'
    output:
        tmpdir + '/HRC.r1-1.GRCh37.gwasCat_only.{version}.vcf.gz'
    shell:
        'python scripts/extract_from_vcf.py '
        '--gwas {input.gwascat} '
        '--vcf {input.vcf} '
        '--out {output}'

rule annotate_gwas_cat_with_variant_ids:
    ''' Annotates rows in the gwas catalog assoc file with variant IDs from
        a VCF
    '''
    input:
        gwascat= tmpdir + '/gwas-catalog-associations_ontology-annotated.{version}.tsv',
        vcf37= tmpdir + '/HRC.r1-1.GRCh37.gwasCat_only.{version}.vcf.gz'
    output:
        tmpdir + '/gwas-catalog-associations_ontology_variantID-annotated.{version}.tsv'
    shell:
        'python scripts/annotate_gwascat_varaintids.py '
        '--gwas {input.gwascat} '
        '--vcf_grch37 {input.vcf37} '
        '--out {output}'

rule convert_gwas_catalog_to_standard:
    ''' Outputs the GWAS Catalog association data into a standard format
    '''
    input:
        tmpdir + '/gwas-catalog-associations_ontology_variantID-annotated.{version}.tsv'
    output:
        out_assoc = tmpdir + '/gwas-catalog-associations_ot-format.{version}.tsv',
        log = 'logs/gwas-cat-assocs.{version}.log'
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
        GSRemoteProvider().remote(config['credset'], keep_local=True) # DEBUG
    output:
        tmpdir + '/nealeUKB-associations_ot-format.{version}.tsv'
    shell:
        'python scripts/format_nealeUKB_assoc.py '
        '--inf {input} '
        '--outf {output}'

rule merge_gwascat_and_nealeUKB_toploci:
    ''' Merges association files from gwas_cat and neale UKB
    '''
    input:
        gwas = tmpdir + '/gwas-catalog-associations_ot-format.{version}.tsv',
        neale = tmpdir + '/nealeUKB-associations_ot-format.{version}.tsv'
    output:
        'output/ot_genetics_toploci_table.{version}.tsv'
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
        'output/ot_genetics_toploci_table.{version}.tsv'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/ot_genetics_toploci_table.{{version}}.tsv'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
