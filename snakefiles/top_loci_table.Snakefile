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

rule get_variant_index:
    ''' Download variant index site list
    '''
    input:
        GSRemoteProvider().remote(config['var_index_sitelist'], keep_local=KEEP_LOCAL)
    output:
        tmpdir + '/variant-annotation.sitelist.tsv.gz'
    shell:
        'cp {input} {output}'

rule extract_gwascat_rsids_from_variant_index:
    ''' Makes set of GWAS Catalog rsids and chrom:pos strings. Then reads
        these from the variant index. Takes ~2 mins.
    '''
    input:
        gwascat= tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv',
        invar= tmpdir + '/variant-annotation.sitelist.tsv.gz'
    output:
        tmpdir + '/{version}/variant-annotation.sitelist.GWAScat.tsv.gz'
    shell:
        'pypy3 scripts/extract_from_variant-index.py '
        '--gwas {input.gwascat} '
        '--vcf {input.invar} '
        '--out {output}'

rule annotate_gwas_cat_with_variant_ids:
    ''' Annotates rows in the gwas catalog assoc file with variant IDs from
        a VCF
    '''
    input:
        gwascat= tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv',
        invar= tmpdir + '/{version}/variant-annotation.sitelist.GWAScat.tsv.gz'
    output:
        tmpdir + '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv'
    shell:
        'python scripts/annotate_gwascat_varaintids.py '
        '--gwas {input.gwascat} '
        '--invar {input.invar} '
        '--out {output}'

rule convert_gwas_catalog_to_standard:
    ''' Outputs the GWAS Catalog association data into a standard format
    '''
    input:
        gwas=tmpdir + '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv',
        id_lut=tmpdir + '/{version}/gwas-catalog_study_id_lut.tsv'
    output:
        out_assoc = tmpdir + '/{version}/gwas-catalog-associations_ot-format.tsv',
        log = 'logs/{version}/gwas-cat-assocs.log'
    shell:
        'python scripts/format_gwas_assoc.py '
        '--inf {input.gwas} '
        '--id_lut {input.id_lut} '
        '--outf {output.out_assoc} '
        '--log {output.log}'

rule cluster_gwas_catalog:
    ''' GWAS Catalog contains none independent loci from curation error result-
        ing in summary statistics being bulk imported from supplementary tables.
        To correct this, an additional clustering step has been added.

        In addition to clustering, this script also filters by p-value.
    '''
    input:
        tmpdir + '/{version}/gwas-catalog-associations_ot-format.tsv',
    output:
        out_assoc=tmpdir + '/{version}/gwas-catalog-associations_ot-format.clustered.tsv',
        log='logs/{version}/gwas-cat-assocs_clustering.log'
    params:
        min_p=config['gwas_cat_min_pvalue'],
        dist=config['gwas_cat_cluster_dist_kb'],
        min_loci=config['gwas_cat_cluster_min_loci'],
        multi_prop=config['gwas_cat_cluster_multi_proportion']
    shell:
        'python scripts/cluster_gwas_catalog_associations.py '
        '--inf {input} '
        '--outf {output.out_assoc} '
        '--min_p {params.min_p} '
        '--cluster_dist_kb {params.dist} '
        '--cluster_min_loci {params.min_loci} '
        '--cluster_multi_prop {params.multi_prop} '
        '--log {output.log}'

# rule convert_nealeUKB_to_standard:
#     ''' Converts the credible set results into a table of top loci in standard
#         format.
#     '''
#     input:
#         credset=GSRemoteProvider().remote(config['credset'], keep_local=KEEP_LOCAL),
#         study_info=tmpdir + '/{version}/nealeUKB_study_table.tsv'
#     output:
#         tmpdir + '/{version}/nealeUKB-associations_ot-format.tsv'
#     shell:
#         'python scripts/format_nealeUKB_assoc.py '
#         '--inf {input.credset} '
#         '--study_info {input.study_info} '
#         '--outf {output}'

rule make_summarystat_toploci_table:
    ''' Converts the toploci table produce from the finemapping pipeline to
        standardised format. Study table is need to know if a study is
        case-control or not.
    '''
    input:
        toploci=GSRemoteProvider().remote(config['toploci'], keep_local=KEEP_LOCAL),
        study_info = tmpdir + '/{version}/merged_study_table.tsv'
    output:
        tmpdir + '/{version}/sumstat-associations_ot-format.tsv'
    shell:
        'python scripts/format_sumstat_toploci_assoc.py '
        '--inf {input.toploci} '
        '--study_info {input.study_info} '
        '--outf {output}'

rule merge_gwascat_and_sumstat_toploci:
    ''' Merges associations from gwas_cat and sumstat derived
    '''
    input:
        gwascat = tmpdir + '/{version}/gwas-catalog-associations_ot-format.clustered.tsv',
        sumstat = tmpdir + '/{version}/sumstat-associations_ot-format.tsv'
    output:
        'output/{version}/toploci.parquet'
    shell:
        'python scripts/merge_top_loci_tables.py '
        '--in_gwascat {input.gwascat} '
        '--in_sumstat {input.sumstat} '
        '--output {output}'

rule toploci_to_GCS:
    ''' Copy to GCS
    '''
    input:
        'output/{version}/toploci.parquet'
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/toploci.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
