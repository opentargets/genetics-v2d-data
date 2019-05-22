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
        gwascat = rules.get_gwas_cat_assoc.output,
        invar = rules.get_variant_index.output
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
        gwascat = rules.get_gwas_cat_assoc.output,
        invar = rules.extract_gwascat_rsids_from_variant_index.output
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
        gwas = rules.annotate_gwas_cat_with_variant_ids.output,
        id_lut = rule.make_gwas_cat_studies_table.output.lut
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
        rules.convert_gwas_catalog_to_standard.output.out_assoc
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

rule make_summarystat_toploci_table:
    ''' Converts the toploci table produce from the finemapping pipeline to
        standardised format. Study table is need to know if a study is
        case-control or not.
    '''
    input:
        toploci=GSRemoteProvider().remote(config['toploci'], keep_local=KEEP_LOCAL),
        study_info = rules.merge_study_tables.output
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
        gwascat = rules.cluster_gwas_catalog.output.out_assoc,
        sumstat = rules.make_summarystat_toploci_table.output
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
        rules.merge_gwascat_and_sumstat_toploci.output
    output:
        GSRemoteProvider().remote(
            '{gs_dir}/{{version}}/toploci.parquet'.format(gs_dir=config['gs_dir'])
            )
    shell:
        'cp {input} {output}'
