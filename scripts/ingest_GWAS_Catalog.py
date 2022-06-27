### Disclarimer: this is a work in-progress version of the ingest script.
import logging
import sys

import numpy as np

import pyspark.sql
import pyspark.sql.types as t
import pyspark.sql.functions as f
from pyspark import SparkFiles


# Init logging:
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr
)

spark = (
    pyspark.sql.SparkSession
    .builder
    .master("local[*]")
    .getOrCreate()
)

# These parameters will come as command line arguments:
variant_annotation = 'gs://genetics-portal-dev-analysis/dsuveges/variant_annotation/2022-06-22'
# associations = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'  # URL to GWAS Catalog associations
associations = 'gs://ot-team/dsuveges/v2d_files/gwas_associations_2022-06-27.tsv'

##
## Temporary output files:
##
OUTPUT_PATH = 'gs://ot-team/dsuveges/v2d_files/'

# These parameters will be read from a config file:


# Ingesting GWAS Catalog associations:
DROPPED_STUDIES = [
    'GCST005806',  # Sun et al pQTL study
    'GCST005837'  # Huang et al IBD study
]

# This flag will allow to run the code for PPP allowing the ingestion of unpublished studies:
INGEST_UNPUBLISHED_STUDIES = False

# P-value threshold for associations:
PVALCUTOFF = 8e-8


ASSOCIATION_COLUMNS_MAP = {
    # Variant related columns:
    'STRONGEST SNP-RISK ALLELE': 'strongest_snp_risk_allele',  # variant id and the allele is extracted (; separated list)
    'CHR_ID': 'chr_id',  # Mapped genomic location of the variant (; separated list)
    'CHR_POS': 'chr_pos',
    'RISK ALLELE FREQUENCY': 'risk_allele_frequency',
    'CNV': 'cnv',  # Flag if a variant is a copy number variant
    'SNP_ID_CURRENT': 'snp_id_current',  #
    'SNPS': 'snp_ids',  # List of all SNPs associated with the variant

    # Study related columns:
    'STUDY ACCESSION': 'study_accession',

    # Disease/Trait related columns:
    'DISEASE/TRAIT': 'disease_trait',  # Reported trait of the study
    'MAPPED_TRAIT_URI': 'mapped_trait_uri',  # Mapped trait URIs of the study
    'MERGED': 'merged',
    'P-VALUE (TEXT)': 'p_value_text',  # Extra information on the association.

    # Association details:
    'P-VALUE': 'p_value',  # p-value of the association, string: split into exponent and mantissa.
    'PVALUE_MLOG': 'pvalue_mlog',  # -log10(p-value) of the association, float
    'OR or BETA': 'or_beta',  # Effect size of the association, Odds ratio or beta
    '95% CI (TEXT)': 'confidence_interval',  # Confidence interval of the association, string: split into lower and upper bound.

    'CONTEXT': 'context',
}

# Reading and filtering associations:
association_df = (
    spark.read.csv(associations, sep='\t', header=True)

    # Select and rename columns:
    .select(*[f.col(old_name).alias(new_name) for old_name, new_name in ASSOCIATION_COLUMNS_MAP.items()])

    # Cast columns:
    .withColumn('pvalue_mlog', f.col('pvalue_mlog').cast(t.FloatType()))

    # Apply some pre-defined filters on the data:
    .filter(
        (~f.col('chr_id').contains(' x '))  # Dropping associations based on variant x variant interactions
        & (f.col('pvalue_mlog') > -np.log10(PVALCUTOFF))  # Dropping sub-significant associations
        & (f.col('chr_pos').isNotNull() & f.col('chr_id').isNotNull())  # Dropping associations without genomic location
    )
)

# Providing stats on the filtered association dataset:
logging.info(f'Number of associations: {association_df.count()}')
logging.info(f'Number of studies: {association_df.select("study_accession").distinct().count()}')
logging.info(f'Number of variants: {association_df.select("snp_ids").distinct().count()}')
logging.info(f'Assocation: {association_df.show(2, False, True)}')

logging.info('Processing associaitons:')

# Processing associations:
parsed_associations = (

    # spark.read.csv(associations, sep='\t', header=True)
    association_df

    # Adding association identifier for future deduplication:
    .withColumn('association_id', f.monotonically_increasing_id())

    # Processing variant related columns:
    #   - Sorting out current rsID field:
    #   - Removing variants with no genomic mappings -> losing ~3% of all associations
    #   - Multiple variants can correspond to a single association.
    #   - Variant identifiers are stored in the SNPS column, while the mapped coordinates are stored in the CHR_ID and CHR_POS columns.
    #   - All these fields are split into arrays, then they are paired with the same index eg. first ID is paired with first coordinate, and so on
    #   - Then the association is exploded to all variants.
    #   - The risk allele is extracted from the 'STRONGEST SNP-RISK ALLELE' column.

    # The current snp id field is just a number at the moment (stored as a string). Adding 'rs' prefix if looks good.
    .withColumn(
        'snp_id_current',
        f.when(f.col('snp_id_current').rlike('^[0-9]*$'), f.format_string('rs%s', f.col('snp_id_current')))
        .otherwise(f.col('snp_id_current'))
    )

    # Variant notation (chr, pos, snp id) are split into array:
    .withColumn('chr_id', f.split(f.col('chr_id'), ';'))
    .withColumn('chr_pos', f.split(f.col('chr_pos'), ';'))
    .withColumn('strongest_snp_risk_allele', f.split(f.col('strongest_snp_risk_allele'), '; '))
    .withColumn('snp_ids', f.split(f.col('snp_ids'), '; '))

    # Variant fields are joined together in a matching list, then extracted into a separate rows again:
    .withColumn('VARIANT', f.explode(f.arrays_zip('chr_id', 'chr_pos', 'strongest_snp_risk_allele', 'snp_ids')))

    # Updating variant columns:
    .withColumn('snp_ids', f.col('VARIANT.snp_ids'))
    .withColumn('chr_id', f.col('VARIANT.chr_id'))
    .withColumn('chr_pos', f.col('VARIANT.chr_pos').cast(t.IntegerType()))
    .withColumn('strongest_snp_risk_allele', f.col('VARIANT.strongest_snp_risk_allele'))

    # Extracting risk allele:
    .withColumn('risk_allele', f.split(f.col('strongest_snp_risk_allele'), '-').getItem(1))

    # Create a unique set of SNPs linked to the assocition:
    .withColumn(
        'rsid_gwas_catalog',
        f.array_distinct(f.array(
           f.split(f.col('strongest_snp_risk_allele'), '-').getItem(0),
           f.col('snp_id_current'),
           f.col('snp_ids')
        ))
    )

    # Processing EFO terms:
    #   - Multiple EFO terms can correspond to a single association.
    #   - EFO terms are stored as full URIS, separated by semicolons.
    #   - Associations are exploded to all EFO terms.
    #   - EFO terms in the study table is not considered as association level EFO annotation has priority (via p-value text)

    # Process EFO URIs:
    .withColumn('efo', f.explode(f.expr(r"regexp_extract_all(mapped_trait_uri, '([A-Z]+_[0-9]+)')")))

    # Splitting p-value into exponent and mantissa:
    .withColumn('exponent', f.split(f.col('p_value'), 'E').getItem(1).cast('integer'))
    .withColumn('mantissa', f.split(f.col('p_value'), 'E').getItem(0).cast('float'))

    # Cleaning up:
    .drop('mapped_trait_uri', 'strongest_snp_risk_allele', 'VARIANT')
    .persist()
)

# Providing stats on the filtered association dataset:
logging.info(f'Number of associations: {parsed_associations.count()}')
logging.info(f'Number of studies: {parsed_associations.select("study_accession").distinct().count()}')
logging.info(f'Number of variants: {parsed_associations.select("snp_ids").distinct().count()}')
logging.info(f'Assocation: {parsed_associations.show(2, False, True)}')

# Saving data:
# parsed_associations.write.mode('overwrite').parquet(f'{OUTPUT_PATH}/parsed_associations.parquet')

# Loading variant annotation and join with parsed associations:
logging.info('Loading variant annotation:')

# Reading and joining variant annotation:
variants = (
    spark.read.parquet(variant_annotation)
    .select(
        f.col('chrom_b38').alias('chr_id'),
        f.col('pos_b38').alias('chr_pos'),
        f.col('rsid').alias('rsid_gnomad'),
        f.col('ref').alias('ref'),
        f.col('alt').alias('alt')
    )
)

mapped_associations = (
    parsed_associations
    .join(variants, on=['chr_id', 'chr_pos'], how='left_outer')
    .persist()
)

mapped_associations.write.mode('overwrite').parquet(f'{OUTPUT_PATH}/mapped_parsed_associations.parquet')

