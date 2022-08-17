### Disclarimer: this is a work in-progress version of the ingest script.
# Scope of this script:
# - Read and process GWAS Catalog associations and studies.
# - Save TopLoci table.
# - Save study table.

## Summary of the logic:
# - Read GWAS Catalog associations
# - Read GWAS Catalog study table.
# - Excluding a set of studies from the study table.
# - Be flexible to handle unpublished studies if needed.

import logging
import sys

import numpy as np
from functools import reduce

import pyspark.sql
from pyspark.sql import DataFrame
import pyspark.sql.types as T
import pyspark.sql.functions as F
from pyspark import SparkFiles

# These parameters will come as command line arguments:
VARIANT_ANNOTATION = 'gs://ot-team/dsuveges/variant_annotation/2022-08-11'

# associations = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'  # URL to GWAS Catalog associations
GWAS_CATALOG_ASSOCIATIONS = 'gs://ot-team/dsuveges/v2d_files/gwas_associations_2022-06-27.tsv'

# Pointers to study files:
STUDY_FILE = 'gs://ot-team/dsuveges/v2d_files/gwas_studies_v1.03_2022_08_16.tsv'
UNPUBLISHED_STUDY_FILE = 'gs://ot-team/dsuveges/v2d_files/gwas_studies_v1.03_unpublished_2022_08_16.tsv'

##
## Temporary output files:
##
OUTPUT_PATH = 'gs://ot-team/dsuveges/v2d_files/'

# These parameters will be read from a config file:


# List of studies to be dropped:
DROPPED_STUDIES = [
    'GCST005806',  # Sun et al pQTL study
    'GCST005837'  # Huang et al IBD study
]

# This flag will allow to run the code for PPP allowing the ingestion of unpublished studies:
INGEST_UNPUBLISHED_STUDIES = False

# P-value threshold for associations:
PVALCUTOFF = 5e-8

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

STUDY_COLUMNS_MAP = {
    # 'DATE ADDED TO CATALOG': 'date_added_to_catalog',
    'PUBMED ID': 'pubmed_id',
    'FIRST AUTHOR': 'first_author',
    'DATE': 'publication_date',
    'JOURNAL': 'journal',
    'LINK': 'link',
    'STUDY': 'study',
    'DISEASE/TRAIT': 'disease_trait',
    'INITIAL SAMPLE SIZE': 'initial_sample_size',
    # 'REPLICATION SAMPLE SIZE': 'replication_sample_size',
    # 'PLATFORM [SNPS PASSING QC]': 'platform',
    # 'ASSOCIATION COUNT': 'association_count',
    'MAPPED_TRAIT': 'mapped_trait',
    'MAPPED_TRAIT_URI': 'mapped_trait_uri',
    'STUDY ACCESSION': 'study_accession',
    # 'GENOTYPING TECHNOLOGY': 'genotyping_technology',
    'SUMMARY STATS LOCATION': 'summary_stats_location',
    # 'SUBMISSION DATE': 'submission_date',
    # 'STATISTICAL MODEL': 'statistical_model',
    'BACKGROUND TRAIT': 'background_trait',
    'MAPPED BACKGROUND TRAIT': 'mapped_background_trait',
    'MAPPED BACKGROUND TRAIT URI': 'mapped_background_trait_uri',
}

def read_study_table(read_unpublished: bool = False) -> pyspark.sql.DataFrame:
    """
    This function reads the study table and returns a Spark DataFrame with the columns renamed to match the names
    in the `STUDY_COLUMNS_MAP` dictionary

    Args:
      read_unpublished (bool): bool=False. Defaults to False. Indicating if unpublished studies should be read.

    Returns:
      A dataframe with the columns specified in STUDY_COLUMNS_MAP.
    """

    # Reading the study table:
    study_table = spark.read.csv(STUDY_FILE, sep='\t', header=True)

    if read_unpublished:

        # Reading the unpublished study table:
        unpublished_study_table = (
            spark.read.csv(UNPUBLISHED_STUDY_FILE, sep='\t', header=True)

            # The column names are expected to be the same in both tables, so this should be gone at some point:
            .withColumnRenamed('MAPPED TRAIT', 'MAPPED_TRAIT')
            .withColumnRenamed('MAPPED TRAIT URI', 'MAPPED_TRAIT_URI')
        )

        # Expression to replace "not yet cureated" string with nulls:
        expressions = map(
            lambda column: (column, F.when(F.col(column) == 'not yet curated', F.lit(None)).otherwise(F.col(column))),
            STUDY_COLUMNS_MAP.keys()
        )
        unpublished_study_table = reduce(lambda DF, value: DF.withColumn(*value), expressions, unpublished_study_table)

        study_table = study_table.union(unpublished_study_table)

    # Selecting and renaming relevant columns:
    return (
        study_table
        .select(*[F.col(old_name).alias(new_name) for old_name, new_name in STUDY_COLUMNS_MAP.items()])
        .persist()
    )


def extract_sample_sizes(df: pyspark.sql.DataFrame) -> pyspark.sql.DataFrame:
    """
    It takes a dataframe with a column called `initial_sample_size` and returns a dataframe with columns
    `n_cases`, `n_controls`, and `n_samples`

    Args:
      df (pyspark.sql.DataFrame): pyspark.sql.DataFrame

    Returns:
      A dataframe with the columns:
        - n_cases
        - n_controls
        - n_samples
    """

    columns = df.columns

    return (
        df
        # .withColumn('samples', F.explode_outer(F.col('initial_sample_size')))
        .withColumn('samples', F.explode(F.split(F.col('initial_sample_size'), r',\s+')))

        # Extracting the sample size from the string:
        .withColumn('sample_size', F.regexp_extract(F.regexp_replace(F.col('samples'), ',', ''), r'[0-9,]+', 0).cast(T.IntegerType()))

        # Extracting number of cases:
        .withColumn('n_cases', F.when(F.col('samples').contains('cases'), F.col('sample_size')).otherwise(F.lit(0)))
        .withColumn('n_controls', F.when(F.col('samples').contains('controls'), F.col('sample_size')).otherwise(F.lit(0)))

        .groupBy(columns)
        .agg(
            F.sum('n_cases').alias('n_cases'),
            F.sum('n_controls').alias('n_controls'),
            F.sum('sample_size').alias('n_samples')
        )
        .persist()
    )

def read_GWAS_associations() -> DataFrame:
    """
    It reads the GWAS Catalog association dataset, selects and renames columns, casts columns, and
    applies some pre-defined filters on the data

    Returns:
      A dataframe with the GWAS Catalog associations.
    """

    # Reading and filtering associations:
    association_df = (
        spark.read.csv(GWAS_CATALOG_ASSOCIATIONS, sep='\t', header=True)

        # Select and rename columns:
        .select(*[F.col(old_name).alias(new_name) for old_name, new_name in ASSOCIATION_COLUMNS_MAP.items()])

        # Cast columns:
        .withColumn('pvalue_mlog', F.col('pvalue_mlog').cast(T.FloatType()))

        # Apply some pre-defined filters on the data:
        .filter(
            (~F.col('chr_id').contains(' x '))  # Dropping associations based on variant x variant interactions
            & (F.col('pvalue_mlog') > -np.log10(PVALCUTOFF))  # Dropping sub-significant associations
            & (F.col('chr_pos').isNotNull() & F.col('chr_id').isNotNull())  # Dropping associations without genomic location
        )
    )

    # Providing stats on the filtered association dataset:
    logging.info(f'Number of associations: {association_df.count()}')
    logging.info(f'Number of studies: {association_df.select("study_accession").distinct().count()}')
    logging.info(f'Number of variants: {association_df.select("snp_ids").distinct().count()}')
    logging.info(f'Assocation: {association_df.show(2, False, True)}')

    logging.info('Processing associaitons:')
    return association_df

def parse_associations(association_df: DataFrame) -> DataFrame:
    # Processing associations:
    parsed_associations = (

        # spark.read.csv(associations, sep='\t', header=True)
        association_df

        # Adding association identifier for future deduplication:
        .withColumn('association_id', f.monotonically_increasing_id())

        # Processing variant related columns:
        #   - Sorting out current rsID field: <- why do we need this? rs identifiers should always come from the GnomAD dataset.
        #   - Removing variants with no genomic mappings -> losing ~3% of all associations
        #   - Multiple variants can correspond to a single association.
        #   - Variant identifiers are stored in the SNPS column, while the mapped coordinates are stored in the CHR_ID and CHR_POS columns.
        #   - All these fields are split into arrays, then they are paired with the same index eg. first ID is paired with first coordinate, and so on
        #   - Then the association is exploded to all variants.
        #   - The risk allele is extracted from the 'STRONGEST SNP-RISK ALLELE' column.

        # The current snp id field is just a number at the moment (stored as a string). Adding 'rs' prefix if looks good.
        .withColumn(
            'snp_id_current',
            F.when(F.col('snp_id_current').rlike('^[0-9]*$'), F.format_string('rs%s', F.col('snp_id_current')))
            .otherwise(F.col('snp_id_current'))
        )

        # Variant notation (chr, pos, snp id) are split into array:
        .withColumn('chr_id', F.split(F.col('chr_id'), ';'))
        .withColumn('chr_pos', F.split(F.col('chr_pos'), ';'))
        .withColumn('strongest_snp_risk_allele', F.split(F.col('strongest_snp_risk_allele'), '; '))
        .withColumn('snp_ids', F.split(F.col('snp_ids'), '; '))

        # Variant fields are joined together in a matching list, then extracted into a separate rows again:
        .withColumn('VARIANT', F.explode(F.arrays_zip('chr_id', 'chr_pos', 'strongest_snp_risk_allele', 'snp_ids')))

        # Updating variant columns:
        .withColumn('snp_ids', F.col('VARIANT.snp_ids'))
        .withColumn('chr_id', F.col('VARIANT.chr_id'))
        .withColumn('chr_pos', F.col('VARIANT.chr_pos').cast(T.IntegerType()))
        .withColumn('strongest_snp_risk_allele', F.col('VARIANT.strongest_snp_risk_allele'))

        # Extracting risk allele:
        .withColumn('risk_allele', F.split(F.col('strongest_snp_risk_allele'), '-').getItem(1))

        # Create a unique set of SNPs linked to the assocition:
        .withColumn(
            'rsid_gwas_catalog',
            F.array_distinct(F.array(
                F.split(F.col('strongest_snp_risk_allele'), '-').getItem(0),
                F.col('snp_id_current'),
                F.col('snp_ids')
            ))
        )

        # Processing EFO terms:
        #   - Multiple EFO terms can correspond to a single association.
        #   - EFO terms are stored as full URIS, separated by semicolons.
        #   - Associations are exploded to all EFO terms.
        #   - EFO terms in the study table is not considered as association level EFO annotation has priority (via p-value text)

        # Process EFO URIs:
        .withColumn('efo', F.explode(F.expr(r"regexp_extract_all(mapped_trait_uri, '([A-Z]+_[0-9]+)')")))

        # Splitting p-value into exponent and mantissa:
        .withColumn('exponent', F.split(F.col('p_value'), 'E').getItem(1).cast('integer'))
        .withColumn('mantissa', F.split(F.col('p_value'), 'E').getItem(0).cast('float'))

        # Cleaning up:
        .drop('mapped_trait_uri', 'strongest_snp_risk_allele', 'VARIANT')
        .persist()
    )

    # Providing stats on the filtered association dataset:
    logging.info(f'Number of associations: {parsed_associations.count()}')
    logging.info(f'Number of studies: {parsed_associations.select("study_accession").distinct().count()}')
    logging.info(f'Number of variants: {parsed_associations.select("snp_ids").distinct().count()}')
    logging.info(f'Assocation: {parsed_associations.show(2, False, True)}')

    return parsed_associations

def map_associations(parsed_associations: DataFrame) -> DataFrame:
    # Loading variant annotation and join with parsed associations:
    logging.info('Loading variant annotation:')

    # Reading and joining variant annotation:
    variants = (
        spark.read.parquet(VARIANT_ANNOTATION)
        .select(
            F.col('chrom_b38').alias('chr_id'),
            F.col('pos_b38').alias('chr_pos'),
            F.col('rsid').alias('rsid_gnomad'),
            F.col('ref').alias('ref'),
            F.col('alt').alias('alt'),
            F.col('id').alias('variant_id')
        )
    )

    mapped_associations = (
        parsed_associations
        .join(variants, on=['chr_id', 'chr_pos'], how='left')
        .persist()
    )

    return mapped_associations

def main():

    # Read GWAS associations and process:
    gwas_associations = read_GWAS_associations()
    parsed_gwas_associations = parse_associations(gwas_associations)

    # Debug: save data parsed associations:
    # parsed_gwas_associations.write.mode('overwrite').parquet(f'{OUTPUT_PATH}/parsed_associations.parquet')

    # Mapping GWAS Catalog variants to GnomAD variants:
    mapped_associations = map_associations(parsed_gwas_associations)

    # Debug: save mapped data:
    # mapped_associations.write.mode('overwrite').parquet(f'{OUTPUT_PATH}/mapped_parsed_associations.parquet')

    # Debug: read mapped associations:
    # mapped_associations = spark.read.parquet(f'{OUTPUT_PATH}/mapped_parsed_associations.parquet').persist()

    # Drop association/variant mappings, where the alleles are discordant:
    mapped_associations = (
        mapped_associations

        # Dropping associations with no mapped variants:
        .filter(F.col('alt').isNotNull())

        # Adding column with the reverse-complement of the risk allele:
        .withColumn(
            'risk_allele_reverse_complement',
            F.when(
                F.col('risk_allele').rlike(r'^[ACTG]+$'),
                F.reverse(F.translate(F.col('risk_allele'), "ACTG", "TGAC"))
            )
            .otherwise(F.col('risk_allele'))
        )

        # Adding columns flagging concordance:
        .withColumn(
            'is_concordant',
            # If risk allele is found on the positive strand:
            F.when((F.col('risk_allele') == F.col('ref')) | (F.col('risk_allele') == F.col('alt')), True)
            # If risk allele is found on the negative strand:
            .when(
                (F.col('risk_allele_reverse_complement') == F.col('ref'))
                | (F.col('risk_allele_reverse_complement') == F.col('alt')),
                True
            )
            # If risk allele is ambiguous, still accepted:
            .when(F.col('risk_allele') == '?', True)
            # Allele is discordant:
            .otherwise(False)
        )

        # Dropping discordant associations:
        .filter((F.col('is_concordant') == True) & (F.col('risk_allele') == '?'))
        .drop('is_concordant', 'risk_allele_reverse_complement')

        # Adding column for variant id:
        .withColumn('variant_id', F.concat_ws('_', F.col('chr_id'), F.col('chr_pos'), F.col('ref'), F.col('alt')))

        .show(1, False, True)
    )

    # Debug: saving concordant mappings:
    mapped_associations.write.mode('overwrite').parquet(f'{OUTPUT_PATH}/concordant_associations.parquet')

    # Reading study table:
    study_table = read_study_table(INGEST_UNPUBLISHED_STUDIES)

    # Extract sample sizes:
    study_table = extract_sample_sizes(study_table)


if __name__ == '__main__':
    # Init logging:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr
    )

    # Init spark
    global spark
    spark = (
        pyspark.sql.SparkSession
        .builder
        .master("local[*]")
        .getOrCreate()
    )

    main()
