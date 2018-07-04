Variant-disease tables
======================

Tables:
  1. Variant-disease associations
    - GWAS Catalog
    - Neale et al UK Biobank top loci
  2. Study information table
  3. Lead variant to tag variant table
  4. Summary statistic table

#### Notes
- Do we only include UKB studies with an EFO

#### Top loci table
  1. Download all associations from: ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest
  2. Map rsid and chr:pos to variant_ids (GRCh37)
    - This is complicated by:
      (i) some GWAScat loci have mulitple rsid assignments,
      (ii) some rows don't have rsids, use chrom:pos (GRCh37 assumed) instead
      (iii) rsids/chrom:pos can map to multiple variant ids
            (e.g. rs773154059)
    - Notes:
      - Two associations don't have the same number of SNPs in SNPS, CHR_ID and
        CHR_POS. These have been dropped as they are ambiguous.
      - I am assuming that variants with id e.g. chr21:43693789 are from build GRCh37.
      - Associations that do not have a rsid or chr:pos will be dropped
      - Associations that do not have a match in the Ensembl VCF will be dropped.
    - Steps:
      - Load GWAS Catalog data and split multiSNPs into separate rows
      - Load Ensembl VCF and split multiallelic sites into multiple rows
      - For each assoc in gwascat, find best RSID to use from (i) SNP_ID_CURRENT, (ii) SNPS
      - Where RSID starts with "rs", join to VCF on rsid
      - Where RSID starts with "chr", split into chrom and pos then join to VCF on these fields
      - Make variant IDs chr_pos_ref_alt
      - Use "risk allele" to see if it is concordant with alt and ref. Keep those that are concordant or ambiguous (missing)
      - Collapse multiple variants into one row per locus
  3. (Not doing this, see below). Extract beta/OR:
    - Note:
      - I will only keep p-values
      - I have decided not to extract effect sizes/OR and SE. This this can't
        be done in a consistent way across studies, because:
          - Both alleles are not reported so the direction of the effect
            cannot be inferred unambiguously
          - The 95% CI column str in poorly formatted making it very difficult
            to parse
          - Many do not report effects
  4. (Not doing this, see above). Harmonise effects where possible (mark)
  5. Create standard format containing all required columns
    Steps:
      - Remove rows with no mapped `variant_id_b37`
      - Extract p-value mantissa and exponenet
      - Were there are multiple p-values when grouped by (study, variant), keep only the lowest p-value.
    Notes:
      - There can be multiple P-values for each (study, variant). This occurs where different subphenotypes (`P-VALUE (TEXT)`) are reported (e.g. different sexs, subpopulations). Occasionally there are multiple P-values for the same subphenotye also For the time being I will drop these.
      - Question: How to store p-values?

  6. Append UKB top loci from conditional analysis
    Steps:
      - Download from GS
      - Select rows where `locus_index_varid == varid` as these will be the top loci
      - Make study_id from trait name
      - Set any p-values that == 0 to minimum pvalue
      - Remove any rows where nominal p-value is > 1e-5
      - Extract pvalue mantissa and exponent
      - Append

  ? Define locus interval +-1cM. Could be useful for clustering lead SNPs.

#### Study information table
  1. Parse EFOs
  2. N
  3. Number of cases
  4. Ancestry
  5. Whether fine mapping is available
  5. Whether summary statistics

#### Required columns
TODO

#### Usage

```
# Install dependencies into isolated environment
conda env create -n v2d_data --file environment.yaml

# Activate environment
source activate v2d_data

# Alter configuration file
nano config.yaml

# Authenticate google cloud storage
gcloud auth application-default login

# Execute workflow (locally)
snakemake
```
