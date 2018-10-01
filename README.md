Variant-disease tables
======================

This repositroy contains scripts to produce variant-disease association tables for the OT Genetics Portal.

### Contents

  1. Top loci associations table
  2. Study information table
  3. Finemapping (credible set) results table
  4. LD table
  5. Locus overlap table

### Outputs

#### Top loci table

List of loci associated with disease. Currently this data comes from two sources: (i) GWAS Catalog, (ii) Neale et al UK Biobank summary statistics conditional analysis.

Note: GWAS Catalog RSIDs are mapped to variant IDs using the HRC sitelist. Therefore, if the variant is not in HRC it will be dropped from the top loci table.

Columns:
- `study_id`: unique identifier for study
- `variant_id_b37`: chrom_pos_ref_alt (build 37) identifier for variant. RSID to variant ID mapping is non-unique, therefore multiple IDs may exist separated by ';'
- `rsid`: RSID for each corresponding variant ID, separated by ';'
- `pval_mantissa`: the p-value coefficient when written in scientfic notation
- `pval_exponent`: the p-value exponent (base 10)

#### Study table

Information about each study found in the top loci table.

Columns:
  - `study_id`: unique identifier for study
  - `pmid`: pubmed ID (GWAS Catalog studies only)
  - `pub_date`: publication date
  - `pub_journal`: publication journal
  - `pub_title`: publication title
  - `pub_author`: publication author
  - `trait_reported`: trait reported in publication
  - `trait_mapped`: ontology trait label
  - `trait_efos`: EFO short_form
  - `ancestry_initial`: ancestry of initial GWAS sample, separated by ';'
  - `ancestry_replication`: ancestry of replication GWAS sample, separated by ';'
  - `n_initial`: GWAS initial sample size
  - `n_replication`: GWAS replication sample size
  - `n_cases`: number of cases. Warning: there is currently no easy way to get this information from GWAS Catalog, therefore it is set to null

Merges to top loci table using `study_id`.

#### Finemapping table

Finemapping results for linking lead variants (from top loci table) to tagging variants using credible set analysis. Methods: https://github.com/opentargets/finemapping.

Columns:
  - `study_id`: unique identifier for study
  - `index_variantid_b37`: unique variant identifier for index variant, chrom_pos_ref_alt (build 37)
  - `tag_variantid_b37`: unique variant identifier for tagging variant, chrom_pos_ref_alt (build 37)
  - `log10_ABF`: log10 of the approximate Bayes factor for this tagging variant
  - `posterior_prob`: posterior probability of causality for this tagging variant compared to other variants at this locus

Table is pre-filtered to contain only 95% credible set variants. Merges to top loci table using `study_id` and `index_variantid_b37`.

#### LD table

Table of LD values linking index and tagging variants. This information is currently parsed directly from the POSTGAP output but we plan to recalculate these values taking into account the ancestry of the original GWAS.

Columns:
  - `study_id`: unique identifier for study
  - `index_variantid_b37`: unique variant identifier for index variant, chrom_pos_ref_alt (build 37)
  - `tag_variantid_b37`: unique variant identifier for tagging variant, chrom_pos_ref_alt (build 37)
  - `overall_r2`: overall R<sup>2</sup> averaged of superpopulations
  - `AFR_1000G_prop`: proportion of sample from AFR superpopulation
  - `AMR_1000G_prop`: proportion of sample from AMR superpopulation
  - `EAS_1000G_prop`: proportion of sample from EAS superpopulation
  - `EUR_1000G_prop`: proportion of sample from EUR superpopulation
  - `SAS_1000G_prop`: proportion of sample from SAS superpopulation

Table is pre-filtered to only contain R<sup>2</sup> > 0.7. All samples are currently assumed to be European (`EUR_1000G_prop` == 1.0).

#### Locus overlap table

Table showing the number of overlapping tag variants for each (study_id, index_variant) found within 5Mb of each other. This calculated using: (i) only finemapping sets, (ii) only LD sets, (iii) a combination of finemapping and LD, prefering finemapping over LD where available.

Columns:
  - `study_id_A`: locusA unique identifier for study
  - `index_variantid_b37_A`: locusA unique variant identifier, chrom_pos_ref_alt (build 37)
  - `study_id_B`: locusB unique identifier for study
  - `index_variantid_b37_B`: locusB unique variant identifier, chrom_pos_ref_alt (build 37)
  - `set_type`: one of 'finemapping', 'ld' or 'combined' (see above)
  - `distinct_A`: number of tag variants distinct to locusA
  - `overlap_AB`: number of tag variants overlapping locusA and locusB
  - `distinct_B`: number of tag variants distinct to locusB

Only loci with >=1 overlapping tag variants are stored.

### Methods

#### Top loci table
  1. Download all associations from: ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest
  2. Map rsid and chr:pos to variant_ids (GRCh37)
    - This is complicated by:
      - some GWAScat loci have mulitple rsid assignments
      - some rows don't have rsids, use chrom:pos (GRCh37 assumed) instead
      - rsids/chrom:pos can map to multiple variant ids (e.g. rs773154059)
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

  6. Append UKB top loci from conditional analysis
    Steps:
      - Download from GS
      - Select rows where `locus_index_varid == varid` as these will be the top loci
      - Make study_id from trait name
      - Set any p-values that == 0 to minimum float
      - Remove any rows where nominal p-value is > 1e-5
      - Extract pvalue mantissa and exponent
      - Append

  ? Define locus interval +-1cM. Could be useful for clustering lead SNPs. (not doing this, this can be done more easily when loaded into a DB).

#### Study information table

###### GWAS Catalog studies
Get GWAS Catalog studies. This is easier to get from the API as ancestry is
difficult to parse from the flat file.

To extract:
  1. EFO and mapped trait
  2. Reported trait
  2. N
  3. Number of cases (?). This is currently not possible from the API. This field will be set to missing.
  4. Ancestries
  5. Publication info
    - Pubmedid
    - Year
    - Author
    - Title
    - Journal

Notes:
  - There are 18 studies where sample size == 0 because ancestry in API is null. Annalisa has recommended these studies be removed, as they are old curations that are not consistent with the rest of the catalog.

###### Neale UKB studies

Steps
  - Get and merge self-reported and ICD10 EFO curations
  - Map EFOs to mapped_trait
  - Load Nelae et al manifest
  - Merge EFOs with manifest
  - Extract required columns

Todo:
  - Fix field `S72	EFO_000855` in EFO curation (short_form not recognised)

#### Finemapping table

Steps
  - Get from GCS
  - Make study ID
  - Filter to only keep those in 95% credible sets
  - Put into a standard format

Notes:
  - Only keep variants in 95% credible sets

#### LD table

LD calculated for each (study_id, index_variant) pair using 1000 Genomes p3.

Methods:
  1. GWAS Catalog populations mapped to 1000 Genomes super populations
    - Manually curated. Suppl fig 6 from here https://arxiv.org/pdf/1805.03233.pdf is useful.
    - Curations in `configs/gwascat_superpopulation_lut.curated.tsv`
  2. Output list of index variants with super populaiton proportions for each study
  3. Use plink to calculate correaltion coefficients (R) for each index study in each 1000 superpopulation.
  4. Fisher-Z transform R coefficients.
  5. For each study take weighted average across populations weighting by sample size.
  6. Inverse-Fisher-Z transform to get R values
  7. Take R-squared and filter rows R2 < 0.7.

Notes:
  - Studies with missing or NR ancestries will be assumed to be European


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
version_date=`date +%y%m%d`
cores=3
snakemake -s 1_make_tables.Snakefile --config version=$version_date --cores $cores
snakemake -s 2_calculate_LD_table.Snakefile --config version=$version_date --cores $cores
snakemake -s 3_make_overlap_table.Snakefile --config version=$version_date --cores $cores

# Dry-run gsutil rsync to copy from staging to live
gsutil -m rsync -rn gs://genetics-portal-staging/v2d/180904 gs://genetics-portal-data/v2d

# Temp
#snakemake -s 1_make_tables.Snakefile output/$version_date/ld_analysis_input.tsv.gz --cores $cores
#snakemake -s 1_make_tables.Snakefile tmp/$version_date/ld/variantID_to_1000Gp3_lut.tsv.gz --cores $cores

```

### Test

```
version_date=180928
cores=3
echo snakemake -s 1_make_tables.Snakefile --config version=$version_date --cores $cores -np
echo snakemake -s 3_make_overlap_table.Snakefile --config version=$version_date --cores $cores -np
```
