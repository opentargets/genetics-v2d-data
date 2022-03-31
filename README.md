Variant-disease tables
======================

This repository contains scripts to produce variant-to-disease (V2D) association tables for Open Targets Genetics.

Changes made (Jan 2019):
- All outputs in Apache Parquet format
- Top loci table contains harmonised effect size, 95% CI and direction
- GWAS Catalog sub-phenotypes (in `P-VALUE (TEXT)` column) are split into separate `study_ids`
- All `variant_id`s are decomposed into `chrom`, `pos`, `ref`, `alt`
- GWAS Catalog curated associations are clustered to remove non-independent loci

### Contents

- [Usage](#usage)
- Tables produced by workflow
  1. [Top loci associations table](#top-loci-table)
  2. [Study information table](#study-table)
  3. [Finemapping (credible set) results table](#finemapping-table)
  4. [LD table](#ld-table)
  5. [Locus overlap table](#locus-overlap-table)

### Usage

```bash
# Set parameters.
export INSTANCE_NAME=v2d_data
export INSTANCE_ZONE=europe-west1-d

# Create the instance and SSH.
gcloud compute instances create \
  ${INSTANCE_NAME} \
  --project=open-targets-genetics-dev \
  --zone=${INSTANCE_ZONE} \
  --machine-type=n1-highmem-32 \
  --service-account=426265110888-compute@developer.gserviceaccount.com \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --create-disk=auto-delete=yes,boot=yes,device-name=${INSTANCE_NAME},image=projects/ubuntu-os-cloud/global/images/ubuntu-2004-focal-v20210927,mode=rw,size=2000,type=projects/open-targets-eu-dev/zones/europe-west1-d/diskTypes/pd-balanced
gcloud compute ssh --zone ${INSTANCE_ZONE} ${INSTANCE_NAME}

# Setup on gcloud if needed
bash setup_gcloud.sh

# Authenticate google cloud storage if needed
gcloud auth application-default login

# Set up the instance.
sudo apt update
sudo apt install -yf \
  openjdk-13-jre-headless \
  python3-pip \
  jq
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

# Install dependencies into isolated environment
conda env create -n v2d_data --file environment.yaml
conda activate v2d_data

# Alter configuration file
nano config.yaml

# Remove downloads from any previous runs
rm -r www.ebi.ac.uk/gwas/

# Execute workflows
# May want to start tmux session to avoid problems if connection is lost
# (I've gotten snakemake problems on subsequent attempts when this happens too)
tmux


# May want to use a smaller machine for step 1, then scale up to more
# cores for step 2, and back down to a small machine for step 3
export PYSPARK_SUBMIT_ARGS="--driver-memory 100g pyspark-shell"

export VERSION_DATE=`date +%y%m%d`
mkdir -p logs/$VERSION_DATE

# Run workflow
time snakemake -s 1_make_tables.Snakefile --config version=$VERSION_DATE --cores all | tee logs/$VERSION_DATE/1_make_tables.log 2>&1 # Takes a couple hours
time snakemake -s 2_calculate_LD_table.Snakefile --config version=$VERSION_DATE --cores all 2>&1 | tee logs/$VERSION_DATE/2_calculate_LD_table.log 2>&1 # Takes ~7 hrs on 31 cores

# Reduce machine size in Google VM instance
# This step only uses 1 core actually - but I'm not sure how much memory
cores=3
time snakemake -s 3_make_overlap_table.Snakefile --config version=$VERSION_DATE --cores all 2>&1 | tee logs/$VERSION_DATE/3_make_overlap_table.log 2>&1 # Takes a couple hours

# Upload output dir to google cloud storage
gsutil -m rsync -r output/$VERSION_DATE gs://genetics-portal-dev-staging/v2d/$VERSION_DATE
```

### Tables

#### Top loci table

List of loci associated with disease. Currently this data comes from two sources: (i) [GWAS Catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads), (ii) [Neale *et al.* UK Biobank summary statistics (version 1)](http://www.nealelab.is/uk-biobank).

##### Parquet meta info

```
file schema:    schema 
--------------------------------------------------------------------------------
study_id:       OPTIONAL BINARY L:STRING R:0 D:1
chrom:          OPTIONAL BINARY L:STRING R:0 D:1
pos:            OPTIONAL INT64 R:0 D:1
ref:            OPTIONAL BINARY L:STRING R:0 D:1
alt:            OPTIONAL BINARY L:STRING R:0 D:1
rsid:           OPTIONAL BINARY L:STRING R:0 D:1
direction:      OPTIONAL BINARY L:STRING R:0 D:1
beta:           OPTIONAL DOUBLE R:0 D:1
beta_ci_lower:  OPTIONAL DOUBLE R:0 D:1
beta_ci_upper:  OPTIONAL DOUBLE R:0 D:1
odds_ratio:     OPTIONAL DOUBLE R:0 D:1
oddsr_ci_lower: OPTIONAL DOUBLE R:0 D:1
oddsr_ci_upper: OPTIONAL DOUBLE R:0 D:1
pval_mantissa:  OPTIONAL DOUBLE R:0 D:1
pval_exponent:  OPTIONAL INT64 R:0 D:1
```

##### Top loci columns
- `study_id`: unique identifier for study
- `variant_id_b37`: chrom_pos_ref_alt (build 37) identifier for variant. RSID to variant ID mapping is non-unique, therefore multiple IDs may exist separated by ';'
- `rsid`: RSID for each corresponding variant ID, separated by ';'
- `direction` ([+, -, null]): direction of the effect wrt to the alt allele
- `beta` (float, nullable): beta effect size for quantitative trait
- `oddsr` (float, nullable): OR effect size for binary trait
- `ci_lower` (float, nullable): Lower 95% confidence interval of effect estimate
- `ci_upper` (float, nullable): Upper 95% confidence interval of effect estimate
- `pval_mantissa`: the p-value coefficient when written in scientfic notation
- `pval_exponent`: the p-value exponent (base 10)

##### Top loci methods


1. Download [all associations](ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest)
2. Map rsIDs and chr:pos fields to variant_ids (GRCh37)
  - This is complicated by the fact that:
    - some GWAScat loci have mulitple rsid assignments
    - some rows don't have rsids, will use chrom:pos (GRCh37 assumed) instead
    - rsIDs/chrom:pos can map to multiple variant ids (e.g. rs773154059)
  - Steps:
    - Load GWAS Catalog data and explode multiSNPs into separate rows
    - Load HRC site list + Ensembl VCF and split multiallelic sites into multiple rows
    - For each assoc in gwascat, extract rsID from SNP_ID_CURRENT or SNPS
    - Variants are first joined to HRC sitelist VCF, then unmerged rows are joined to the Ensembl VCF. Join on:
      1. Where RSID starts with "rs", join to VCF on rsid
      2. Where RSID starts with "chr", split into chrom and pos then join to VCF on these fields
    - Make variant IDs in the format chr_pos_ref_alt
    - Use "risk allele" from GWAS Catalog to see if it is concordant with alt and ref. Keep those that are concordant or ambiguous (missing), remove those that are not concordant.
    - Identify studies with an abnormally large number of reported loci compared to number of loci after distance-based clustering (±500kb). Criteria for "abnormal" are set to number of loci > 10, and decrease in loci count after clustering of >10%. For "abnormal" studies (N≈120), apply distance based clustering.
  - Notes:
    - Two associations don't have the same number of SNPs in SNPS, CHR_ID and
      CHR_POS. These have been dropped as they are ambiguous.
    - I am assuming that variants with id e.g. chr21:43693789 are from build GRCh37.
    - Associations that do not have an rsID or chr:pos will be dropped
    - Associations that do not have a match in the Ensembl VCF will be dropped.
3. Create standard format containing all required columns
  Steps:
    - Remove rows with no mapped `variant_id_b37`
    - Extract p-value mantissa and exponenet
4. Append UKB top loci from the finemapping conditional analysis
  Steps:
    - Download finemapping results from GCS (`gs://genetics-portal-input/uk_biobank_analysis/em21/neale_summary_statistics_20170915/finemapping/results/neale_ukb_uk10kref_180502.crediblesets.long.varids.tsv.gz`)
    - Select rows where `locus_index_varid == varid` as these will be the top loci
    - Make study_id from trait name (prepend "NEALEUKB_")
    - Set any p-values that == 0 to minimum float
    - Remove any rows where nominal p-value is > 1e-5
    - Extract pvalue mantissa and exponent.

###### Top loci notes
- GWAS Catalog RSIDs are mapped to variant IDs using (i) HRC sitelist VCF, (ii) Ensembl VCF. Therefore, if the variant is not in either of these it will be dropped from the top loci table.
- There can be multiple P-values for each (study, variant). This occurs where different subphenotypes (`P-VALUE (TEXT)`) are reported (e.g. different sexs, subpopulations). Occasionally there are multiple P-values for the same subphenotye also For the time being I will drop these.
- We have decided not to extract effect sizes/OR and SE. This this can't be done in a consistent way across studies, because:
  - Both alleles are not reported so the direction of the effect
    cannot be inferred unambiguously
  - The 95% CI column str in poorly formatted making it very difficult to parse
  - Many do not report effects

#### Study table

Information about each study found in the top loci table.

##### Parquet meta info
```
file schema:          schema 
--------------------------------------------------------------------------------
study_id:             OPTIONAL BINARY L:STRING R:0 D:1
pmid:                 OPTIONAL BINARY L:STRING R:0 D:1
pub_date:             OPTIONAL BINARY L:STRING R:0 D:1
pub_journal:          OPTIONAL BINARY L:STRING R:0 D:1
pub_title:            OPTIONAL BINARY L:STRING R:0 D:1
pub_author:           OPTIONAL BINARY L:STRING R:0 D:1
trait_reported:       OPTIONAL BINARY L:STRING R:0 D:1
trait_efos:           OPTIONAL F:1 
.list:                REPEATED F:1 
..item:               OPTIONAL BINARY L:STRING R:1 D:3
ancestry_initial:     OPTIONAL F:1 
.list:                REPEATED F:1 
..item:               OPTIONAL BINARY L:STRING R:1 D:3
ancestry_replication: OPTIONAL F:1 
.list:                REPEATED F:1 
..item:               OPTIONAL BINARY L:STRING R:1 D:3
n_initial:            OPTIONAL INT64 R:0 D:1
n_replication:        OPTIONAL INT64 R:0 D:1
n_cases:              OPTIONAL INT64 R:0 D:1
trait_category:       OPTIONAL BINARY L:STRING R:0 D:1
num_assoc_loci:       OPTIONAL INT64 R:0 D:1
```

##### Study table columns
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
  - `num_assoc_loci` (int): total number of associated loci for this study in the top loci table

Merges to top loci table using `study_id`.

##### Study table methods

###### GWAS Catalog studies
Steps:
1. Get [GWAS Catalog study and ancestry information](https://www.ebi.ac.uk/gwas/docs/file-downloads).
2. Combined ancestry information into a single field.
3. Make trait_code hash using 16 characters from `hasher.hexdigest()[0:16]`
4. Remove rows with no sample size information (see note below)
5. For studies with multiple different values in `P-VALUE (TEXT)` (multi-trait studies), split each phenotype into its own study ID.

Notes:
  - There are 18 studies where sample size == 0 because ancestry in API is null. Annalisa has recommended these studies be removed, as they are old curations that are not consistent with the rest of the catalog.

###### Neale UKB studies
Steps
  - Get and merge self-reported and ICD10 manual EFO curations
  - Map EFOs to mapped_trait
  - Load Neale et al manifest
  - Merge EFOs with manifest
  - Extract required columns
  - Create trait code using field code from Neale et al.
  - Add trait categories field, acquired from [Oxford BIG server](http://big.stats.ox.ac.uk/)

Todo:
  - Fix field `S72	EFO_000855` in EFO curation (short_form not recognised)

#### Finemapping table

Credible set analysis results used to link index variants to tag variants. Full finemapping methods can be seen here: https://github.com/opentargets/finemapping

##### Parquet meta info

```
file schema:    schema 
--------------------------------------------------------------------------------
study_id:       OPTIONAL BINARY L:STRING R:0 D:1
lead_chrom:     OPTIONAL BINARY L:STRING R:0 D:1
lead_pos:       OPTIONAL INT64 R:0 D:1
lead_ref:       OPTIONAL BINARY L:STRING R:0 D:1
lead_alt:       OPTIONAL BINARY L:STRING R:0 D:1
tag_chrom:      OPTIONAL BINARY L:STRING R:0 D:1
tag_pos:        OPTIONAL INT64 R:0 D:1
tag_ref:        OPTIONAL BINARY L:STRING R:0 D:1
tag_alt:        OPTIONAL BINARY L:STRING R:0 D:1
log10_ABF:      OPTIONAL DOUBLE R:0 D:1
posterior_prob: OPTIONAL DOUBLE R:0 D:1
```

##### Finemapping table columns
  - `study_id`: unique identifier for study
  - `index_variantid_b37`: unique variant identifier for index variant, chrom_pos_ref_alt (build 37)
  - `tag_variantid_b37`: unique variant identifier for tagging variant, chrom_pos_ref_alt (build 37)
  - `log10_ABF`: log10 of the approximate Bayes factor for this tagging variant
  - `posterior_prob`: posterior probability of causality for this tagging variant compared to other variants at this locus

Table is pre-filtered to contain only 95% credible set variants. Merges to top loci table using `study_id` and `index_variantid_b37`.

##### Finemapping table methods

Steps
  - Downlaod from from GCS
  - Make study ID frrom the field code
  - Filter to only keep those in 95% credible sets
  - Put into a standardised format (above)

#### LD table

Table of LD values linking index varaints to tag variants. I

##### Parquet meta info

```
root
 |-- study_id: string (nullable = true)
 |-- lead_chrom: string (nullable = true)
 |-- lead_pos: integer (nullable = true)
 |-- lead_ref: string (nullable = true)
 |-- lead_alt: string (nullable = true)
 |-- tag_chrom: string (nullable = true)
 |-- tag_pos: integer (nullable = true)
 |-- tag_ref: string (nullable = true)
 |-- tag_alt: string (nullable = true)
 |-- overall_r2: double (nullable = true)
 |-- pics_mu: double (nullable = true)
 |-- pics_postprob: double (nullable = true)
 |-- pics_95perc_credset: boolean (nullable = true)
 |-- pics_99perc_credset: boolean (nullable = true)
 |-- AFR_1000G_prop: double (nullable = true)
 |-- AMR_1000G_prop: double (nullable = true)
 |-- EAS_1000G_prop: double (nullable = true)
 |-- EUR_1000G_prop: double (nullable = true)
 |-- SAS_1000G_prop: double (nullable = true)
```

##### LD table columns
  - `study_id`: unique identifier for study
  - `index_variantid_b37`: unique variant identifier for index variant, chrom_pos_ref_alt (build 37)
  - `tag_variantid_b37`: unique variant identifier for tagging variant, chrom_pos_ref_alt (build 37)
  - `overall_r2`: overall R<sup>2</sup> averaged of superpopulations
  - `pics_mu`: Mu statistic from PICS calculation
  - `pics_postprob`: PICS posterior probability that this variant is causal
  - `pics_95perc_credset`: If this variant is in 95% credible set
  - `pics_99perc_credset`: If this variant is in 99% credible set
  - `AFR_1000G_prop`: proportion of sample from AFR superpopulation
  - `AMR_1000G_prop`: proportion of sample from AMR superpopulation
  - `EAS_1000G_prop`: proportion of sample from EAS superpopulation
  - `EUR_1000G_prop`: proportion of sample from EUR superpopulation
  - `SAS_1000G_prop`: proportion of sample from SAS superpopulation

Table is pre-filtered to only contain R<sup>2</sup> > 0.5.

##### LD table methods

LD calculated for each (study_id, index_variant) pair using 1000 Genomes p3.

Methods:
  1. GWAS Catalog populations mapped to 1000 Genomes super populations
    - Manually curated. Suppl fig 6 from here https://arxiv.org/pdf/1805.03233.pdf is useful.
    - Curations in `configs/gwascat_superpopulation_lut.curated.tsv`
  2. Output list of index variants with super population proportions for each study (`ld_analysis_input.tsv.gz`). Proportions are calulcated across both initial and replication samples. If a compound ancestry is provided by GWAS Catalog e.g. "popA, popB=100", then the sample size is split equally amoungst them, i.e. popA=50, popB=50. This file is used as the LD analysis input manifest.
  3. Prepare 1000 Genomes haplotypes. Script [available here](https://github.com/opentargets/genetics-backend/tree/master/reference_data/1000Genomes_phase3):
    1. Download
    2. Normalise to reference genome
    3. Extract samples separately for each super population
    4. Filter to remove MAF < 1% and genotyping-rate < 5%
    5. Convert to bed, bim, fam
  3. Use plink to calculate correlation coefficients (R) for each index study in each 1000G superpopulation.
  4. Merge index variant correlation tables to input manifest
  4. Fisher-Z transform R coefficients.
  5. For each study take weighted average across populations weighting by sample size.
  6. Inverse-Fisher-Z transform to get R values
  7. Conduct [PICS finemapping analysis](https://www.ncbi.nlm.nih.gov/pubmed/25363779)

Notes:
  - Studies with missing or NR ancestries will be assumed to be European

#### Locus overlap table

Table containing metrics on the degree of (tag variant) overlap between index variants with 5Mb of each other.

##### Locus overlap table columns
  - `study_id_A`: locusA unique identifier for study
  - `index_variantid_b37_A`: locusA unique variant identifier, chrom_pos_ref_alt (build 37)
  - `study_id_B`: locusB unique identifier for study
  - `index_variantid_b37_B`: locusB unique variant identifier, chrom_pos_ref_alt (build 37)
  - `set_type`: one of 'finemapping', 'ld' or 'combined' (see below)
  - `distinct_A`: number of tag variants distinct to locusA
  - `overlap_AB`: number of tag variants overlapping locusA and locusB
  - `distinct_B`: number of tag variants distinct to locusB

Only loci with >=1 overlapping tag variants are stored.

##### Locus overlap methods

Table showing the number of overlapping tag variants for each (study_id, index_variant) found within 5Mb of each other. This calculated using: (i) only finemapping sets, (ii) only LD sets, (iii) a combination of finemapping and LD, prefering finemapping over LD where available. The 'combined' set is used for the Genetics Portal.


##### Effect directions to check in release
- GCST006612 1_55505647_G_T effect allele=T -0.325
- GCST002898 1_55505647_G_T effect allele=T -0.53
- GCST005194_1 1_55505647_G_T effect allele=T -0.282