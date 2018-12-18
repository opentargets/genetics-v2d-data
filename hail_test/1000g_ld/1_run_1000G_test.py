#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import hail as hl
import os
import sys
import pandas as pd

def main():

    # Args
    min_maf = 0.01
    input_manifest = '/Users/em21/Projects/ot_genetics/genetics-v2d_data/output/181217/ld_analysis_input.tsv.gz'
    window = 500 #kb
    reference_genome = 'GRCh37'

    #
    # Read data ----------------------------------------------------------------
    #

    # Download 1000G example data
    os.makedirs('data', exist_ok=True)
    hl.utils.get_1kg('data')

    # Read the dataset (droppoing duplicate rows)
    mt = ( hl.read_matrix_table('data/1kg.mt')
             .distinct_by_row() )

    #
    # Merge phenotype (superpopulation) annotation -----------------------------
    #

    # Load annotations
    table = (hl.import_table('data/1kg_annotations.txt', impute=True)
               .key_by('Sample'))

    # Merge annotations
    mt = mt.annotate_cols(pheno = table[mt.s])

    #
    # Filter based on intervals ------------------------------------------------
    #

    # Load manifest
    manifest = (pd.read_csv(input_manifest, sep='\t', header=0)
                  .sort_values(['chrom', 'pos', 'ref', 'alt']))

    # Write temp interval file
    interval_file = 'tmp/interval_file.tsv'
    make_interval_file(manifest[['chrom', 'pos']].drop_duplicates().values.tolist(),
                       interval_file,
                       reference_genome=reference_genome,
                       window=window)

    # Load intervals into hail and filter
    interval_table = hl.import_locus_intervals(interval_file,
        reference_genome='GRCh37')
    print('Num variants before interval filter: {0}'.format(mt.count_rows()))
    mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]))
    print('Num variants after interval filter: {0}'.format(mt.count_rows()))

    #
    # Perform QC ---------------------------------------------------------------
    #

    # Calc sample QC
    mt = hl.sample_qc(mt)

    # Apply sample filters
    mt = mt.filter_cols((mt.sample_qc.dp_stats.mean >= 4) & (mt.sample_qc.call_rate >= 0.97))
    print('After filter, %d/284 samples remain.' % mt.count_cols())

    # Apply genotype filter (recommended in hail docs)
    ab = mt.AD[1] / hl.sum(mt.AD)
    filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                           (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                           (mt.GT.is_hom_var() & (ab >= 0.9)))
    mt = mt.filter_entries(filter_condition_ab)

    # Calc variant QC stats
    mt = hl.variant_qc(mt)

    # Apply row filters
    mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)
    mt = mt.filter_rows(mt.variant_qc.AF[1] > min_maf)
    print('Samples: %d  Variants: %d' % (mt.count_cols(), mt.count_rows()))


    #
    # Get variant row indexes --------------------------------------------------
    #

    # Make variant key table using variants in the manifest
    var_id = (
        manifest.loc[:, ['chrom', 'pos', 'ref', 'alt']]
                .drop_duplicates()
                .apply(lambda row: ':'.join([str(x)
                                             for x
                                             in row.values.tolist()]), axis=1) )
    data = [{'v': entry} for entry in var_id.values.tolist()]
    var_table = hl.Table.parallelize(data, hl.dtype('struct{v: str}'))
    var_table = var_table.transmute(**hl.parse_variant(var_table.v))
    var_table = var_table.key_by('locus', 'alleles')

    # Inner merge with mt, to get a list of variants in ld dataset
    print('Total number of unique variants in manifest: ', var_table.count())
    var_table = var_table.join(mt.rows().select(), how='inner')
    print('Number of variant in 1000G: ', var_table.count())

    # Add row_index to the variant table
    mt = mt.add_row_index()
    var_table = var_table.annotate(
        mt_var_index = mt.index_rows(var_table.locus, var_table.alleles).row_idx
    )

    # DEBUG - e.g. row 5420
    var_table.export('tmp/var_table.table.tsv')
    mt.rows().export('tmp/mt.table.tsv')

    # Extract row vertices
    var_indexes = var_table.mt_var_index.collect()

    #
    # Calc LD ------------------------------------------------------------------
    #

    print('Calculating LD...')
    ld = hl.ld_matrix(mt.GT.n_alt_alleles(), mt.locus, radius=window*1000)

    # Filter rows to only return those in the manifest
    print(list(var_indexes))
    ld_filt = ld.filter_rows(var_indexes)

    print('Full LD matrix shape: ', ld.shape)
    print('Filtered LD matrix shape: ', ld_filt.shape)

    return 0

def make_interval_file(locus_list, outf, reference_genome, window=500):
    ''' Writes a file of intervals
    Args:
        locus_list (list): list of (chrom, pos) tuples
        outf (str): output file
        window (int): kb window to create locus
    Returns:
        None
    '''
    # Make out dir
    os.makedirs(os.path.dirname(outf), exist_ok=True)

    # Get chrom lengths
    if reference_genome == 'GRCh37':
        chrom_lengths = {'1': 249250621,
                         '10': 135534747,
                         '11': 135006516,
                         '12': 133851895,
                         '13': 115169878,
                         '14': 107349540,
                         '15': 102531392,
                         '16': 90354753,
                         '17': 81195210,
                         '18': 78077248,
                         '19': 59128983,
                         '2': 243199373,
                         '20': 63025520,
                         '21': 48129895,
                         '22': 51304566,
                         '3': 198022430,
                         '4': 191154276,
                         '5': 180915260,
                         '6': 171115067,
                         '7': 159138663,
                         '8': 146364022,
                         '9': 141213431,
                         'X': 155270560,
                         'Y': 59373566}
    else:
        sys.exit('Error: only supports GRCh37 currently')

    # Write file
    with open(outf, 'w') as out_h:
        for chrom, pos in locus_list:
            start = max(int(pos) - window * 1000, 1)
            end = min(int(pos) + window * 1000, chrom_lengths[chrom])
            out_row = [chrom, start, end]
            out_h.write('\t'.join([str(x) for x in out_row]) + '\n')

    return 0


if __name__ == '__main__':

    main()
