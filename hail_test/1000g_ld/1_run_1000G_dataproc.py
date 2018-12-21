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
    input_manifest = 'gs://genetics-portal-data/v2d/extras/ld_analysis_input.tsv.gz'
    out_path = 'gs://genetics-portal-analysis/hail_ld_test/ld_filt.table.tsv.gz'
    reference_genome = 'GRCh37'
    min_r2 = 0.7
    min_maf = 0.01
    window = 500 #kb

    #
    # Read data ----------------------------------------------------------------
    #

    # Load
    mt = hl.experimental.load_dataset(
        name='1000_genomes',
        version='phase3',
        reference_genome='GRCh37'
    )

    # Remove duplicated
    mt = mt.distinct_by_row()

    # Downsample
    mt = mt.sample_rows(0.01)

    # # Download 1000G example data
    # os.makedirs('data', exist_ok=True)
    # hl.utils.get_1kg('data')
    #
    # # Read the dataset (droppoing duplicate rows)
    # mt = ( hl.read_matrix_table('data/1kg.mt')
    #          .distinct_by_row()
    #          # .sample_rows(0.01) # DEBUG
    # )
    #
    # #
    # # Merge phenotype (superpopulation) annotation -----------------------------
    # #
    #
    # # Load annotations
    # table = (hl.import_table('data/1kg_annotations.txt', impute=True)
    #            .key_by('Sample'))
    #
    # # Merge annotations
    # mt = mt.annotate_cols(pheno = table[mt.s])

    #
    # Filter based on intervals ------------------------------------------------
    #

    # Load manifest
    with hl.utils.hadoop_open(input_manifest, 'r') as in_h:
        manifest = (pd.read_csv(in_h, sep='\t', header=0)
                      .sort_values(['chrom', 'pos', 'ref', 'alt']))

    # Write temp interval file
    interval_file = 'file:///tmp/interval_file.tsv' # ERROR here
    make_interval_file(manifest[['chrom', 'pos']].drop_duplicates()
                                                 .values
                                                 .tolist(),
                       interval_file,
                       reference_genome=reference_genome,
                       window=window)

    # Load intervals into hail and filter
    interval_table = hl.import_locus_intervals(interval_file,
        reference_genome='GRCh37')
    # interval_table = interval_table.order_by(hl.asc('interval')).key_by('interval')
    print('Num variants before interval filter: {0}'.format(mt.count_rows()))
    mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]))
    print('Num variants after interval filter: {0}'.format(mt.count_rows()))

    #
    # Perform QC ---------------------------------------------------------------
    #

    # Calc sample QC
    mt = hl.sample_qc(mt)

    # Apply sample filters
    mt = mt.filter_cols((mt.sample_qc.dp_stats.mean >= 4) &
                        (mt.sample_qc.call_rate >= 0.97))
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
    # var_table.describe()
    # sys.exit()

    # DEBUG - e.g. row 5420
    # var_table.export('tmp/var_table.table.tsv')
    # mt.rows().export('tmp/mt.table.tsv')

    # Extract row vertices
    var_indexes = var_table.mt_var_index.collect()

    #
    # Calc LD ------------------------------------------------------------------
    #

    # Perform LD calculation on each superpopulation separately
    ld_tables = []
    superpops = list(mt.aggregate_cols(
        hl.agg.collect_as_set(mt.pheno.SuperPopulation)))
    for superpop in superpops:
        print('Calculating LD for: ', superpop)
        # Filter to keep only columns matching the superpop
        mt_pop = mt.filter_cols(mt.pheno.SuperPopulation == superpop)
        # Calculate ld
        ld_pop = ( hl.ld_matrix(mt_pop.GT.n_alt_alleles(),
                                mt_pop.locus,
                                radius=window*1000)
                     .filter_rows(var_indexes)
                     .entries()
                     .rename({'entry': superpop}) )
        ld_tables.append(ld_pop)

    # Join all superpopulations together
    ld_merged = ld_tables[0]
    for ld_table in ld_tables[1:]:
        ld_merged = ld_merged.join(ld_table)

    # Filter to keep any with
    ld_filt = ld_merged.filter(
        (ld_merged['AFR'] ** 2 > min_r2) |
        (ld_merged['EUR'] ** 2 > min_r2) |
        (ld_merged['SAS'] ** 2 > min_r2) |
        (ld_merged['EAS'] ** 2 > min_r2) |
        (ld_merged['AMR'] ** 2 > min_r2)
    )
    print('LD filter rows: ', ld_filt.count())

    # Annotate i with names from var_table
    i_anno = var_table.add_index().key_by('idx')
    ld_filt = ld_filt.annotate(
        i_locus = i_anno[ld_filt.i].locus,
        i_alleles = i_anno[ld_filt.i].alleles
    )

    # Annotate j with names from mt
    j_anno = mt.rows().key_by('row_idx')
    ld_filt = ld_filt.annotate(
        j_locus = j_anno[ld_filt.j].locus,
        j_alleles = j_anno[ld_filt.j].alleles
    )

    # DEBUG
    print('Writing...')
    ld_filt.export(out_path)
    # ld_filt.export('tmp/ld_filt.table.tsv.gz')

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


#   File "/tmp/bf7daf429b9b47b48852b57b9d1df7cf/1_run_1000G_dataproc.py", line 261, in <module>
#     main()
#   File "/tmp/bf7daf429b9b47b48852b57b9d1df7cf/1_run_1000G_dataproc.py", line 84, in main
#     print('Num variants after interval filter: {0}'.format(mt.count_rows()))
#   File "/tmp/bf7daf429b9b47b48852b57b9d1df7cf/hail-0.2-07b91f4bd378.zip/hail/matrixtable.py", line 2092, in count_rows
#   File "/tmp/bf7daf429b9b47b48852b57b9d1df7cf/hail-0.2-07b91f4bd378.zip/hail/backend/backend.py", line 36, in interpret
#   File "/usr/lib/spark/python/lib/py4j-0.10.4-src.zip/py4j/java_gateway.py", line 1133, in __call__
#   File "/tmp/bf7daf429b9b47b48852b57b9d1df7cf/hail-0.2-07b91f4bd378.zip/hail/utils/java.py", line 210, in deco
# hail.utils.java.FatalError: FileNotFoundException: File file:/tmp/interval_file.tsv does not exist
#
# Java stack trace:
# org.apache.spark.SparkException: Job aborted due to stage failure: Task 0 in stage 1.0 failed 20 times, most recent failure: Lost task 0.19 in stage 1.0 (TID 320, em-cluster-w-0.c.open-targets-genetics.internal, executor 4): java.io.FileNotFoundException: File file:/tmp/interval_file.tsv does not exist
# 	at org.apache.hadoop.fs.RawLocalFileSystem.deprecatedGetFileStatus(RawLocalFileSystem.java:635)
# 	at org.apache.hadoop.fs.RawLocalFileSystem.getFileLinkStatusInternal(RawLocalFileSystem.java:861)
# 	at org.apache.hadoop.fs.RawLocalFileSystem.getFileStatus(RawLocalFileSystem.java:625)
# 	at org.apache.hadoop.fs.FilterFileSystem.getFileStatus(FilterFileSystem.java:442)
# 	at org.apache.hadoop.fs.ChecksumFileSystem$ChecksumFSInputChecker.<init>(ChecksumFileSystem.java:146)
# 	at org.apache.hadoop.fs.ChecksumFileSystem.open(ChecksumFileSystem.java:347)
# 	at org.apache.hadoop.fs.FileSystem.open(FileSystem.java:787)
# 	at org.apache.hadoop.mapred.LineRecordReader.<init>(LineRecordReader.java:109)
# 	at org.apache.hadoop.mapred.TextInputFormat.getRecordReader(TextInputFormat.java:67)
# 	at org.apache.spark.rdd.HadoopRDD$$anon$1.liftedTree1$1(HadoopRDD.scala:251)
# 	at org.apache.spark.rdd.HadoopRDD$$anon$1.<init>(HadoopRDD.scala:250)
# 	at org.apache.spark.rdd.HadoopRDD.compute(HadoopRDD.scala:208)
# 	at org.apache.spark.rdd.HadoopRDD.compute(HadoopRDD.scala:94)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.scheduler.ResultTask.runTask(ResultTask.scala:87)
# 	at org.apache.spark.scheduler.Task.run(Task.scala:108)
# 	at org.apache.spark.executor.Executor$TaskRunner.run(Executor.scala:338)
# 	at java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1149)
# 	at java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:624)
# 	at java.lang.Thread.run(Thread.java:748)
#
# Driver stacktrace:
# 	at org.apache.spark.scheduler.DAGScheduler.org$apache$spark$scheduler$DAGScheduler$$failJobAndIndependentStages(DAGScheduler.scala:1517)
# 	at org.apache.spark.scheduler.DAGScheduler$$anonfun$abortStage$1.apply(DAGScheduler.scala:1505)
# 	at org.apache.spark.scheduler.DAGScheduler$$anonfun$abortStage$1.apply(DAGScheduler.scala:1504)
# 	at scala.collection.mutable.ResizableArray$class.foreach(ResizableArray.scala:59)
# 	at scala.collection.mutable.ArrayBuffer.foreach(ArrayBuffer.scala:48)
# 	at org.apache.spark.scheduler.DAGScheduler.abortStage(DAGScheduler.scala:1504)
# 	at org.apache.spark.scheduler.DAGScheduler$$anonfun$handleTaskSetFailed$1.apply(DAGScheduler.scala:814)
# 	at org.apache.spark.scheduler.DAGScheduler$$anonfun$handleTaskSetFailed$1.apply(DAGScheduler.scala:814)
# 	at scala.Option.foreach(Option.scala:257)
# 	at org.apache.spark.scheduler.DAGScheduler.handleTaskSetFailed(DAGScheduler.scala:814)
# 	at org.apache.spark.scheduler.DAGSchedulerEventProcessLoop.doOnReceive(DAGScheduler.scala:1732)
# 	at org.apache.spark.scheduler.DAGSchedulerEventProcessLoop.onReceive(DAGScheduler.scala:1687)
# 	at org.apache.spark.scheduler.DAGSchedulerEventProcessLoop.onReceive(DAGScheduler.scala:1676)
# 	at org.apache.spark.util.EventLoop$$anon$1.run(EventLoop.scala:48)
# 	at org.apache.spark.scheduler.DAGScheduler.runJob(DAGScheduler.scala:630)
# 	at org.apache.spark.SparkContext.runJob(SparkContext.scala:2029)
# 	at org.apache.spark.SparkContext.runJob(SparkContext.scala:2050)
# 	at org.apache.spark.SparkContext.runJob(SparkContext.scala:2069)
# 	at org.apache.spark.SparkContext.runJob(SparkContext.scala:2094)
# 	at org.apache.spark.rdd.RDD$$anonfun$collect$1.apply(RDD.scala:936)
# 	at org.apache.spark.rdd.RDDOperationScope$.withScope(RDDOperationScope.scala:151)
# 	at org.apache.spark.rdd.RDDOperationScope$.withScope(RDDOperationScope.scala:112)
# 	at org.apache.spark.rdd.RDD.withScope(RDD.scala:362)
# 	at org.apache.spark.rdd.RDD.collect(RDD.scala:935)
# 	at is.hail.sparkextras.ContextRDD.collect(ContextRDD.scala:153)
# 	at is.hail.rvd.RVD$.getKeyInfo(RVD.scala:1061)
# 	at is.hail.rvd.RVD$.makeCoercer(RVD.scala:1125)
# 	at is.hail.rvd.RVD$.coerce(RVD.scala:1083)
# 	at is.hail.rvd.RVD.changeKey(RVD.scala:121)
# 	at is.hail.rvd.RVD.changeKey(RVD.scala:118)
# 	at is.hail.rvd.RVD.enforceKey(RVD.scala:113)
# 	at is.hail.expr.ir.TableKeyBy.execute(TableIR.scala:227)
# 	at is.hail.expr.ir.MatrixAnnotateRowsTable.execute(MatrixIR.scala:1725)
# 	at is.hail.expr.ir.CastMatrixToTable.execute(TableIR.scala:1293)
# 	at is.hail.expr.ir.TableFilter.execute(TableIR.scala:301)
# 	at is.hail.expr.ir.TableMapRows.execute(TableIR.scala:625)
# 	at is.hail.expr.ir.CastTableToMatrix.execute(MatrixIR.scala:2278)
# 	at is.hail.expr.ir.MatrixMapRows.execute(MatrixIR.scala:1145)
# 	at is.hail.expr.ir.CastMatrixToTable.execute(TableIR.scala:1293)
# 	at is.hail.expr.ir.Interpret$$anonfun$apply$1.apply$mcJ$sp(Interpret.scala:710)
# 	at is.hail.expr.ir.Interpret$$anonfun$apply$1.apply(Interpret.scala:710)
# 	at is.hail.expr.ir.Interpret$$anonfun$apply$1.apply(Interpret.scala:710)
# 	at scala.Option.getOrElse(Option.scala:121)
# 	at is.hail.expr.ir.Interpret$.apply(Interpret.scala:710)
# 	at is.hail.expr.ir.Interpret$.apply(Interpret.scala:92)
# 	at is.hail.expr.ir.Interpret$.apply(Interpret.scala:62)
# 	at is.hail.expr.ir.Interpret$.interpretJSON(Interpret.scala:21)
# 	at is.hail.expr.ir.Interpret.interpretJSON(Interpret.scala)
# 	at sun.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
# 	at sun.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:62)
# 	at sun.reflect.DelegatingMethodAccessorImpl.invoke(DelegatingMethodAccessorImpl.java:43)
# 	at java.lang.reflect.Method.invoke(Method.java:498)
# 	at py4j.reflection.MethodInvoker.invoke(MethodInvoker.java:244)
# 	at py4j.reflection.ReflectionEngine.invoke(ReflectionEngine.java:357)
# 	at py4j.Gateway.invoke(Gateway.java:280)
# 	at py4j.commands.AbstractCommand.invokeMethod(AbstractCommand.java:132)
# 	at py4j.commands.CallCommand.execute(CallCommand.java:79)
# 	at py4j.GatewayConnection.run(GatewayConnection.java:214)
# 	at java.lang.Thread.run(Thread.java:748)
#
# java.io.FileNotFoundException: File file:/tmp/interval_file.tsv does not exist
# 	at org.apache.hadoop.fs.RawLocalFileSystem.deprecatedGetFileStatus(RawLocalFileSystem.java:635)
# 	at org.apache.hadoop.fs.RawLocalFileSystem.getFileLinkStatusInternal(RawLocalFileSystem.java:861)
# 	at org.apache.hadoop.fs.RawLocalFileSystem.getFileStatus(RawLocalFileSystem.java:625)
# 	at org.apache.hadoop.fs.FilterFileSystem.getFileStatus(FilterFileSystem.java:442)
# 	at org.apache.hadoop.fs.ChecksumFileSystem$ChecksumFSInputChecker.<init>(ChecksumFileSystem.java:146)
# 	at org.apache.hadoop.fs.ChecksumFileSystem.open(ChecksumFileSystem.java:347)
# 	at org.apache.hadoop.fs.FileSystem.open(FileSystem.java:787)
# 	at org.apache.hadoop.mapred.LineRecordReader.<init>(LineRecordReader.java:109)
# 	at org.apache.hadoop.mapred.TextInputFormat.getRecordReader(TextInputFormat.java:67)
# 	at org.apache.spark.rdd.HadoopRDD$$anon$1.liftedTree1$1(HadoopRDD.scala:251)
# 	at org.apache.spark.rdd.HadoopRDD$$anon$1.<init>(HadoopRDD.scala:250)
# 	at org.apache.spark.rdd.HadoopRDD.compute(HadoopRDD.scala:208)
# 	at org.apache.spark.rdd.HadoopRDD.compute(HadoopRDD.scala:94)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.rdd.MapPartitionsRDD.compute(MapPartitionsRDD.scala:38)
# 	at org.apache.spark.rdd.RDD.computeOrReadCheckpoint(RDD.scala:323)
# 	at org.apache.spark.rdd.RDD.iterator(RDD.scala:287)
# 	at org.apache.spark.scheduler.ResultTask.runTask(ResultTask.scala:87)
# 	at org.apache.spark.scheduler.Task.run(Task.scala:108)
# 	at org.apache.spark.executor.Executor$TaskRunner.run(Executor.scala:338)
# 	at java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1149)
# 	at java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:624)
# 	at java.lang.Thread.run(Thread.java:748)
#
#
#
#
# Hail version: 0.2.6-07b91f4bd378
# Error summary: FileNotFoundException: File file:/tmp/interval_file.tsv does not exist
# ERROR: (gcloud.dataproc.jobs.submit.pyspark) Job [bf7daf429b9b47b48852b57b9d1df7cf] entered state [ERROR] while waiting for [DONE].
