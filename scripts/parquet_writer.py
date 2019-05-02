#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from numpy import nan
from collections import OrderedDict


def main():

    print('Pandas version: ', pd.__version__)
    print('Pyarrow version: ', pa.__version__)

    # Create df with nan columns
    data = [ ['a' , 1, nan],
             ['b' , nan, 0.2],
             [nan , 3, 0.3] ]
    df = pd.DataFrame(data, columns=['str', 'int', 'float'])

    # Set data types
    dtypes = { 'str':'object',
               'int':'Int64',
               'bool':'bool',
               'float':'float64' }
    df = df.astype(dtype=dtypes)

    # Write file
    path = 'test.parquet'
    write_parquet(
        df,
        path,
        compression='snappy',
        flavor='spark'
    )

    return 0

def write_parquet(df, path, preserve_index=False, str_list_cols=None, **kwargs):
    ''' Writes a pandas df to parquet, preserving Int64 datatype.
    Params:
        df (pd.df): pandas dataframe
        path (file): path to write parquet to
        preserve_index (bool): whether to keep the pd.df index
        str_list_cols (list of str): column names that contain list of strings
    '''

    # Get pandas.df schema
    pd_dtypes = OrderedDict(df.dtypes)

    # Make pa.Table. Convert pd.df ints to floats in pandas, but specify them
    # as ints in the pyarrow schema
    table = pa.Table.from_pandas(
        df=df.astype(dtype=pd_dtype_int_to_float(pd_dtypes)),
        preserve_index=preserve_index,
        schema=pd_dtype_to_pa_schema(pd_dtypes, str_list_cols=str_list_cols)
    )

    # Write output
    pq.write_table(
        table,
        path,
        **kwargs
    )

    return 0

def pd_dtype_to_pa_schema(dtypes, str_list_cols=None):
    ''' Converts a pandas dtype to a pyarrow schema
    '''
    # Create dict to map pandas to pa type
    # https://arrow.apache.org/docs/python/api/datatypes.html
    #
    type_map = {
        'str': pa.string(),
        'object': pa.string(),
        'Int64': pa.int64(),
        'int64': pa.int64(),
        'bool': pa.bool_(),
        'float64': pa.float64()
    }

    # Construct pyarrow schema
    sc_list = []
    for key, val in dtypes.items():
        if not str_list_cols or not key in str_list_cols:
            sc_list.append( (key, type_map[str(val)]) )
        else:
            sc_list.append( (key, pa.list_(pa.string())) )
    pa_sc = pa.schema(sc_list)

    return pa_sc

def pd_dtype_int_to_float(dtypes):
    ''' Return a copy of the pandas dtype with Int64 type switched to float64
    '''
    dtypes_new = dtypes.copy()
    for key, val in dtypes_new.items():
        if str(val) == 'Int64':
            dtypes_new[key] = 'float64'
    return dtypes_new

if __name__ == '__main__':

    main()
