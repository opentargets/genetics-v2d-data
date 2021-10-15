
# coding: utf-8

from collections import OrderedDict

import argparse
import pandas as pd


def main():

    # Parse args
    args = parse_args()

    # read manifest + EFO mappings

    FINNGEN_EFO=pd.read_csv(args.in_EFO, sep='\t',
                               header=0)

    FINNGEN_manifest = (pd.read_json(args.in_manifest, lines=True).rename(
            columns={
                'phenocode': 'study_id',
                'phenostring': 'trait',
                'category': 'trait_category',
                'num_cases': 'n_cases',
                'num_controls': 'n_controls'
            }
        )
    )

    keep_columns = [
        'study_id',
        'trait',
        'trait_category',
        'n_cases',
        'n_controls'
    ]
    FINNGEN_manifest = FINNGEN_manifest[keep_columns]

    FINNGEN_manifest['study_id'] = "FINNGEN_R5_" + FINNGEN_manifest['study_id']
    FINNGEN_manifest['n_total'] = FINNGEN_manifest['n_cases'] + FINNGEN_manifest['n_controls']


    # Keep valid EFOs

    FINNGEN_valid_EFO=FINNGEN_EFO.loc[FINNGEN_EFO['valid'] == True]
    #FINNGEN_valid_EFO=FINNGEN_EFO
    #FINNGEN_valid_EFO.head()

    FINNGEN_valid_EFO['FINNGEN_ID']="FINNGEN_R5_"+FINNGEN_valid_EFO['NAME']
    to_keep= OrderedDict([
            ('FINNGEN_ID', 'study_id'),
            ('LONGNAME', 'trait_reported'),
            ('efo_cls', 'trait EFOs')])

    FINNGEN_valid_EFO_trimmed=FINNGEN_valid_EFO.loc[:, to_keep.keys()].rename(columns=to_keep)

    # Map manifest to efo trait and ID:


    FINNGEN_manifest['trait_efos']=FINNGEN_manifest.apply(lambda x: Get_EFO_mapping_list(x['study_id'], FINNGEN_valid_EFO_trimmed), axis=1)


    FINNGEN_manifest['trait_reported']=FINNGEN_manifest.apply(lambda x: Get_EFO_trait_list(x['study_id'], FINNGEN_valid_EFO_trimmed), axis=1)
    FINNGEN_manifest.head()

    # Format table:

    FINNGEN_manifest['pmid']=""
    FINNGEN_manifest['pub_date']='2021-5-11'
    FINNGEN_manifest['pub_author']='FINNGEN_R5'
    FINNGEN_manifest['ancestry_initial']="European="+FINNGEN_manifest['n_total'].astype(str)
    FINNGEN_manifest['n_replication']=0
    FINNGEN_manifest['ancestry_replication']=""
    FINNGEN_manifest['pub_journal']=""
    FINNGEN_manifest['pub_title']=""

    cols = OrderedDict([
        ('study_id', 'study_id'),
        ('pmid', 'pmid'),
        ('pub_date', 'pub_date'),
        ('pub_journal', 'pub_journal'),
        ('pub_title', 'pub_title'),
        ('pub_author', 'pub_author'),
        ('trait', 'trait_reported'),
        ('trait_efos', 'trait_efos'),
        ('ancestry_initial', 'ancestry_initial'),
        ('ancestry_replication', 'ancestry_replication'),
        ('n_total', 'n_initial'),
        ('n_cases', 'n_cases'),
        ('n_replication', 'n_replication')
        ])

    FINNGEN_manifest=FINNGEN_manifest.loc[:, list(cols.keys())].rename(columns=cols)

    Additional_row=pd.Series(["GCST90013791", "",
                              "2021-02-22", "", "", "Crouch D", "Type 1 diabetes", 
                              ["EFO_0001359"], "European=7977", "", "7977", "3983", "0"   ])

    row_df=pd.DataFrame([Additional_row])
    row_df=row_df.set_axis(list(FINNGEN_manifest.columns.values), axis=1, inplace=False)

    FINNGEN_manifest_new=pd.concat([FINNGEN_manifest, row_df], ignore_index=True)
    FINNGEN_manifest_new['trait_reported'][FINNGEN_manifest_new['study_id']=="FINNGEN_R5_I9_HEARTFAIL_AND_CHD"]="cardiovascular disease"

    # Load no finngen study table
    gwas = pd.read_json(args.in_study_table, 
                        orient='records', lines=True)
    # Merge
    merged = pd.concat([gwas, FINNGEN_manifest], sort=False)
    print(merged)
    # Write
    merged.to_json(args.outf, orient='records', lines=True)

def Get_EFO_mapping_list(study_ids, EFOs):
    EFO_mappings=[]
    Found_EFOs=EFOs.loc[EFOs['study_id'].str.contains(study_ids)]['trait EFOs'].unique()
    return(Found_EFOs)

def Get_EFO_trait_list(study_ids, EFOs):
    EFO_mappings=[]
    Found_EFOs=EFOs.loc[EFOs['study_id'].str.contains(study_ids)]['trait_reported'].unique()
    return(Found_EFOs)


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_manifest', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_EFO', metavar="<str>", type=str, required=True)
    parser.add_argument('--in_study_table', metavar="<str>", type=str, required=True)
    parser.add_argument('--outf', metavar="<str>", help=("Output"), type=str, required=True)
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()

#    FINNGEN_manifest_path="/home/xg1/r5_finngen.json"
#    FINNGEN_EFO_path="/home/xg1/finngen_df5_efo_mapping.lastest.nonull.tsv"
# "/home/xg1/genetics-v2d-data/tmp/210526/merged_study_table.old.json"
# "/home/xg1/genetics-v2d-data/tmp/210526/merged_study_table.json"

#python Manual_FINNGENs --in_manifest r5_finngen.json --in_EFO finngen_df5_efo_mapping.lastest.nonull.tsv 
#    --in_study_table tmp/version_date/merged_study_table.old.json --outf tmp/version_date/merged_study_table.json

