import csv
import json
import sys
import gzip

csvfile = gzip.open(sys.argv[1], 'rt')
jsonfile = gzip.open(sys.argv[2], 'wt')

fields = ["study_id_A", "index_variantid_b37_A", "study_id_B", "index_variantid_b37_B", "set_type", "distinct_A", "overlap_AB", "distinct_B"]

fields_array = []

fields_int = ["distinct_A", "overlap_AB", "distinct_B"]
def remove_if_empty(line, fields):
    for field in fields:
        if field in line and len(line[field]) == 0:
            del(line[field])
    return line

def string_to_int(line, fields):
    for field in fields:
        if field in line and len(line[field]) > 0:
            line[field] = int(round(float(line[field])))
    return line

def string_to_list(line, fields, sep=';'):
    for field in fields:
        if field in line and len(line[field]) > 0:
            tokens = line[field].split(sep)
            line[field] = list(map(lambda el: el.strip(), tokens))
        else:
            line[field] = []
    return line

reader = csv.DictReader(csvfile, delimiter="\t")
for row in reader:
    cleaned_line = string_to_int(
            string_to_list(
                remove_if_empty(row,fields+fields_array+fields_int),fields_array), fields_int)
    json.dump(row, jsonfile)
    jsonfile.write('\n')

csvfile.close()
jsonfile.close()
