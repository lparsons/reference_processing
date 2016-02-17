#!/usr/bin/env python
# encoding: utf-8

import tempfile
import argparse
import gffutils
import sys

dbxref = 'Dbxref'


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_file')
    parser.add_argument('--in_memory', action='store_true')
    parser.add_argument('--temp_file', action='store_true')
    parser.add_argument('-f', '--feature_type',
                        default='all',
                        help='Feature type to parse, default: %(default)s')
    parser.add_argument('--attributes', nargs='+',
                        default=['ID', 'Name', 'Dbxref'],
                        help="Attributes to extract, default: %(default)s")
    parser.add_argument('--dbxref_attributes', nargs='+',
                        default=['GeneID'])
    args = parser.parse_args()
    print_attributes(args.gff_file, sys.stdout, args.feature_type,
                     args.attributes, args.dbxref_attributes,
                     args.in_memory, args.temp_file)


def print_attributes(gff_file, output_file, feature_type, attributes,
                     dbxref_attributes, in_memory=False, temp_file=False):
    """TODO: Docstring for extract_ids.

    :gff_file: TODO
    :output_file: TODO
    :feature_type: TODO
    :attributes: TODO
    :dbxref_attributes: TODO
    :in_memory: TODO
    :temp_file: TODO
    :returns: TODO

    """
    if in_memory:
        filename = ':memory:'
    elif temp_file:
        tmp = tempfile.NamedTemporaryFile()
        filename = tmp.name
    else:
        filename = "%s.db" % gff_file
    db = None
    try:
        db = gffutils.FeatureDB(filename)
    except:
        db = gffutils.create_db(
            gff_file, filename,
            merge_strategy="create_unique")

    # Print attributes as header
    all_attributes = list()
    for a in attributes:
        if a == dbxref:
            for dbxa in dbxref_attributes:
                all_attributes.append(dbxa)
        else:
            all_attributes.append(a)
    output_file.write("%s\n" % "\t".join(all_attributes))

    # Print features
    features = db.all_features()
    if feature_type is not "all":
        features = db.features_of_type(feature_type)
    for f in features:
        values = parse_feature_attributes(f, attributes,
                                          dbxref_attributes)
        output_file.write("%s\n" % "\t".join(values))


def parse_feature_attributes(feature, attributes, dbxref_attributes):
    values = list()
    for a in attributes:
            # import pdb; pdb.set_trace()
            if a == dbxref:
                # print(feature.attributes[a])
                result = {}
                # import pdb; pdb.set_trace()
                for item in feature.attributes[a]:
                    key, val = item.split(":", 1)
                    if key in result:
                        result[key] = ",".join(result[key], val)
                    else:
                        result[key] = val
                # print(result)
                for dbxa in dbxref_attributes:
                    value = ''
                    if dbxa in result:
                        value = result[dbxa]
                    values.append(value)
            else:
                value = ''
                try:
                    value = ",".join(feature.attributes[a])
                except:
                    pass
                values.append(value)
    return values

# When called as script from snakemake
if 'snakemake' in globals():
    fh = open(snakemake.output[0], 'w')
    print_attributes(snakemake.input[0], fh,
                     snakemake.config["gff_feature_type"],
                     snakemake.config["attributes"],
                     snakemake.config["dbxref_attributes"])
    fh.close()

if __name__ == '__main__':
    main()
