#!/usr/bin/env python
# encoding: utf-8
""" Extract selected GFF attributes as columns

Flatten children and roll attributes up by default

"""

import tempfile
import argparse
import gffutils
import sys

__author__ = "Lance Parsons <lparsons@princeton.edu>"
__license__ = "BSD 2-Clause"

dbxref = 'Dbxref'


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_file')
    parser.add_argument('--in_memory', action='store_true')
    parser.add_argument('--temp_file', action='store_true')
    parser.add_argument('--do-not-flatten', action='store_true',
                        help='Do not flatten child attributes')
    parser.add_argument('-f', '--feature_type',
                        default='all',
                        help='Feature type to parse, default: %(default)s')
    parser.add_argument('--attributes', nargs='+',
                        default=['ID', 'Name', 'Product'],
                        help="Attributes to extract, default: %(default)s")
    parser.add_argument('--dbxref_attributes', nargs='+',
                        default=['GeneID'])
    args = parser.parse_args()

    print_attributes(args.gff_file, sys.stdout, args.feature_type,
                     args.attributes, args.dbxref_attributes,
                     args.in_memory, args.temp_file,
                     not(args.do_not_flatten))


def print_attributes(gff_file, output_file, feature_type, attributes,
                     dbxref_attributes, in_memory=False, temp_file=False,
                     flatten=True):
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
    gff_fields = ("seqid", "source", "type", "start",
                  "end", "score", "strand", "phase")
    all_attributes = list()
    for a in attributes:
        all_attributes.append(a)
    for dbxa in dbxref_attributes:
        all_attributes.append(dbxa)
    output_file.write("%s\t%s\n" % ("\t".join(gff_fields),
                                    ("\t".join(all_attributes))))

    features = None
    # Flatten to top level features
    if flatten:
        if feature_type == "all":
            features = top_level_features(db)
        else:
            features = db.features_of_type(feature_type)
        for f in features:
            values = flatten_attributes(db, dict(), f,
                                        attributes, dbxref_attributes)
            write_attributes_as_tsv(output_file, f, values, all_attributes)

    # Else don't flatten, just print features
    else:
        if feature_type == "all":
            features = db.all_features()
        else:
            features = db.features_of_type(feature_type)
        for f in features:
            values = get_attributes_as_dict(f, attributes, dbxref_attributes)
            write_attributes_as_tsv(output_file, f, values, all_attributes)


def write_attributes_as_tsv(output_file, feature, values, attribute_list):
    """Outuput attributes values in order of attribute_list

    :output_file: TODO
    :feature: TODO
    :values: TODO
    :attribute_list: TODO
    :returns: TODO

    """
    values_list = list()
    for a in attribute_list:
        v = values.get(a)
        if v is None:
            v = ""
        values_list.append(v)
    output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                      (feature.seqid, feature.source, feature.featuretype,
                       feature.start, feature.end, feature.score,
                       feature.strand, feature.frame,
                       "\t".join(values_list)))


def flatten_attributes(db, feature_attrs, feature, attributes,
                       dbxref_attributes):
    """TODO: Docstring for flatten_attributes.

    :db: TODO
    :feature_attrs: TODO
    :feature: TODO
    :attributes: TODO
    :dbxref_attributes: TODO
    :returns: TODO

    """
    result = merge_two_dicts(
            get_attributes_as_dict(feature, attributes, dbxref_attributes),
            feature_attrs)
    for child in db.children(feature.id):
        result = merge_two_dicts(
            get_attributes_as_dict(child, attributes, dbxref_attributes),
            result)
    return result


def get_attributes_as_dict(feature, attributes, dbxref_attributes):
    """TODO: Docstring for get_attributes_as_dict.

    :feature: TODO
    :attributes: TODO
    :dbxref_attributes: TODO
    :returns: TODO

    """
    feature_attrs = dict()
    for attr in attributes:
        if attr in feature.attributes:
            feature_attrs[attr] = ",".join(feature.attributes.get(attr))
    if 'Dbxref' in feature.attributes:
        feature_dbxref_attributes = \
            dict(item.split(":")
                 for item in feature.attributes['Dbxref'])
        for attr in dbxref_attributes:
            if attr in feature_dbxref_attributes:
                feature_attrs[attr] = feature_dbxref_attributes.get(attr)
    return feature_attrs


def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z


def top_level_features(db):
    """Generator that returns top level features from db

    :db: gffutils.FeatureDB
    :returns: Iterator of top-level features (having no parent)

    """
    # Get all top level features
    top_level_features = db.execute(
            "select * "
            "from features "
            "where id not in (select distinct(child) from relations)")
    for row in top_level_features:
        yield(gffutils.Feature(**row))


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
