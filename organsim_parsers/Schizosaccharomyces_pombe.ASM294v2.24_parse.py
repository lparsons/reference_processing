#!/usr/bin/env python
# encoding: utf-8

'''
Parse ENSEMBL GFF file for Schizosaccharomyces pombe
'''

import argparse
import re
import gffutils
import tempfile
# import sys


def add_seqid(d):
    ''' Add seqid to ID and Parent attributesa
    Useful when combining multiple NCBI GFF files
    (e.g. bacteria plus plasmid)
    '''
    cleanseq = re.sub('-', '_', d.seqid)
    d.attributes['ID'][0] = "%s_%s" % (cleanseq, d.attributes['ID'][0])
    if 'Parent' in d.attributes:
        d.attributes['Parent'][0] = "%s_%s" % (cleanseq,
                                               d.attributes['Parent'][0])
    return d


def transform(f):
    """Transform feature
    Add rank to CDS ID to make them unique

    :f: gffutils feature to transform
    :returns: transformed gffutils feature

    """
    if f.featuretype == 'CDS':
        f.attributes['ID'][0] += "_" + f.attributes['rank'][0]
        f.id = f.attributes['ID'][0]
    if 'PARENT' in f.attributes:
        f.attributes['Parent'] = f.attributes['PARENT']
        f.attributes.pop('PARENT', None)
    return f


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_file')
    parser.add_argument('--in_memory', action='store_true')
    parser.add_argument('--output_format', default='original')
    args = parser.parse_args()

    if args.in_memory:
        filename = ':memory:'
    else:
        tmp = tempfile.NamedTemporaryFile()
        filename = tmp.name

    db = gffutils.create_db(args.gff_file,
                            filename,
                            # merge_strategy="create_unique",
                            # merge_strategy="merge",
                            id_spec=['ID'],
                            transform=transform,
                            keep_order=True)

    if args.output_format == 'original':
        output_original(db)
    elif args.output_format == 'simple_named_bed':
        output_simple_named_bed(db)
    elif args.output_format == 'pombase':
        output_pombase(db)


def output_original(db):
    """Print original file (with parsing modifications)

    :db: gffutils database
    :returns: None

    """
    for f in db.all_features():
        print f


def output_simple_named_bed(db):
    """Print named bed file, skipping exons and CDS regions

    :db: gffutils database
    :returns: None

    """
    for f in db.all_features():
        if f.featuretype not in ['exon', 'CDS', 'five_prime_UTR',
                                 'three_prime_UTR', 'transcript',
                                 'repeat_region']:
            feature_name = f.id
            if 'external_name' in f.attributes:
                feature_name = "%s %s" % (feature_name,
                                          f.attributes['external_name'])
            print("%s\t%s\t%s\t%s\t%s\t%s" %
                  (f.seqid, f.start - 1, f.end,
                   feature_name, '0', f.strand))


def output_pombase(db):
    """Print original file with only PomBase features and children

    :db: gffutils database
    :returns: None

    """
    # gffwriter = gffutils.gffwriter.GFFWriter(sys.stdout)
    printed_ids = dict()
    for f in db.all_features():
        if f.id not in printed_ids:
            if f.source == 'PomBase':
                # gffwriter.write_gene_recs(db, f.id)
                print f.id
                print f
                printed_ids[f.id] = 1
                for c in db.children(f.id):
                    print c
                    printed_ids[c.id] = 1


def print_with_children(f):
    """Print feature along with it's children

    :f: gffutils Feature
    :returns: List of ids printed

    """
    pass
if __name__ == '__main__':
    main()
