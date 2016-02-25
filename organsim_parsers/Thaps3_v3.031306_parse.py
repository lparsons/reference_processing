#!/usr/bin/env python
# encoding: utf-8

'''
Parse JGI GFF file for Thaps3_v3.031306 and output GTF formatted file
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
    Convert attributes to GTF standard attributes

    :f: gffutils feature to transform
    :returns: transformed gffutils feature

    """
    if f.featuretype == 'CDS':
        f.attributes['gene_id'] = f.attributes['proteinId']
        f.attributes.pop('proteinId', None)
        if 'transcriptId' in f.attributes:
            f.attributes['transcript_id'] = f.attributes['transcriptId']
            f.attributes.pop('transcriptId', None)
        else:
            f.attributes['transcript_id'] = f.attributes['gene_id']

    if f.featuretype == 'exon':
        f.attributes['transcript_id'] = f.attributes['transcriptId']
        f.attributes.pop('transcriptId', None)
        if 'proteinId' in f.attributes:
            f.attributes['gene_id'] = f.attributes['proteinId']
            f.attributes.pop('proteinId', None)
        else:
            f.attributes['gene_id'] = f.attributes['transcript_id']
    #     f.attributes['ID'][0] += "_" + f.attributes['rank'][0]
    #     f.id = f.attributes['ID'][0]
    # if 'PARENT' in f.attributes:
    #     f.attributes['Parent'] = f.attributes['PARENT']
    #     f.attributes.pop('PARENT', None)
    return f


def update_start_stop_codons(db):
    """Update gene_id and transcript_id for start and stop codon
    Find matching CDS and transfer the gene_ids and transcript_ids
    Additionally, substract the stop_codon from the CDS region per GTF spec
    """
    # Update start codon transcript_id and gene_id
    start_codon_list = list()
    for start_codon in db.features_of_type("start_codon"):
        start_codon_list.append(start_codon)
    for start_codon in start_codon_list:
        # sys.stderr.write("%s\n" % sc)
        for f in db.region(start_codon):
            if (f.featuretype == 'CDS' and
                    start_codon.attributes['name'] == f.attributes['name']):
                start_codon.attributes['gene_id'] = f.attributes['gene_id']
                start_codon.attributes['transcript_id'] = \
                    f.attributes['transcript_id']
                # sys.stderr.write("%s\n" % f)
        # sys.stderr.write("%s\n" % start_codon)
    db.update(start_codon_list, make_backup=False,
              merge_strategy="replace")

    # Update stop codon transcript_id and gene_id
    # Remove stop codon coordinates from CDS region (per GTF spec)
    sc_list = list()
    for sc in db.features_of_type("stop_codon"):
        sc_list.append(sc)
    for sc in sc_list:
        # sys.stderr.write("1 %s: %s\n" % (sc.id, sc))
        for f in db.region(sc):
            if (f.featuretype == 'CDS' and
                    sc.attributes['name'] == f.attributes['name']):
                # sys.stderr.write("2 %s: %s\n" % (f.id, f))
                sc.attributes['gene_id'] = f.attributes['gene_id']
                sc.attributes['transcript_id'] = f.attributes['transcript_id']
                if f.strand == "-":
                    f.start = sc.end + 1
                else:
                    f.end = sc.start - 1
                # sys.stderr.write("3 %s: %s\n" % (f.id, f))
                db.update([f], make_backup=False,
                          merge_strategy="replace")
        # sys.stderr.write("4 %s: %s\n" % (sc.id, sc))
    db.update(sc_list, make_backup=False,
              merge_strategy="replace")


def main():
    """Main
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('gff_file')
    parser.add_argument('--in_memory', action='store_true')
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
                            # id_spec=['ID'],
                            transform=transform,
                            gtf_transcript_key="transcriptId",
                            gtf_gene_key="proteinId",
                            keep_order=True)
    update_start_stop_codons(db)
    output_original(db)


def output_original(db):
    """Print original file (with parsing modifications)

    :db: gffutils database
    :returns: None

    """
    for f in db.all_features():
        print f


if __name__ == '__main__':
    main()
