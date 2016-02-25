#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import re
import sys


def main():
    """Add gene_id attribute if it's missing to ensure file conforms to
    specifications.
    If gene_id is missing, use transcript_id else use empty string.
    The gene_id is sometimes missing from gffread output.

    :returns: None

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    parser.add_argument('--use_transcript_id_as_gene_id', action="store_true",
                        help='If gene_id is missing, use transcript_id if '
                        'set, otherwise use empty string')
    args = parser.parse_args()
    fix_gtf_records(args.gtf_file, sys.stdout)


def fix_gtf_records(gtf_file, output_file):
    """Add missing gene_id attributes to gtf records

    :gtf_file: TODO
    :output_file: TODO

    :returns: None

    """
    gene_id_regex = re.compile(r'(.*gene_id ")([^"]+)(".*)')
    transcript_id_regex = re.compile(r'(.*transcript_id ")([^"]+)(".*)')
    for line in open(gtf_file, 'rU'):
        gene_id_match = re.match(gene_id_regex, line)
        if gene_id_match:
            output_file.write(line)
        else:
            transcript_id_match = re.match(transcript_id_regex, line)
            if transcript_id_match:
                gene_id = transcript_id_match.group(2)
            else:
                gene_id = ""
            output_file.write('%s gene_id "%s";\n' % (line.strip(), gene_id))

# When called as script from snakemake
if 'snakemake' in globals():
    fh = open(snakemake.output[0], 'w')
    fix_gtf_records(snakemake.input[0], fh)
    fh.close()

if __name__ == "__main__":
    main()
