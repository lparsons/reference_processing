#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import re
import sys


def main():
    """Add empty gene_id attribute if it's missing to ensure file conforms to
    specifications. The gene_id is sometimes missing from gffread output.

    :returns: None

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    fix_gtf_records(args.gtf_file, sys.stdout)


def fix_gtf_records(gtf_file, output_file):
    """Add missing gene_id attributes to gtf records

    :gtf_file: TODO
    :output_file: TODO

    :returns: None

    """
    regex = re.compile(r'(.*gene_id ")([^"]+)(".*)')
    for line in open(gtf_file, 'rU'):
        m = re.match(regex, line)
        if m:
            output_file.write(line)
        else:
            output_file.write('%s gene_id "";\n' % line.strip())

# When called as script from snakemake
if 'snakemake' in globals():
    fh = open(snakemake.output[0], 'w')
    fix_gtf_records(snakemake.input[0], fh)
    fh.close()

if __name__ == "__main__":
    main()
