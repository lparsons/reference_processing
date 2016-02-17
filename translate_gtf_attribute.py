#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import re
import sys


def main():
    """Translate gtf attribute in based on a lookup table
    :returns: TODO

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    parser.add_argument('id_file')
    parser.add_argument('--gtf_attribute', default="gene_id",
                        help="gtf attribute to translate "
                        "(default %(default)s)")
    parser.add_argument('-f', '--from-field',  type=int, default=1,
                        help="Column of id file to translate from "
                        "(default %(default)s)")
    parser.add_argument('-t', '--to-field',  type=int, default=2,
                        help="Column of id file to translate to "
                        "(default %(default)s)")
    args = parser.parse_args()
    translate_gtf_attribute(args.gtf_file, args.id_file, sys.stdout,
                            args.gtf_attribute,
                            args.from_field,
                            args.to_filed)


def translate_gtf_attribute(gtf_file, id_file, output_file,
                            gtf_attribute, from_field, to_field):
    """TODO: Docstring for translate_gtf_attribute.

    :gtf_file: TODO
    :id_file: TODO
    :output_file: TODO
    :gtf_attribute: TODO
    :from_field: TODO
    :to_field: TODO
    :returns: TODO

    """
    pass
    translations = dict()
    for line in open(id_file, 'rU'):
        fields = line.strip().split("\t")
        translations.setdefault(fields[from_field-1],
                                fields[to_field-1])

    regex = re.compile(r'(.*%s ")([^"]+)(".*)' % gtf_attribute)
    for line in open(gtf_file, 'rU'):
        m = re.match(regex, line)
        if m:
            output_file.write("%s%s%s\n" % (m.group(1),
                              translations[m.group(2)],
                              m.group(3)))
        else:
            output_file.write("%s" % line)

# When called as script from snakemake
if 'snakemake' in globals():
    fh = open(snakemake.output[0], 'w')
    translate_gtf_attribute(snakemake.input['gtf'],
                     snakemake.input['lookup'],
                     fh, "gene_id", 1, 3)
    fh.close()

if __name__ == "__main__":
    main()
