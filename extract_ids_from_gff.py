#!/usr/bin/env python
# encoding: utf-8

import tempfile
import re
import argparse
import gffutils


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_file')
    parser.add_argument('--in_memory', action='store_true')
    parser.add_argument('--attributes', nargs='+',
                        default=['ID', 'Name', 'Dbxref', 'locus_tag'])
    args = parser.parse_args()

    if args.in_memory:
        filename = ':memory:'
    else:
        tmp = tempfile.NamedTemporaryFile()
        filename = tmp.name
    db = gffutils.create_db(args.gff_file, filename)

    for f in db.all_features():
        values = []
        for a in args.attributes:
            value = ''
            try:
                value = ",".join(f.attributes[a])
            except:
                pass
            values.append(value)
        print "\t".join(values)


if __name__ == '__main__':
    main()
