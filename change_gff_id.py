#!/usr/bin/env python
# encoding: utf-8

import argparse
import sys
import gffutils
import tempfile

# Global ID Map
id_map = dict()


def main():
    """Main
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('gff_file')
    parser.add_argument('--in_memory', action='store_true')
    parser.add_argument('--id_attribute',
                        help="Attribute to use for new id")
    args = parser.parse_args()

    build_id_map(args.gff_file, args.in_memory, args.id_attribute)

    # Change ID
    if args.in_memory:
        filename = ':memory:'
    else:
        tmp = tempfile.NamedTemporaryFile()
        filename = tmp.name
    db = gffutils.create_db(args.gff_file,
                            filename,
                            transform=map_id,
                            merge_strategy="create_unique",
                            keep_order=True)

    for f in db.all_features():
        print f


def map_id(f):
    """Change id using id_map

    :f: TODO
    :returns: TODO

    """
    f.attributes['ID'][0] = id_map[f.attributes['ID'][0]]
    if 'Parent' in f.attributes:
        f.attributes['Parent'][0] = id_map[f.attributes['Parent'][0]]
    return f


def build_id_map(gff_file, in_memory, id_attribute):
    """Get the map from original ids to new id_attribute

    :gff_file: gff filename to parse
    :filename: database filename
    :id_attribute: attribute to use a new id
    :returns: id map dictionary

    """
    pass
    # Get ID map
    if in_memory:
        filename = ':memory:'
    else:
        tmp = tempfile.NamedTemporaryFile()
        filename = tmp.name
    db = gffutils.create_db(gff_file,
                            filename,
                            merge_strategy="create_unique")
    all_ids = dict()
    for f in db.all_features():
        new_id = f.id
        if id_attribute in f.attributes:
            if f.attributes[id_attribute][0] in all_ids:
                sys.stderr.write("Duplicate %s found: %s, using %s" %
                                 id_attribute, f.attributes[id_attribute][0],
                                 f.id)
                new_id = f.id
            else:
                new_id = f.attributes[id_attribute][0]
        all_ids[new_id] = 1
        id_map[f.id] = new_id


if __name__ == '__main__':
    main()
