#!/usr/bin/env python
# encoding: utf-8

import sys
import re
import gffutils


def add_chr(d):
    cleanseq = re.sub('-', '_', d.seqid)
    d.attributes['ID'][0] = "%s_%s" % (cleanseq, d.attributes['ID'][0])
    if 'Parent' in d.attributes:
        d.attributes['Parent'][0] = "%s_%s" % (cleanseq,
                                               d.attributes['Parent'][0])
    return d

db = gffutils.create_db(sys.argv[1], ':memory:', force=True,
                        merge_strategy="merge", transform=add_chr)

for f in db.all_features():
    print f
