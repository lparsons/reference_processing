#!/usr/bin/env python3
'''
gtf2bed.py converts GTF file to BED file.
Usage: gtf2bed.py {OPTIONS} [.GTF file]

History
    http://allaboutbioinfo.blogspot.se/2011/08/converting-cufflinks-gtf-predictions-to.html
    Nov.5th 2012:
        1. Allow conversion from general GTF files (instead of only Cufflinks
           supports).
        2. If multiple identical transcript_id exist, transcript_id will be
           appended a string like "_DUP#" to separate.
    2014-11-14 Lance Parsons
        Cleaned up code style
        Added option to specify feature types recognized
'''

import argparse
import sys
import re

allids = {}


def main():
    parser = argparse.ArgumentParser(epilog='''
    Note:
    1. Only "exon" and "transcript" are recognized in the feature
    field (3rd field).\n
    2. In the attribute list of .GTF file, the script tries to find
    "gene_id", "transcript_id" and "FPKM" attribute, and convert
    them as name and score field in .BED file.\n
    Author: Wei Li (li.david.wei AT gmail.com)
    ''')
    parser.add_argument('gtf_file')
    parser.add_argument('-c', '--color', default='255,0,0', help='Specify the '
                        'color of the track. This is a RGB value represented '
                        'as "r,g,b". Default 255,0,0 (red)')
    parser.add_argument('--exon_features', default=['exon'],
                        nargs='+', help='Specify the gtf '
                        'features that are recognized as exons')
    args = parser.parse_args()

    estart = []
    eend = []
    # read lines one to one
    nline = 0
    prevfield = []
    prevtransid = ''
    for lines in open(args.gtf_file):
        field = lines.strip().split('\t')
        nline = nline+1
        if len(field) < 9:
            print('Error: the GTF should has at least 9 fields at line ' +
                  str(nline), file=sys.stderr)
            continue
        if field[1] != 'Cufflinks':
            pass
            # print('Warning: the second field is expected to be '
            #       '\'Cufflinks\' at line '+str(nline), file=sys.stderr)
        if field[2] not in ['transcript'] + args.exon_features:
            # print(field[2])
            # print('Error: the third field is expected to be one of: '
            #       '%s at line %i, saw %s' %
            #       (['transcript'] + args.exon_features, nline, field[2]),
            #       file=sys.stderr)
            continue
        transid = re.findall(r'transcript_id \"([\w\.]+)\"', field[8])
        if len(transid) > 0:
            transid = transid[0]
        else:
            transid = ''
        if field[2] == 'transcript' or (prevtransid != '' and transid != '' and
                                        transid != prevtransid):
            # print('prev:'+prevtransid+', current:'+transid)
            # print(estart, eend, prevfield)
            # A new transcript record, write
            if len(estart) != 0:
                printbedline(estart, eend, prevfield, nline, args.color)
            estart = []
            eend = []
        prevfield = field
        prevtransid = transid
        if field[2] in ['exon', 'CDS']:
            try:
                est = int(field[3])
                eed = int(field[4])
                estart += [est]
                eend += [eed]
            except ValueError:
                print('Error: non-number fields at line ' + str(nline),
                      file=sys.stderr)
    # the last record
    if len(estart) != 0:
        printbedline(estart, eend, field, nline, args.color)


def printbedline(estart, eend, field, nline, color):
    try:
        estp = estart[0]-1
        eedp = eend[-1]
        # use regular expression to get transcript_id, gene_id and expression
        # level
        geneid = re.findall(r'gene_id \"([\w\.]+)\"', field[8])
        transid = re.findall(r'transcript_id \"([\w\.]+)\"', field[8])
        fpkmval = re.findall(r'FPKM \"([\d\.]+)\"', field[8])
        if len(geneid) == 0:
            print('Warning: no gene_id field', file=sys.stderr)
        else:
            geneid = geneid[0]
        if len(transid) == 0:
            print('Warning: no transcript_id field', file=sys.stderr)
            transid = 'Trans_'+str(nline)
        else:
            transid = transid[0]
        if transid in allids.keys():
            transid2 = transid+'_DUP'+str(allids[transid])
            allids[transid] = allids[transid]+1
            transid = transid2
        else:
            allids[transid] = 1
        if len(fpkmval) == 0:
            # print('Warning: no FPKM field',file=sys.stderr)
            fpkmval = '100'
        else:
            fpkmval = fpkmval[0]
        fpkmint = round(float(fpkmval))
        print(field[0] + '\t' + str(estp) + '\t' + str(eedp) + '\t' +
              transid + '\t' + str(fpkmint) + '\t' + field[6] + '\t' +
              str(estp) + '\t' + str(eedp) + '\t' + color + '\t' +
              str(len(estart))+'\t', end='')
        seglen = [eend[i]-estart[i]+1 for i in range(len(estart))]
        segstart = [estart[i]-estart[0] for i in range(len(estart))]
        strl = str(seglen[0])
        for i in range(1, len(seglen)):
            strl += ',' + str(seglen[i])
        strs = str(segstart[0])
        for i in range(1, len(segstart)):
            strs += ',' + str(segstart[i])
        print(strl + '\t' + strs)
    except ValueError:
        print('Error: non-number fields at line '+str(nline),
              file=sys.stderr)


if __name__ == '__main__':
    main()
