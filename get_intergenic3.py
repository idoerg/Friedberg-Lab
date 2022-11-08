#!/usr/bin/env python
import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
# Copyright(C) 2009-2019 Iddo Friedberg
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment
# This programs reads a GenBank file, find intergenic regions that are
# longer than intergene_length,
# and writes out a multi-fasta file with the intergenic regions.
def get_interregions(genbank_path,intergene_length=100,stats=False):
    seq_record = next(SeqIO.parse(open(genbank_path), "genbank"))
    
    cds_list = []
    intergenic_records = []
    plus_lens = []
    minus_lens = []
    # Loop over the genome file, get the 
    for fnum, feature in enumerate(seq_record.features):
        if feature.type == 'gene':  
            mystart = feature.location._start.position
            myend = feature.location._end.position
            cds_list.append((mystart,myend,feature.strand))
    for i,pospair in enumerate(cds_list[1:]):
        # Compare current start position to previous end position
        last_end = cds_list[i][1]
        this_start = pospair[0]
        strand = pospair[2]
        if this_start - last_end >= intergene_length:
            intergene_seq = seq_record.seq[last_end:this_start]
            if strand == -1:
                intergene_seq = intergene_seq.reverse_complement()
                strand_string = "-"
                minus_lens.append(abs(this_start-last_end)+1)
            else:
                strand_string = "+"
                plus_lens.append(abs(this_start-last_end)+1)
            intergenic_records.append( 
                  SeqRecord(intergene_seq,id="%s-ign-%d" % (seq_record.name,i),
                  description="%s %d %d %s" % (seq_record.name, last_end+1,
                                                        this_start,strand_string)))
    outpath = os.path.splitext(os.path.basename(genbank_path))[0] + ".ign"
    SeqIO.write(intergenic_records, open(outpath,"w"), "fasta")
    if stats:
        intergenic_stats(intergenic_records,genbank_path,seq_record.name)

def intergenic_stats(intergenic_records,genbank_path,seq_name):
    ign_lens = []
    raw_stats_out = open("%s.%s" % (os.path.splitext(os.path.basename(genbank_path))[0],
                         "rst"),"w")
    sum_stats_out = open("%s.%s" % (os.path.splitext(os.path.basename(genbank_path))[0],
                         "sst"),"w")
    for seq_record in intergenic_records:
        seqid = seq_record.id 
        end, start, strand = seq_record.description.split()[1:]
        ign_len = int(start) - int(end) + 1
        raw_stats_out.write("%s\t%d\t%s\n" % (seqid, ign_len, strand))
        ign_lens.append(ign_len)
    lens_mean = np.mean(ign_lens) 
    lens_std = np.std(ign_lens) 
    sum_stats_out.write("%s\t%.4f\t%.4f\n" % (seq_name, lens_mean, lens_std))
    raw_stats_out.close()
    sum_stats_out.close()
        
def multi_interregions(infiles, minlen, stats):
    for inpath in infiles:
        print("Processing:", inpath)
        get_interregions(inpath, minlen, stats)

def inpath_interregions(pathfile, minlen, stats):
    infiles = []
    for inline in open(pathfile[0]).readlines():
        infiles.append(inline.strip())
    multi_interregions(infiles, minlen, stats)

def violinplots(infiles):
    ign_df = pd.read_table('Pseudomonas_stutzeri_CCUG_29243_177.rst',
                           header=None,names=['seqid','seqlen','strand'])
    

if __name__ == '__main__':
    lineparse = argparse.ArgumentParser()
    lineparse.add_argument('infiles',nargs="+")
    lineparse.add_argument('-s','--stats', action="store_true",  
                           help='write stats of intergenic lengths')
    lineparse.add_argument('-m','--minlen', type=int, 
                            help='minimum intergenic length',default=0)
    lineparse.add_argument('-p','--pathnames', action="store_true",
                           help='file contains pathnames of genbank files')
    args = lineparse.parse_args()
    if args.pathnames:
        inpath_interregions(args.infiles, args.minlen, args.stats)
    else:
        multi_interregions(args.infiles, args.minlen, args.stats)

