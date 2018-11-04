#!/usr/bin/env python

import os
import sys
import argparse
import multiprocessing
import dendropy
from collections import defaultdict
from Bio.Seq import translate
from Bio.Seq import reverse_complement

fastafile = open(sys.argv[1], 'r')

blastfile = open(sys.argv[2], 'r')

out = open(sys.argv[3], 'w')

fastadict = dict()

for line in fastafile:
	if '>' in line:
		name = line.strip().split(' ')[0][1:]

		fastadict[name] = next(fastafile).strip()


for line in blastfile:
	info = line.strip().split('\t')
	query = info[0]
	db = info[20]
	start = int(info[8])
	end = int(info[9])

	if end > start:
		seq = fastadict[db]
		seq = seq[start-1:end]
		outname = '>' + query + '\n'
		outseq = seq + '\n'
		out.write(outname)
		out.write(outseq)
	else:
		seq = fastadict[db]
		seq = seq[end-1:start]
		seq = reverse_complement(seq)
		outname = '>' + query + '\n'
		outseq = seq + '\n'
		out.write(outname)
		out.write(outseq)		

















	

