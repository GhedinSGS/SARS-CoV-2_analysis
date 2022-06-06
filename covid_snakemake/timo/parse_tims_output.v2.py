# Written by: Mohammed Khalfan

# Downloaded from git: https://github.com/gencorefacility/MAD/blob/master/bin/parse_tims_output.py

# to run:
import sys
from pyfaidx import Fasta
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r',required=True,help='Indicate path and reference fasta file')
parser.add_argument('--var','-v',required=True,help='Indicate path to variant files')
parser.add_argument('--ignore_binom','-b',default = False, action="store_true",help='Indicate whether to use binomial test column')
parser.add_argument('--output','-o',required=True,help='Indicate output vcf')
args = parser.parse_args()


fa = args.ref # input fasta name
path = args.var
out = args.output
ignore_binom = args.ignore_binom

contigs = Fasta(fa)

out=args.output

with open(out, 'w') as vcf:
	vcf.write("##fileformat=VCFv4.2\n")
	vcf.write("##source=timo\n")
	vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">\n')
	vcf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	vcf.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">\n')
	vcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	vcf.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">\n')
	vcf.write("##reference=file://{}\n".format(fa))

	vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
	for chrom in contigs.keys():
		csv = glob.glob(path)[0]
		ref_seq = contigs[chrom]
		with open(csv, 'r') as infile:
			# read/skip the first line
			# saving it but don't really need it
			header_line = infile.readline()
			for line in infile:
				do_write = False
				columns = line.split(',')
				pos = int(columns[2])
				major = columns[3]
				majorfreq = float(columns[4])
				minor = columns[5]
				minorfreq = float(columns[6]) if columns[6] is not '' else ''
				totalcount = int(columns[13])
				# Binom must be True to accept a minor allele
				binom = columns[7] == "True"
				if ignore_binom: binom = True

				ref = ref_seq[pos - 1]

				if major == ref and minor is '': 
					continue
				elif major == "N":
					continue
				elif major == ref and minor != '' and binom:
					if minor != "-":
						alt = minor
						af = minorfreq
						dp = totalcount
						ad = "{},{}".format(round(dp*(1-af)), round(dp*af))
						do_write = True
				elif major != ref and (minor == '' or minor == ref):
					if major != "-":
						alt = major
						af = majorfreq
						dp = totalcount
						ad = "{},{}".format(round(dp*(1-af)), round(dp*af))
						do_write = True
				 # elif major != ref and minor != '' and minor != ref and binom:
					# alt = "{},{}".format(major, minor)
					# af = "{},{}".format(majorfreq, minorfreq)
					# dp = totalcount
					# ad = "{},{},{}".format(round(dp*(1 - majorfreq - minorfreq)), round(majorfreq*dp), round(minorfreq*dp))
					# do_write = True

				if do_write:
					info = "AF={};DP={}".format(af,dp)
					f = "GT:AD:DP:GQ:PL	1:{}:{}:1:1".format(ad, dp)
					vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
						.format(chrom, pos, '.', ref, alt, '.', '.', info, f))
