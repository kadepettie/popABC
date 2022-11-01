# -*- coding: utf-8 -*-
"""
Exctracts snps from a vcf file with AA as the only INFO field where
the alt allele is the ancestral allele.

Outputs statistics on numbers and % of SNPs where this is true as well
as SNPs where no ancestral allele information is available.
"""

from os import path
from argparse import ArgumentParser
import pysam
import random
import sys
import numpy as np
import random
from collections import defaultdict
from collections import Counter
import pandas as pd
import pickle as pkl
import gzip as _gzip
import re


def open_zipped(infile, mode='r'):
    """return file handle of file regardless of compressed or not.
    also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # return already open files
    if hasattr(infile, 'write'):
        return infile
    # make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    # refuse to handle non-strings that aren't files.
    if not isinstance(infile, str):
        raise ValueError("i cannot open a filename that isn't a string.")
    # treat '-' appropriately
    if infile == '-':
        if 'w' in mode:
            return sys.stdout
        return sys.stdin
    # if possible open zipped files
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        return _bz2.bz2file(infile, mode)
    # fall back on regular open
    return open(infile, mode)

def parse_vcf_ancestralRef( in_file, out_file, verbose=True, debug=False ):
    """
    Read through lines in vcf file, writing out only those where ancestral allele
    matches the reference allele.
    Return counts of different types of ancestral allele information.
    """

    if verbose:
        print("## Parsing VCF file " + in_file + "\nfor ancestral/derived allele info...")

    vcf_handle = open_zipped(in_file)
    out_handle = open_zipped(out_file,'w')
    header = []
    samples = []

    total_vars = 0
    total_snps = 0
    ref_ancestral = 0
    alt_ancestral = 0
    ambig_ancestral = 0
    other_ancestral = 0
    no_ancestral = 0
    not_biallelic_snp = 0
    lineage_specific_ins = 0
    more_confident_ra = 0
    less_confident_ra = 0
    unknown = 0

    for line in vcf_handle:
        # keep newline character
        # line = line.rstrip()
        ## for now we don't care about the header
        if line.startswith('##'):
            # write header info to file
            # out_handle.write(line)
            continue
        elif line.startswith('#'):
            header = line.split('\t')
            header[0] = header[0][1:] # strip "#" char from front of header
            samples = [ s.split('.')[0] for s in header[9:] ]
            print("Number of samples in {} = {}".format(in_file, len(samples)))
            # out_handle.write(line)
            headline = "chr\tpos\tref\talt\tancestral\tderived\n"
            out_handle.write(headline)
            continue
        else:
            fields = line.split('\t',9)
            total_vars += 1
            n = len(fields)
            chrom = fields[0]
            start = int(fields[1]) ## 1-based, read coordinates converted for matching
            ref = fields[3]
            alt = fields[4]
            qfilter = fields[6]
            aa_info = fields[7]

            if chrom.startswith('chr'):
                outline_pre = "{}\t{}\t{}\t{}\t".format(chrom,start,ref,alt)
            else:
                outline_pre = "chr{}\t{}\t{}\t{}\t".format(chrom,start,ref,alt)

            ## Check format for first variant
            if total_vars == 1:
                format = fields[8] if n>8 else None
                if format.split(':')[0] != "GT":
                    sys.stderr.write("Error : Invalid format - GT not detected at first position in " + format)
                    sys.exit(-1)
                if not aa_info.startswith("AA"):
                    sys.stderr.write("Error : Ancestral allele information not detected in {}".format(in_file))
                    sys.exit(-1)

            if len(ref)!=1 or len(alt)!=1:
                not_biallelic_snp += 1
                if debug:
                    print("not biallelic")
                    print("ref = {}, alt = {}, {}".format(ref, alt, aa_info))
                continue

            aa = aa_info.replace("AA=","").split("|")[0]
            if len(aa) < 1:
                unknown += 1
                if debug and unknown < 5:
                    print("no ancestral info reported")
                    print("ref = {}, alt = {}, {}".format(ref, alt, aa_info))
                outline = outline_pre + "NA\tNA\n"
                out_handle.write(outline)
                continue
            if aa == ".":
                no_ancestral += 1
                if debug and no_ancestral < 5:
                    print("no ancestral coverage")
                    print("ref = {}, alt = {}, {}".format(ref, alt, aa_info))
                outline = outline_pre + "NA\tNA\n"
                out_handle.write(outline)
            elif aa == "-":
                lineage_specific_ins += 1
                if debug:
                    print("lineage specific insertion")
                    print("ref = {}, alt = {}, {}".format(ref, alt, aa_info))
                outline = outline_pre + "NA\tNA\n"
                out_handle.write(outline)
            elif aa.upper() == alt.upper():
                alt_ancestral += 1

                if aa.isupper():
                    more_confident_ra += 1
                else:
                    less_confident_ra += 1
                outline = outline_pre + "{}\t{}\n".format(alt.upper(),ref.upper())
                out_handle.write(outline)
            elif aa.upper() == ref.upper():
                ref_ancestral += 1

                if aa.isupper():
                    more_confident_ra += 1
                else:
                    less_confident_ra += 1
                outline = outline_pre + "{}\t{}\n".format(ref.upper(),alt.upper())
                out_handle.write(outline)
            elif aa=="N":
                # can't call from ancestral seqs
                ambig_ancestral += 1
                if debug and ambig_ancestral < 5:
                    print("Ancestral seqs disagree")
                    print("ref = {}, alt = {}, {}".format(ref, alt, aa_info))
                outline = outline_pre + "NA\tNA\n"
                out_handle.write(outline)
            else:
                other_ancestral += 1
                if debug and other_ancestral < 5:
                    print("Non ref or alt ancestral allele")
                    print("ref = {}, alt = {}, {}".format(ref, alt, aa_info))
                outline = outline_pre + "NA\tNA\n"
                out_handle.write(outline)

            total_snps += 1

        if (total_vars % 100000 == 0 and verbose):
                print("## " + str(total_vars) + " variants")

    vcf_handle.close()
    out_handle.close()
    if verbose:
        print("# Total variants processed = {} ({}%)".format(total_vars, round((total_vars/total_vars)*100,2)))
        print("#    Not biallelic SNP = {} ({}%)".format(not_biallelic_snp, round((not_biallelic_snp/total_snps)*100,2)))
        print("#    No ancestral info reported = {} ({}%)".format(unknown, round((unknown/total_vars)*100,2)))
        print("# Total biallelic SNPs processed = {} ({}%)".format(total_snps, round((total_snps/total_vars)*100,2)))
        print("#    Ancestral seqs disagree = {} ({}%)".format(ambig_ancestral, round((ambig_ancestral/total_snps)*100,2)))
        print("#    Non ref/alt ancestral allele = {} ({}%)".format(other_ancestral, round((other_ancestral/total_snps)*100,2)))
        print("#    No ancestral coverage = {} ({}%)".format(no_ancestral, round((no_ancestral/total_snps)*100,2)))
        print("#    Lineage-specific inserstions = {} ({}%)".format(lineage_specific_ins, round((lineage_specific_ins/total_vars)*100,2)))
        print("#    Ref ancestral alleles = {} ({}%)".format(ref_ancestral, round((ref_ancestral/total_snps)*100,2)))
        print("#       All ancestral seqs agree = {} ({}%)".format(more_confident_ra, round((more_confident_ra/ref_ancestral)*100,2)))
        print("#       Some ancestral seqs agree = {} ({}%)".format(less_confident_ra, round((less_confident_ra/ref_ancestral)*100,2)))
        print("#    Alt ancestral alleles = {} ({}%)".format(alt_ancestral, round((alt_ancestral/total_snps)*100,2)))
        print("#       All ancestral seqs agree = {} ({}%)".format(more_confident_ra, round((more_confident_ra/alt_ancestral)*100,2)))
        print("#       Some ancestral seqs agree = {} ({}%)".format(less_confident_ra, round((less_confident_ra/alt_ancestral)*100,2)))

    return [total_vars,
            total_snps,
            ref_ancestral,
            alt_ancestral,
            other_ancestral,
            no_ancestral,
            not_biallelic_snp,
            lineage_specific_ins,
            more_confident_ra,
            less_confident_ra,
            unknown]

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--verbose','-v', action='store_true')
    parser.add_argument('invcf', help="Vcf file with ancestral allele information to parse")
    parser.add_argument('outname', help="Name of file to output with vcf records filtered by ancestral-ref/alt matching")

    return parser.parse_args()

def main():
    args = parse_args()
    parse_vcf_ancestralRef( args.invcf, args.outname, verbose=True, debug=False)

if __name__ == '__main__':
    main()
