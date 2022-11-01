#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adapted from the `prep` subcommand of cisvar.py at
https://github.com/TheFraserLab/cisVar/blob/master/cisVar.py
"""

from __future__ import print_function
import os
import sys
import bz2
import gzip
import argparse
from argparse import ArgumentParser
import operator
import subprocess
import multiprocessing as mp

from time import sleep
from datetime import datetime

from collections import defaultdict

# psutil is used to print memory usage if installed, otherwise ignored
try:
    import psutil
except ImportError:
    psutil = None

CORES = mp.cpu_count()

__version__ = '2.0.0b3'

def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.
    Returns already opened files unchanged
    Text mode automatic for compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile

    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'

    # Refuse to handle non-strings that aren't files.
    if not isinstance(infile, str):
        raise ValueError("i cannot open a filename that isn't a string.")

    # Treat '-' appropriately
    if infile == '-':
        if 'w' in mode:
            return sys.stdout
        return sys.stdin

    # If possible open zipped files
    if infile.endswith('.gz'):
        return gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(bz2, 'open'):
            return bz2.open(infile, mode)
        return bz2.bz2file(infile, mode)

    # Fall back on regular open
    return open(infile, mode)

def prep_files(vcf_files, prefix_name, inds=None, limit_file=None,
               alleles_file=None, chrom_format=0, skip_indels=True, cores=4):
    """Create genotype, bed, and allele files for pipeline.
    Bed is used by mpileup to define SNPs to pileup.
    Alleles file is used to set default ref and alt in POST.
    Note, this only needs to be created once per experiment, not once per
    group/pool. The same files can be used for each population/group/pool,
    they are further filtered in the genotype step.
    This code first loads all locations from the limit file (if present) and
    any alleles to override. These two lists are kept in memory. It then loops
    through each vcf_file in parallel, line-by-line, writing temp files, before
    concatenating the results into a single file.
    Params
    ------
    vcf_files : list
        A list of VCF files to parse, should contain genotypes, ref and alt.
    prefix_name : str
        A name for the output files.
    inds : list or path, optional
        Either a list of individuals or a list of paths to individuals files.
        Should include all individuals in this experiment, in all pools.
    limit_file : path, optional
        Path to a bed/vcf/simple file of SNPs to consider. Should be a subset
        of SNPs in the VCF files. BED3 format is fine, simple format is just
        chromosome and position (1-based) to exclude, tab separated.
    alleles_file : path, optional
        Path to a VCF/bed/txt file with alleles, it text, should be
        chr position (base-1). Can be same file as limit file if desired.
    chrom_format : {0, 1, 2}, optional
        0 = leave
        1 = force chr#
        2 = force #
    skip_indels : bool, optional
        Don't include rows with indels
    cores : int, optional
        Number of parallel cores to use on multiple vcf files.
    Writes
    ------
    prefix.genotypes.txt.gz
        A tab separated file, including header with individuals, of:
        chrom, position (base-1), ref (upper), alt (upper), genotypes (0/1/2)
    prefix.locations.bed.gz
        A BED file of all of the locations to be tested, including names
    prefix.individuals.txt.gz
        A newline separated file of individuals, replicates inds if provided,
        overwrites even if already exists.
    """
    if isinstance(vcf_files, str):
        vcf_files = [vcf_files]
    for fl in vcf_files:
        if not os.path.isfile(fl):
            raise OSError('Could not find file {}'.format(fl))
    assert isinstance(prefix_name, str)
    par = False if len(vcf_files) == 1 or cores < 2 else True
    if limit_file:
        sys.stderr.write(
            "Loading test locations from {}\n".format(limit_file)
        )
        limit_locs = parse_loc_file(limit_file)
    else:
        limit_locs = None
    if alleles_file:
        sys.stderr.write(
            "Loading alleles from {}\n".format(alleles_file)
        )
        alleles = parse_loc_file(alleles_file, mode='alleles')
    else:
        alleles = None
    if par:
        pool = mp.Pool(cores)

    # Individuals
    if not inds:
        inds = geno_file(vcf_files[0], get_inds=True)
    elif os.path.isfile(inds):
        with open_zipped(inds) as fl:
            inds = fl.read().strip().split('\n')

    final_inds = list(inds)

    # Write the final individuals, order preserved
    ind_output = '{}.individuals.txt.gz'.format(prefix_name)
    with open_zipped(ind_output, 'w') as fout:
        fout.write('\n'.join(final_inds))

    # Primary parsing
    sys.stderr.write("Begining genotype file parse\n")

    # Header
    header = geno_file(vcf_files[0], get_header=True)

    count = 0
    jobs = []
    geno_output = '{}.genotypes.txt.gz'.format(prefix_name)
    bed_output = '{}.locations.bed.gz'.format(prefix_name)
    for fl in vcf_files:
        ofl = '{}.{}.gz'.format(geno_output, count)
        bfl = '{}.{}.gz'.format(bed_output, count)
        args = (
            fl, ofl, header, bfl, final_inds, limit_locs,
            alleles, skip_indels, chrom_format
        )
        if par:
            jobs.append([ofl, bfl, pool.apply_async(parse_file, args)])
        else:
            parse_file(*args)
            jobs.append((ofl, bfl, None))
        count += 1

    ofls = [i[0] for i in jobs]
    bfls = [i[1] for i in jobs]
    jobs = [i[2] for i in jobs]

    if par:
        for job in jobs:
            job.wait()
        sleep(2)

    for ofl in ofls:
        done_file = ofl + '.done'
        if not os.path.isfile(done_file):
            sleep(1)
            if not os.path.isfile(done_file):
                print('Cannot find {}, assuming failure'.format(done_file))
                raise Exception('Missing done file')
        os.remove(done_file)

    for bfl in bfls:
        done_file = bfl + '.done'
        if not os.path.isfile(done_file):
            sleep(1)
            if not os.path.isfile(done_file):
                print('Cannot find {}, assuming failure'.format(done_file))
                raise Exception('Missing done file')
        os.remove(done_file)

    # Merge
    print('Merging files')
    with open_zipped(geno_output, 'w') as fout:
        fout.write(
            'chrom\tposition\tref\talt\t{}\n'.format('\t'.join(final_inds))
        )
    subprocess.check_call(
        'zcat {} | gzip >> {}'.format(' '.join(ofls), geno_output),
        shell=True
    )
    subprocess.check_call(
        'zcat {} | gzip > {}'.format(' '.join(bfls), bed_output), shell=True
    )

    print('Removing tmp files')
    for temp_file in ofls + bfls:
        if os.path.isfile(temp_file):
            os.remove(temp_file)


def parse_file(geno, outfile, header, bedfile, inds=None, limit_locs=None,
               alleles=None, skip_indels=False, chrom_format=0, log=None,
               fail_if_too_few_inds=True):
    """File parser, see parse_files for information."""
    if not log:
        log = outfile + '.parse.log'
    kept = 0
    indels = 0
    skipped = 0
    superbad = []
    print('Creating {} and {} from {}'.format(outfile, bedfile, geno))
    if chrom_format == 0:
        def chromit(x):
            """Return x."""
            return x
    elif chrom_format == 1:
        def chromit(x):
            """Force x to start with chr."""
            return x if x.startswith('chr') else 'chr' + x
    elif chrom_format == 2:
        def chromit(x):
            """Force x to not start with chr."""
            return x if not x.startswith('chr') else x[3:]
    else:
        raise ValueError('chrom_format must be one of 0, 1, 2')

    if not limit_locs:
        limit_locs = dict()
    if not alleles:
        alleles = dict()

    lg = open_zipped(log, 'a')
    with open_zipped(outfile, 'w') as gn, open_zipped(bedfile, 'w') as bd:
        lg.write('Parse initiated\n')
        for f in geno_file(geno, check_header=header, ind_list=inds, log=lg):
            # Write core records
            if skip_indels and (len(f[2]) > 1 or len(f[3]) > 1):
                indels += 1
                continue
            chrom = chromit(f[0])
            pos = int(f[1])
            # Drop excluded positions
            try:
                if pos not in limit_locs[chrom]:
                    skipped += 1
                    continue
            except KeyError:
                pass

            # Override alleles
            try:
                if pos in alleles[chrom]:
                    ref, alt = alleles[chrom][pos]
                else:
                    ref = f[2]
                    alt = f[3]
            except KeyError:
                ref = f[2]
                alt = f[3]

            name = f[4]
            gt_inds = f[5:]

            outstr = (
                '{chr}\t{pos}\t{ref}\t{alt}'.format(
                    chr=chrom, pos=pos, ref=ref, alt=alt
                )
            )
            # Write genotypes
            if inds:
                gt_ind_l = len(gt_inds)
                ind_l = len(inds)
                if gt_ind_l < ind_l:
                    lg.write(
                        '{chr}.{pos} had too few ({ind}) individuals\n'
                        .format(chr=chrom, pos=f[1], ind=len(gt_inds))
                    )
                    if fail_if_too_few_inds:
                        superbad.append(f)
                    continue
                if gt_ind_l > ind_l:
                    lg.write(
                        '{chr}.{pos} had too many ({ind}) individuals\n'
                        .format(chr=chrom, pos=f[1], ind=len(gt_inds))
                    )
                    superbad.append(f)
                    continue
            gn.write('{}\t{}\n'.format(outstr, '\t'.join(gt_inds)))
            bd.write('{}\t{}\t{}\t{}\n'.format(chrom, pos-1, pos, name))
            kept += 1
    lg.write(
        'Skipped: {}\nIndels: {}\nKept: {}\nParse Complete\n'
        .format(skipped, indels, kept)
    )
    print(
        'Finished {}. Kept {} SNPs. Dropped {} indels and skipped {}'
        .format(outfile, kept, indels, skipped)
    )
    if superbad:
        msg = 'Length failure despite length filtering:\n{}\n'.format(superbad)
        lg.write(msg)
        print(msg)

    print('Writing done files')
    with open(outfile + '.done', 'w') as fout:
        fout.write('Done\n')
    with open(bedfile + '.done', 'w') as fout:
        fout.write('Done\n')

    lg.close()
    return

# Simple location file handling

def parse_loc_file(infile, mode='locations'):
    """Return a dict of SNPs from a bed/vcf/text of chr postition (base-1)."""
    if mode not in ['locations', 'alleles', 'names']:
        raise ValueError(
            'mode must be one of locations, alleles, names. Is {}'
            .format(mode)
        )
    names = infile.split('.')
    if 'vcf' in names:
        fl_tp = 'vcf'
    elif 'bed' in names:
        fl_tp = 'bed'
    else:
        fl_tp = 'txt'
    splt_char = '\t ' if fl_tp == 'txt' else '\t'
    with open_zipped(infile) as fl:
        for line in fl:
            # loc = fl.tell()
            # Comment lines
            if line.startswith('#'):
                # loc = fl.tell()
                continue
            # First non-comment
            else:
                header = line.strip().split(splt_char)
                if line[1].isdigit():
                    # Not a header
                    # fl.seek(loc)
                    pass
                else:
                    continue
                break
        # All non-head lines
        if mode == 'locations':
            res = _get_locs(fl, fl_tp, splt_char)
        elif mode == 'alleles':
            res = _get_alleles(fl, fl_tp, splt_char)
        elif mode == 'names':
            res = _get_names(fl, fl_tp, splt_char)
        else:
            raise ValueError('Invalid mode')
        # Generalize chromosome name
        for k, v in res.items():
            res[flip_chrom(k)] = v
        return res


def _get_locs(fl, tp, sp):
    result = defaultdict(set)
    for line in fl:
        f = line.rstrip().split(sp)
        if tp == 'bed':
            loc = int(f[2])
        else:
            loc = int(f[1])
        result[f[0]].add(loc)
    return result


def _get_alleles(fl, tp, sp):
    if tp == 'bed':
        raise ValueError("Can't get alleles from a bed file")
    result = defaultdict(dict)
    for line in fl:
        f = line.rstrip().split(sp)
        loc = int(f[1])
        if tp == 'vcf':
            ref = f[3]
            alt = f[4]
        else:
            ref = f[2]
            alt = f[3]
        result[f[0]][loc] = (ref, alt)
    return result


def _get_names(fl, tp, sp):
    result = defaultdict(dict)
    for line in fl:
        f = line.rstrip().split(sp)
        if tp == 'bed':
            name = f[3]
            loc = int(f[2])
        else:
            name = f[2]
            loc = int(f[1])
        result[f[0]][loc] = name
    return result


def flip_chrom(chrom):
    """Return opposite of chr# or # for chromosome name."""
    if isinstance(chrom, int):
        chrom = str(chrom)
    if chrom.startswith('chr'):
        return chrom[3:]
    return 'chr' + chrom


# Genotype file handling


def geno_file(infile, check_header=None, ind_list=None,
              get_header=False, get_inds=False, log=sys.stderr):
    """Wrapper for genotype files, currently only vcf
    Parameters
    ----------
    infile : str
        Path to a vcf file, can be zipped.
    check_header : list, optional
        If a list is provided, assert that header matches before iterating.
    ind_list : list of str, optional
        If list provided, filter output to only include these individuals.
    get_header : bool, optional
        Return a header instead of yielding a row
    get_inds : bool, optional
        Return a list of individuals instead of yielding a row
    Yields
    ------
    chr : str
    pos : str (base 1)
    ref : str
    alt : str
    name : str
    inds : list of str
    """
    name = infile.split('.')
    if 'vcf' in name:
        return vcf_file_gt(
            infile, check_header, ind_list, get_header, get_inds, log
        )
    else:
        raise NotImplementedError(
            'Currently only VCFs supported, filetype required in name'
        )


def vcf_file_gt(vcf, check_header=None, ind_list=None,
                get_header=False, get_inds=False, log=sys.stderr):
    """VCF file iterator with goodies
    Parameters
    ----------
    vcf : str
        Path to a vcf file, can be gzipped
    check_header : list, optional
        If a list is provided, assert that header matches before iterating.
    ind_list : list of str, optional
        If list provided, filter output to only include these individuals.
    get_header : bool, optional
        Return a header instead of yielding a row
    get_inds : bool, optional
        Return a list of individuals instead of yielding a row
    bad : file
        File or handle to write errors
    Yields
    ------
    chr : str
    pos : str (base 1)
    ref : str
    alt : str
    name : str
    inds : list of str
    """
    # Get header at genotype loc
    with open_zipped(vcf) as fin:
        line = fin.readline()
        while True:
            if not line.startswith('#'):
                break
            header = line
            pos    = fin.tell()
            line   = fin.readline()

    # Get genotype location
    try:
        gti = line.split('\t')[8].split(':').index('GT')
    except ValueError:
        raise ValueError(
            'GT column missing from VCF format string. In need '
            'the genotype info location to know where to look '
            'for the genotype (looks like 0|1)'
        )

    headers = header.strip().split('\t')
    inds = headers[9:]
    if check_header and header != check_header:
        sys.stderr.write('{} bad header\n'.format(vcf))
        sys.stderr.write('Should be:\n{}\nInd:\n{}\n'.format(check_header, header))
        raise Exception('{} file had an invalid header'.format(vcf))
    if get_header:
        return header
    if get_inds:
        return inds

    return _vcf_file(vcf, inds, ind_list, gti, pos, log)


def _vcf_file(vcf, inds, ind_list, gti, pos, log):
    """Iterator for VCF."""

    # Get individual locations for lookup/filter
    ind_pos = []
    total_inds = len(inds)
    if ind_list:
        ind_count = len(ind_list)
        if sorted(ind_list) != sorted(inds):
            for ind in ind_list:
                ind_pos.append(inds.index(ind))
    else:
        ind_count = len(inds)

    short_lines = 0
    bad_gts = 0
    count = 0

    # Actually parse file
    with open_zipped(vcf) as fin:
        fin.seek(pos)
        count = 0
        for line in fin:
            count += 1
            f = line.strip().split('\t')
            out = [f[0], int(f[1]), f[3], f[4], f[2]]
            gti = f[8].split(':').index('GT')
            ind_data = f[9:]
            if not len(ind_data) == total_inds:
                short_lines += 1
                continue
            if ind_pos:
                these_inds = []
                for loc in ind_pos:
                    these_inds.append(ind_data[loc])
                if len(these_inds) != ind_count:
                    short_lines += 1
                    continue
            else:
                these_inds = ind_data
            for ind in these_inds:
                gt = parse_gt(ind.split(':')[gti])
                if not gt:
                    break
                out.append(gt)
            if (len(out) - 5) != ind_count:
                #  sys.stderr.write(
                    #  'Bad line: {}\t{}\t{}\nLen: {}\tExpected: {}\n'.format(
                        #  count, f[0], f[1], len(out)-4, ind_count
                    #  )
                #  )
                #  raise Exception('Bugger')
                bad_gts += 1
                continue
            yield out
    log.write(
        'Total: {}\nBad Genotypes: {}\nShort Lines (not enough inds): {}\n'
        .format(count, bad_gts, short_lines)
    )
    print(
        '{}: Total: {}, Bad Genotypes: {}, Short Lines (not enough inds): {}'
        .format(vcf, count, bad_gts, short_lines)
    )
    return


def parse_gt(gt):
    """Convert common genotypes into 0,1,2."""
    if '|' in gt:
        if gt == '0|0':
            gt = '0'
        elif gt in ['0|1', '1|0']:
            gt = '1'
        elif gt == '1|1':
            gt = '2'
        else:
            gt = None
    elif '/' in gt:
        if gt == '0/0':
            gt = '0'
        elif gt in ['0/1', '1/0']:
            gt = '1'
        elif gt == '1/1':
            gt = '2'
        else:
            gt = None
    elif gt.isdigit():
        gt = gt
    else:
        gt = None
    return gt

def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        '-F', '--SampleName', dest='prefix_name', default='AA_admixed',
        help='sample/population name (default: AA_admixed)'
    )
    parser.add_argument(
        '-i', '--all-inds',
        help='File of individuals in all groups, one per line'
    )
    parser.add_argument(
        '-l', '--limit-file',
        help='BED/VCF/txt file of SNPs to consider'
    )
    parser.add_argument(
        '-a', '--allele-file',
        help='BED/VCF/txt file of alleles to override VCF allels (subset of vcf)'
    )
    parser.add_argument(
        '--chrom-format', choices={'chr', 'num', 'ignore'}, default='ignore',
        help='chr: make format "chr#", num: make format "#", ' +
        'ignore: do nothing (default: ignore)'
    )
    parser.add_argument(
        '--include-indels', action='store_false',
        help='Do not skip indels'
    )
    parser.add_argument(
        '-c', '--cores', type=int, default=CORES,
        help='Number of cores to use (default: all)'
    )
    parser.add_argument(
        'vcf_files', nargs='+', help='VCF files with genotypes'
    )

    return parser.parse_args()

def main(argv=None):
    args = parse_args()
    if args.chrom_format == 'ignore':
            chrom_format = 0
    elif args.chrom_format == 'chr':
        chrom_format = 1
    elif args.chrom_format == 'num':
        chrom_format = 2
    prep_files(
        args.vcf_files, args.prefix_name, inds=args.all_inds,
        limit_file=args.limit_file, alleles_file=args.allele_file,
        chrom_format=chrom_format, skip_indels=args.include_indels,
        cores=args.cores
    )

if __name__ == '__main__':
    sys.exit(main())
