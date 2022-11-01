#!/usr/bin/env nextflow
// Copyright (C) 2019 Kade Pettie

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Nextflow pipeline to calculate allele frequencies in population subsets of
// vcf files for variants within pre-defined regions.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

// process 1000G genotypes for allele freqs by pop

Channel.of( 1..22 )
   .map{ n -> tuple("chr"+n, file("~/data/gzvcfs/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")) } // vcf_input must have ${n} in the name
   .combine( Channel.fromList(params.inds) )
   .set{ FRQ }

process allele_freq{
  label 'vcf'
  publishDir "${params.outdir}/frq_all/"

  when:
  params.mode =~ /(all)/

  input:
  tuple chr, path(vcffile), path(inds) from FRQ

  output:
  tuple chr, indname, file(frqout) into FRQ_ALL

  script:
  indname = inds.simpleName
  frqout = "${indname}.${chr}"
  """
  vcftools \
    --gzvcf $vcffile \
    --keep $inds \
    --recode \
    --recode-INFO AA \
    --stdout \
  | vcftools \
    --vcf - \
    --freq \
    --stdout \
  | grep -v 'CHROM' \
  | awk '\$1 ~ /^chr/ {print \$0;next} {print "chr" \$0}' \
  | awk -F \$'\t' ' {print \$1, \$2-1, \$2, \$5, \$6 }' OFS='\t' \
  > $frqout
  """
}

Channel.fromPath(params.region_bed)
    .set{ BED_ALL }

process split_bed{
  label 'base'
  publishDir "${params.outdir}/bed_chr/"

  when:
  params.mode =~ /(all)/

  input:
  path(bed) from BED_ALL

  output:
  path("${outbase}.bed.*") into BED_CHR_PRE

  script:
  outbase = bed.baseName
  """
  for chr in `cut -f 1 $bed | sort | uniq`; \
  do echo \$chr; \
  grep -w \$chr $bed \
  > ${outbase}.bed.\${chr}; \
  done
  """
}

BED_CHR_PRE
    .flatten()
    .map{ it -> [ it.extension, it ] }
    .set{ BED_CHR }

FRQ_ALL
    .combine( BED_CHR, by: [0] )
    .set{ FRQ_SUB }

process frq_subset{
  label 'bed'

  when:
  params.mode =~ /(all)/

  input:
  tuple chr, indname, path(frq), path(bed) from FRQ_SUB

  output:
  tuple indname, path(outname) into FRQ_SUBBED

  script:
  outname = "${indname}.${bed.simpleName}.frq.${chr}"
  """
  bedtools intersect \
  -u \
  -a $frq \
  -b $bed \
  > $outname
  """
}

FRQ_SUBBED
  .groupTuple(by: [0])
  .set{ FRQ_CONCAT }

process frq_concat{
  label 'bed'
  publishDir "${params.outdir}/frq_sub/"

  when:
  params.mode =~ /(all)/

  input:
  tuple indname, path(frq) from FRQ_CONCAT

  output:
  path(outname) into FRQ_RAW

  script:
  frq_chroms = frq.findAll{ it.toString() }.sort()
  outname = "${indname}.RefAltAllele.frq.gz"
  """
  cat ${frq_chroms.join(' ')} \
  | bgzip -c \
  > $outname
  """
}
