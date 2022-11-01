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

// Nextflow pipeline to find trans eQTL of target gene set
// mean expression.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

// process 1000G genotypes for Fst

Channel.of( 1..22 )
   .map{ n -> tuple("chr"+n, file("~/data/gzvcfs/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")) } // vcf_input must have ${n} in the name
   .set{ FST_VCF }

Channel.fromPath(params.fst_pop1)
   .map{ file -> tuple(file.simpleName, file) }
   .set{ FST_POP1 }
Channel.fromPath(params.fst_pop2)
   .map{ file -> tuple(file.simpleName, file) }
   .set{ FST_POP2 }
FST_POP1
   .merge(FST_POP2)
   .combine(FST_VCF)
   .set{ FST }

process weir_fst{
  label 'vcf'
  publishDir "${params.outdir}/fst_all/"

  when:
  params.mode =~ /(all)/

  input:
  tuple npop1, file(pop1), npop2, file(pop2), chr, file(vcffile) from FST

  output:
  tuple chr, file(fstout) into FST_HG19

  script:
  fstbase = "${npop1}_vs_${npop2}.${chr}"
  fstout = "${fstbase}.fst.bed"
  """
  vcftools \
    --gzvcf $vcffile \
    --weir-fst-pop $pop1 \
    --weir-fst-pop $pop2 \
    --out $fstbase; \
  tail -n +2 ${fstbase}.weir.fst \
  | awk '\$1 ~ /^chr/ {print \$0;next} {print "chr" \$0}' \
  | awk -F \$'\t' '{ if(\$3 != "-nan" && \$3 > 0.0) {print \$1, \$2-1, \$2, \$3 }}' OFS='\t' \
  > $fstout
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

FST_HG19
    .combine( BED_CHR, by: [0] )
    .set{ FST_SUB }

process fst_subset{
  label 'bed'

  when:
  params.mode =~ /(all)/

  input:
  tuple chr, path(fst), path(bed) from FST_SUB

  output:
  path(outname) into FST_SUBBED

  script:
  outname = "${fst.simpleName}.${bed.simpleName}.fst.${chr}"
  """
  bedtools intersect \
  -u \
  -a $fst \
  -b $bed \
  > $outname
  """
}

process fst_concat{
  label 'bed'
  publishDir "${params.outdir}/fst_sub/"

  when:
  params.mode =~ /(all)/

  input:
  path(fst) from FST_SUBBED.collect().dump(tag: 'concat_in')

  output:
  path(outname) into FST_RAW

  script:
  fst_chroms = fst.findAll{ it.toString() }.sort()
  outname = "${fst[0].baseName}.txt.gz"
  """
  cat ${fst_chroms.join(' ')} \
  | bgzip -c \
  > $outname
  """
}

// r process for adding percentile columns to FST
