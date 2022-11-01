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

// Nextflow pipeline to calculate EUR vs AFR 1000 Genomes Fst for a given
// set of SNPs.



Channel.of( 1..22 )
   .map{ n -> tuple("chr"+n, file("/home/kpettie/data/gzvcfs/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")) } // vcf_input must have ${n} in the name
   .set{ VCFIN }

Channel.fromList(params.snps)
  .combine( Channel.fromPath( params.gene_set ) )
  .combine( Channel.fromPath( params.liftover ) )
  .set{ SNP_PREP }

process prep_snps{
 label 'R'

 when:
 params.mode =~ /(all)/

 input:
 tuple path(snps), path(gs), path(ch) from SNP_PREP

 output:
 path("*.positions.txt") into FST_SNPS

 script:
 outbase = snps.baseName
 """
 Rscript ${params.bin}/fst_prep.R \
 --gene_set $gs \
 --snps $snps \
 --liftover $ch \
 --script_dir ${params.bin} \
 --name $outbase
 """
}

FST_SNPS
  .combine(VCFIN)
  .set{ FST_SUB }

process subset_vcf{
  label 'vcf'

  when:
  params.mode =~ /(all)/

  input:
  tuple path(snps), chr, path(vcffile) from FST_SUB

  output:
  tuple chr, path(outname) into FST_CHROM

  script:
    outname = "${snps.simpleName}.vcf.gz"
    """
    vcftools \
      --gzvcf $vcffile \
      --positions-overlap $snps \
      --recode \
      --stdout \
      | bgzip -c \
      > $outname
    """
}

FST_CHROM
  .combine( Channel.fromPath( params.fst_pop1 ) )
  .combine( Channel.fromPath( params.fst_pop2 ) )
  .set{ FST_CALC }

process fst_by {
  label 'vcf'

  when:
  params.mode =~ /(all)/

  input:
  tuple chr, path(vcffile), path(pop1), path(pop2) from FST_CALC

  output:
  path "${outbase}.weir.fst" into FST_CONCAT

  script:
  outbase = "${vcffile.simpleName}_${pop1.simpleName}_${pop2.simpleName}.${chr}"
  """
  vcftools \
    --gzvcf $vcffile \
    --weir-fst-pop $pop1 \
    --weir-fst-pop $pop2 \
    --out ${outbase}_header; \
  cat ${outbase}_header.weir.fst \
  | tail -n +2 \
  > ${outbase}.weir.fst; \
  rm ${outbase}_header.weir.fst
  """
}

process concat_fst_chroms{
  label 'base'
  publishDir "${params.outdir}/fst/"

  when:
  params.mode =~ /all/

  input:
  path(fsts) from FST_CONCAT.collect()

  output:
  path(outname) // into FST_ANAL

  script:
  fst_chroms = fsts.findAll{ it.toString().endsWith('.weir.fst') }.sort()
  outname = fst_chroms[0].simpleName+".fst.gz"
  """
  cat ${fst_chroms.join(' ')} \
  | bgzip -c \
  > $outname
  """
}
