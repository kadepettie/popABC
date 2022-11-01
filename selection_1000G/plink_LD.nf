#!/usr/bin/env nextflow
// Copyright (C) 2020 Kade Pettie

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

// Nextflow pipeline to generate plink bed files based on individual
// genotypes for finding variants in LD within or across populations.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

Channel.of( 1..22 )
   .map{ n -> tuple("chr"+n, file("~/data/gzvcfs/ALL.chr${n}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")) } // vcf_input must have ${n} in the name
   .set{ VCFIN }

Channel.fromPath( params.sample_info )
  .combine( Channel.fromList( params.ld_pops ) )
  .set{ LD_IDS }


process get_LD_inds{
  label 'r'
  publishDir "${params.outdir}/inds/"

  when:
  params.mode =~ /(all)/

  input:
  tuple path(si), pops from LD_IDS

  output:
  path(name) into VCF_INDS

  script:
  name = "${pops}.individuals.txt"
  """
  Rscript ${params.bin}/get_1kg_pop_individuals.R \
  --sample_info $si \
  --pops $pops \
  --name $name
  """

}

VCFIN
   .combine( VCF_INDS )
   .set{ VCF_PARSE }


process filter_vcf{
  label 'vcf'

  when:
  params.mode =~ /all/

  input:
  tuple chr, path(vcf), path(inds) from VCF_PARSE

  output:
  tuple chr, path(outname) into PLINK_BED

  script:
  outname = "${inds.simpleName}.${chr}.vcf.gz"
  """
  vcftools \
    --gzvcf $vcf \
    --keep $inds \
    --recode \
    --stdout \
    | bgzip -c \
    > $outname
  """
}

process plink_makebed{
  label 'vcf'
  publishDir "${params.outdir}/plink/"

  when:
  params.mode =~ /all/

  input:
  tuple chr, path(vcf) from PLINK_BED

  output:
  tuple "${outpre}.bed", "${outpre}.bim", "${outpre}.fam" into OLAPS

  script:
  outpre = "${vcf.simpleName}.${chr}"
  """
  ${params.plink}/plink \
    --vcf $vcf \
    --make-bed \
    --out $outpre
  """
}
