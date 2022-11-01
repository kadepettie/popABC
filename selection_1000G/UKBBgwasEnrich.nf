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

// Nextflow pipeline to perform enrichment test on multiple bQTL types
// in differential CREs for GWAS phenotypes and matching directionality in
// UK biobank

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

// R script for writing select UKBB phenotype codes from manifest file
// to .txt file with one phenotype code per line to use for downloading
// only this subset
// 'getUKBBphenoCodes.R'

Channel.fromPath( params.ukbbmani_fname )
  .set{ GETPHENOS }

process get_phenotypes{
  label 'r'
  publishDir "${params.outdir}/manifest/", pattern: "*.txt"

  when:
  params.mode =~ /(all)/

  input:
  path(mani) from GETPHENOS

  output:
  path("phenotype_codes.txt") into PHENOTXT

  script:
  """
  Rscript ${params.bin}/getUKBBphenoCodes.R \
  --ukbbmani_fname $mani \
  --outdir .
  """
}

// Split lines of phenotype codes .txt file into separate values in channel for
// parallel downloading and subsequent enrichment tests
PHENOTXT
    .splitText( by:1 )
    // splitText by line adds an extra newline character, which is specifically
    // noted in the documentation, but not sure the purpose of this behavior
    .map{ it -> it - ~/\n$/ }
    .combine( Channel.fromPath( params.ukbbmani_fname ) )
    .set{ PHENOCODES }

// R script for downloading phenotype files
// Use storeDir for downloaded phenotype files since 1000 of them will be ~.5T
// 'downloadUKBBphenotype.R'

process download_results{
  label 'r'
  storeDir params.ukbb_storedir

  when:
  params.mode =~ /(all)/

  input:
  tuple pcode, path(mani) from PHENOCODES

  output:
  path("${pcode}${params.ukbb_suffix}") into DOWNLOAD_OUT

  script:
  """
  Rscript ${params.bin}/downloadUKBBphenotype.R \
  --ukbbmani_fname $mani \
  --pheno_code $pcode \
  --outdir .
  """
}

// first process for LD-expanding bQTL in diff-CREs for GWAS overlap and storeDir
// for these since LD-expanding ~10000 bQTL without parallelizing takes a few hours
// older version of script that LD explands, then performs GWAS enrichment test
// on just one bQTL without directionality analysis
// 'diffscore_gwasQTL.R'

// skip LD-expand process for now, and run directly from pre-LD expanded
Channel.fromPath( params.ldbqtl_grob )
  .collect()
  .map{ it -> [it] }
  .set{ LDBQTL }

Channel.fromPath( params.fst_grob )
  .collect()
  .map{ it -> [it] }
  .set{ FST }

Channel.fromPath( params.frq_grob )
  .collect()
  .map{ it -> [it] }
  .set{ FRQ }

DOWNLOAD_OUT
  .combine( LDBQTL )
  .combine( FST )
  .combine( FRQ )
  .combine( Channel.fromPath( params.ukbbmani_fname ) )
  .set{ GWAS_IN }

// R script for aggregating multiple LD-expanded bQTL-CRE overlaps into single
// table, then perform GWAS hit-level enrichment test for phenotype from
// downloadUKBBphenotype.R in diff-CREs. Include hits in LD with diff-CRE bQTL
// as output for aggregatting into directionality and FST enrichment tests
// 'combinedBqtlGwas.R'

process gwas_overlap{
  label 'overlap'
  publishDir "${params.outdir}/overlap/individual_phenotypes/", pattern: "*_CREs.txt.gz"
  publishDir "${params.outdir}/overlap/individual_phenotypes/", pattern: "*_CREs.testResults.txt"
  publishDir "${params.outdir}/overlap/individual_phenotypes/plots/", pattern: "*.png"
  publishDir "${params.outdir}/overlap/individual_phenotypes/plots/", pattern: "*.svg"

  when:
  params.mode =~ /(all)/

  input:
  tuple path(gw), path(cres), path(fst), path(frq), path(mani) from GWAS_IN

  output:
  tuple path("*_CREs.testResults.txt"), path("*_CREs.txt.gz") into GWAS_OUT
  path("*.png")
  // path("*.svg")

  script:
  """
  Rscript ${params.bin}/combinedBqtlGwas.R \
  --fst_dir . \
  --fst_pattern ${params.fst_pattern} \
  --fst_prefix ${params.fst_prefix} \
  --frq_dir . \
  --aprefix ${params.aprefix} \
  --eprefix ${params.eprefix} \
  --d_dir . \
  --d_pattern ${params.ldbqtl_pattern} \
  --gw_fname $gw \
  --ukbbmani_fname $mani \
  --outdir .
  """
}

// R script for concatting diff-CRE bQTL-GWAS enrichment test results a
// and plotting, and also concatting the hits from from phenotypes with nominally
// significant enrichments for testing if risk allele frequencies match any
// ancestry-associated prevalance for these traits
// 'concatPlotGwasRiskAllele.R'

GWAS_OUT
  .collect()
  .map{ it -> [it] }
  .set{ GWAS_CONCAT }

process gwas_concat_risk{
  label 'concat'
  publishDir "${params.outdir}/topcands/", pattern: "*.txt.gz"
  publishDir "${params.outdir}/topcands/", pattern: "*.txt"
  publishDir "${params.outdir}/topcands/plots/", pattern: "*.png"

  when:
  params.mode =~ /(all)/

  input:
  path(gwres) from GWAS_CONCAT

  output:
  path("*.txt.gz")
  path("*.txt")
  path("*.png")
  // path("*.svg")

  script:
  """
  Rscript ${params.bin}/concatPlotGwasRiskAllele.R \
  --t_dir . \
  --t_pattern _CREs.testResults.txt \
  --r_pattern _CREs.txt.gz \
  --outdir .
  """
}
