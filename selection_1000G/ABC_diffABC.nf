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

// Nextflow pipeline to define candidate CREs based on ABC score thresholds
// and filter ABC component score types to these E-G pairs from previous
// ABC score prediction output.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

Channel.of( "HiChIP" )
  .combine( Channel.of(".meanQN") )
  .combine( Channel.fromPath( params.mean_qn_dir ) )
  // fromPath reads in already flattened, so don't need to transpose
  // // flatten list of prediciton files into separate channels
  // // will be reversed before differential_abc
  // .transpose(by: 2) // reverse of groupTuple(by: [0,1])
  .map{ it -> [ it[0], it[1].split('\\.')[1], it[2] ] }
  .dump(tag: "analPreGroup")
  .groupTuple(by: [0,1])
  .dump(tag: "analGrouped")
  .combine( Channel.fromList( params.score_types ) )
  .set{ DIFFABC }

process differential_abc {
  label 'diff'
  publishDir "${params.outdir}/diffABC/${ctype}/", pattern: "*.txt.gz"
  publishDir "${params.outdir}/diffABC/${ctype}/plots", pattern: "*.png"

  when:
  params.mode =~ /(all)/

  input:
  tuple ctype, qn, path(preds), score_col from DIFFABC

  output:
  tuple ctype, qn, score_col, path("${outbase}.${score_col}.txt.gz") into OVERLAP_PRE_FRESH
  path("*.scoresAllSamps.txt.gz")
  path("*.scoresAllSampsFilt.txt.gz")
  path("*.png")

  script:
  outbase = "allZerosFilt.${qn}.${params.min_nonzero_hic}.AFR_EUR.diff"
  """
  Rscript ${params.bin}/diffABC.R \
  --abc_positive_ccres \
  --prediction_dir . \
  --qnorm_suffix $qn \
  --min_nonzero_hic ${params.min_nonzero_hic} \
  --enh_threshold ${params.enhancer_threshold} \
  --promoter_threshold ${params.promoter_threshold} \
  --score_column $score_col \
  --plotdir ${params.plotdir} \
  --name $outbase \
  --outdir .
  """

}
