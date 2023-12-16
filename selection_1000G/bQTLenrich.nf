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

// Nextflow pipeline to perform enrichment test on list of ChIP QTLs

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

// perform enrichment test on list of ChIP QTLs

Channel.fromList(params.qtl_fnames)
  .flatten()
  .combine( Channel.fromPath(params.d_all_fname) )
  .combine( Channel.fromPath(params.gs_fname) )
  .combine( Channel.fromPath(params.frq_fname) )
  .combine( Channel.fromPath(params.fst_fname) )
  .combine( Channel.fromPath(params.ihs_fname) )
  .combine( Channel.fromPath(params.chain_fname) )
  .combine( Channel.fromPath(params.ref_fname) )
  .combine( Channel.fromPath(params.eqtl_lfsr_fname) )
  .combine( Channel.fromPath(params.eqtl_es_fname) )
  .combine( Channel.fromPath(params.de_fname) )
  .combine( Channel.fromPath(params.rde_fname) )
  .combine( Channel.fromPath(params.rdr_fname) )
  .combine( Channel.fromPath(params.teqtl_fname) )
  .set{ QTL }

process qtl_overlap{
  label 'r'
  publishDir "${params.outdir}/overlap/"

  when:
  params.mode =~ /(all)/

  input:
  tuple path(qfname), path(d), path(gs), path(frq), path(fst), path(ihs), path(ch), path(ref), path(elf), path(ees), path(de), path(re), path(rr), path(t) from QTL

  output:
  tuple path("*.fisherEnrichments.txt"), path("*.fstWilcoxData.txt.gz"), path("*.fstWilcoxBinned.txt.gz"), path("*.ihsWilcoxData.txt.gz")  into CONCAT
  path("${name}.${qn}.txt.gz") into TOP_CANDS

  script:
  name = d.baseName
  qn = qfname.simpleName
  """
  Rscript ${params.bin}/diffscore_QTL_overlap.R \
  --qtl_fname $qfname \
  --qtl_name $qn \
  --d_all_fname $d \
  --gs_fname $gs \
  --frq_fname $frq \
  --fst_fname $fst \
  --ihs_fname $ihs \
  --chain_fname $ch \
  --ref_fname $ref \
  --eqtl_lfsr_fname $elf \
  --eqtl_es_fname $ees \
  --de_fname $de \
  --rde_fname $re \
  --rdr_fname $rr \
  --teqtl_fname $t \
  --class_resolution ${params.class_resolution} \
  --name $name \
  --outdir .
  """
}

process concat_plot{
  label 'r'
  publishDir "${params.outdir}/aggregate/", pattern: "*.fisherEnrichmentsAgg.txt"
  publishDir "${params.outdir}/aggregate/", pattern: "*.diffQTLoverlap.fstWilcoxDataAgg.txt.gz"
  publishDir "${params.outdir}/aggregate/", pattern: "*.diffQTLoverlap.fstWilcoxBinnedAgg.txt.gz"
  publishDir "${params.outdir}/aggregate/", pattern: "*.diffQTLoverlap.ihsWilcoxDataAgg.txt.gz"
  publishDir "${params.outdir}/aggregate/plots/", pattern: "*.png"

  when:
  params.mode =~ /(all)/

  input:
  path(efs) from CONCAT.flatten().collect()

  output:
  path("*.fisherEnrichmentsAgg.txt")
  path("*.diffQTLoverlap.fstWilcoxDataAgg.txt.gz")
  path("*.diffQTLoverlap.fstWilcoxBinnedAgg.txt.gz")
  path("*.diffQTLoverlap.ihsWilcoxDataAgg.txt.gz")
  path("*.png")

  script:
  name = efs[0].getName().split(".txt")[0]
  """
  Rscript ${params.bin}/concatPlotDiffQTL.R \
  --plotdir ${params.plotdir} \
  --fisher_pattern fisherEnrichments.txt \
  --wilcox_pattern fstWilcoxData.txt.gz \
  --binned_pattern fstWilcoxBinned.txt.gz \
  --ihs_pattern ihsWilcoxData.txt.gz \
  --agg_dir . \
  --name $name \
  --outdir .
  """
}
