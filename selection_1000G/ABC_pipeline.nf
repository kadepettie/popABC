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

// Nextflow pipeline to make Activity-By-Contact enhancer-gene predictions
// using ATAC-seq + HiChIP data normalized for comparison between groups.

out = file(params.outdir)
if( !out.exists() ) out.mkdir()

///////// DOWNLOAD DATA OPTION ////////

process clone_abc {
  label 'git'
  storeDir params.git_dir

  when:
  params.mode =~ /(all)/

  input:
  val(gl) from Channel.of( params.download.abc_git )

  output:
  path("ABC-Enhancer-Gene-Prediction/**") into BIN1

  script:
  """
  git clone $gl
  """

}

process clone_genrich {
  label 'git'
  storeDir params.git_dir

  when:
  params.mode =~ /(all)/

  input:
  val(gl) from Channel.of( params.download.genrich_git )

  output:
  path("Genrich/**") into BIN2

  script:
  """
  git clone $gl; \
  cd Genrich; \
  make; \
  cd ..
  """

}

process rclone_atac {
  label 'rclone'
  storeDir params.outdir + "/atac_bams/"
  storeDir params.outdir + "/atac_bams/"

  when:
  !params.atac_bamdir

  input:
  val(gdir) from Channel.of( params.download.gdrive_dir )
  val(tdroot) from Channel.of( params.download.tardir_root )

  output:
  path("*.bam") into BAM_FROM_DOWNLOAD
  path("*.bam.bai")

  script:
  """
  rclone copyto \
  --tpslimit 1 \
  --transfers 1 \
  ${params.download.gdrive_name}:${gdir} \
  .; \
  ls *.tgz | xargs -i tar -xzvf {}; \
  mv ${tdroot}*/* .; \
  rmdir ${tdroot}*/
  """

}


///////// PEAK CALLING INPUT //////////

(BAM_LOCAL,
  BAM_DOWNLOAD) = ( params.atac_bamdir==""
    ? [ Channel.empty(),
        BAM_FROM_DOWNLOAD.flatten() ]
    : [ Channel.fromPath( params.atac_bamdir + '/' + params.atac_glob + '.bam' ),
        Channel.empty() ]
  )

BAM_LOCAL
  .mix(BAM_DOWNLOAD)
  .branch{
    g1: it =~ Eval.me( params.group_regex[0] )
      return [ params.group_names[0], it ]
    g2: it =~ Eval.me( params.group_regex[1] )
      return [ params.group_names[1], it ]
  }
  .set{ BAM }


BAM.g1
  .concat(BAM.g2)
  .set{ BAM_TO_NSORT }

// must ensure bams are namesorted before calling peaks with genrich
outbamnsort = file(params.outdir + "/atac_bams/namesort/")
if( !outbamnsort.exists() ) outbamnsort.mkdirs()

process namesort {
  label 'samtools'
  storeDir params.outdir + "/atac_bams/namesort/"

  when:
  params.mode =~ /(all)/

  input:
  tuple gn, path(bam) from BAM_TO_NSORT

  output:
  tuple gn, path(outname) into BAM_NSORTED

  script:
  outname = "${bam.baseName}.nsort.bam"
  """
  samtools sort \
  -n \
  -@ 8 \
  -o $outname \
  $bam
  """

}

BAM_NSORTED
  .tap{ BSI }
  .groupTuple(by: 0)
  .tap{ BAMM }
  .combine( Channel.fromPath( params.chrsize ) )
  .combine( BIN1.combine(BIN2) )
  .map{ it -> [ it[0], it[1], it[2] ] } // add channels with code so it gets cloned and stored but don't link it in each workdir
  .set{ BAMS }

//////// cCRE INPUT ////////

Channel.fromPath( params.chrsize )
  .combine( Channel.fromPath(params.blocklist) )
  .combine( Channel.fromPath(params.promoters) )
  .into{ CRE; COMB_CRE }


/////// PROCESSES ////////

process call_peaks {
  label 'genrich'
  conda params.abc_code + "/macs.yml"
  publishDir "${params.outdir}/peaks/", pattern: "*.narrowPeak.sorted"
  publishDir "${params.outdir}/peaks/", pattern: "*.bedgraph"

  memory { 24.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
  maxRetries 5

  when:
  params.mode =~ /(all)/

  input:
  tuple gn, path(bams), path(sizes) from BAMS

  output:
  tuple gn, path("*.narrowPeak.sorted") into PEAKS, COMB_CCRE_PEAKS
  path("*.bedgraph")

  script:
  bamstring = bams.findAll{ it.toString().endsWith('.bam') }.sort()
  """
  ${params.genrich_bin}/Genrich -t ${bamstring.join(',')} \
  -o ${gn}.narrowPeak \
  -v \
  -f ${gn}.pileup.qvals.bedgraph \
  -k ${gn}.pileup.pvals.bedgraph \
  -y \
  -j \
  -d 151 \
  -p 0.1; \
  bedtools sort -faidx $sizes \
  -i ${gn}.narrowPeak \
  > ${gn}.narrowPeak.sorted
  """

}

process samtools_merge {
  label 'samtools'
  publishDir "${params.outdir}/bams/", pattern: "*.bam"

  when:
  params.mode =~ /(all)/

  input:
  tuple gn, path(bams) from BAMM

  output:
  tuple gn, path("*_merged.bam") into MBAM, BIDX

  script:
  bamstring = bams.findAll{ it.toString().endsWith('.bam') }.sort()
  """
  samtools merge \
  --threads ${params.samcores} \
  ${gn}_merged.unsort.bam \
  ${bamstring.join(' ')}; \
  samtools sort \
  --threads ${params.samcores} \
  -o ${gn}_merged.bam \
  ${gn}_merged.unsort.bam
  """

}

process samtools_index {
  label 'samtools'
  publishDir "${params.outdir}/bams/", pattern: "*.bam.bai"

  when:
  params.mode =~ /(all)/

  input:
  tuple gn, path(bam) from BIDX

  output:
  tuple gn, path("*.bam.bai") into BAI

  script:
  """
  samtools index \
  -@ ${params.samcores} \
  $bam
  """

}

CRE
  .combine( PEAKS
              .combine(MBAM, by: 0)
              .combine(BAI,  by: 0) )
  .set{ MCRE }

process call_ccres {
  label 'ccres'
  publishDir "${params.outdir}/cCREs/", pattern: "*.Counts.bed"
  conda params.abc_code + "/abcenv.yml"

  when:
  params.mode =~ /(all)/

  input:
  tuple path(sizes), path(excl), path(incl), gn, path(peaks), path(bam), path(idx) from MCRE

  output:
  // using interleaving method to get even contribution to top N peaks from both groups
  // candidate regions generated here no longer used
  tuple gn, path("*.candidateRegions.bed") // into GCRE
  tuple gn, path("*.Counts.bed") into COMB_CCRE_COUNTS

  script:
  """
  python ${params.abc_code}/src/makeCandidateRegions.py \
  --narrowPeak $peaks \
  --bam $bam \
  --outDir ./ \
  --chrom_sizes $sizes \
  --regions_blocklist $excl \
  --regions_includelist $incl \
  --peakExtendFromSummit 250 \
  --nStrongestPeaks 150000
  """

}

COMB_CCRE_PEAKS
  .combine(COMB_CCRE_COUNTS, by: 0)
  .collect()
  .combine( COMB_CRE )
  .combine( Channel.fromPath(params.chrsizebed) )
  .set{ COMB_CCRE }

// NEW PROCESS
process combine_group_ccres {
  label 'r'
  publishDir "${params.outdir}/cCREs/", pattern: "*.candidateRegions.bed"
  publishDir "${params.outdir}/cCREs/", pattern: "*.Counts.txt"

  when:
  params.mode =~ /(all)/

  input:
  tuple gn1, path(p1), path(c1), gn2, path(p2), path(c2), path(sizes), path(excl), path(incl), path(sizesbed) from COMB_CCRE

  output:
  path(outname) into COMB_TOP_CCRES
  path("combinedTopCcres.Counts.txt")

  script:
  // TODO: get suffixes from input filenames
  outname = "${gn1}_${gn2}.candidateRegions.bed"
  """
  Rscript ${params.bin}/interleave_topN_peaks.R \
  --genome_build hg19 \
  --peaks_suffix .narrowPeak.sorted \
  --counts_suffix _merged.bam.Counts.bed \
  --group1_name $gn1 \
  --group2_name $gn2 \
  --group_type ancestry \
  --peakExtendFromSummit 250 \
  --nStrongestPeaks 150000 \
  --name combinedTopCcres; \
  bedtools sort \
  -i combinedTopCcres.bed \
  -faidx $sizes \
  | bedtools intersect \
  -v \
  -wa \
  -a stdin \
  -b $excl \
  | cut -f 1-3 \
  | (bedtools intersect -a $incl -b $sizesbed -wa | cut -f 1-3 && cat) \
  | bedtools sort \
  -i stdin \
  -faidx $sizes \
  | bedtools merge \
  -i stdin \
  > $outname
  """

}

process sort_index {
  label 'samtools'
  publishDir "${params.outdir}/bams/", pattern: "*.psort.bam"
  publishDir "${params.outdir}/bams/", pattern: "*.psort.bam.bai"

  when:
  params.mode =~ /(all)/

  input:
  tuple gn, path(bam) from BSI

  output:
  tuple sn, path("*.bam"), path("*.bam.bai") into CBAM // no need to track group name due to combined cCREs

  script:
  outbam = bam.baseName + '.psort.bam'
  sn = bam.simpleName
  """
  samtools sort \
  --threads ${params.samcores} \
  -o $outbam \
  $bam; \
  samtools index \
  -@ ${params.samcores} \
  $outbam
  """

}

COMB_TOP_CCRES
  .combine( CBAM ) // single, unbiased cCRE set for comparison between groups
  .combine( Channel.fromPath(params.genes) )
  .combine( Channel.fromPath(params.expression) )
  .combine( Channel.fromPath(params.chrsize) )
  .combine( Channel.fromPath(params.chrsizebed) )
  .combine( Channel.fromPath(params.ubiquitous) )
  .set{ CTC }

( CTCL,
  CATAC_NO_QN ) = ( params.qnorm_refs==""
                    ? [ [Channel.empty(),
                         Channel.empty()],
                        CTC ]
                    : [ CTC.into(2),
                        Channel.empty() ] )

QN_REFS_UNFILT = CTCL[0]
TO_NORM = CTCL[1]

QN_REFS_UNFILT
  .dump(tag: 'qn_pre')
  .filter{ it[1] in params.qnorm_refs } // test if this works
  .dump(tag: 'qn_filt')
  .set{ QN_REFS }

process atac_activity_qnorm_ref {
  label 'ccres'
  publishDir "${params.outdir}/neighborhoods/${snpdx}/"
  conda params.abc_code + "/abcenv.yml"

  memory { 2.GB * task.attempt }
  time { 1.hour * task.attempt }
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 5

  when:
  params.mode =~ /(all)/

  input:
  tuple path(ccres), sn, path(bam), path(idx), path(g), path(e), path(s), path(sb), path(u) from QN_REFS

  output:
  tuple snpdx, path("EnhancerList.txt"), path("GeneList.txt") into PDXQN
  tuple sn, path("EnhancerList.txt") into MAKE_QN

  script:
  snpdx = "${sn}.${sn}QN"
  """
  python ${params.abc_code}/src/run.neighborhoods.py \
  --candidate_enhancer_regions $ccres \
  --genes $g \
  --ATAC $bam \
  --expression_table $e \
  --chrom_sizes $s \
  --ubiquitously_expressed_genes $u \
  --gene_name_annotations symbol \
  --primary_gene_identifier symbol \
  --cellType $snpdx \
  --outdir ./
  """

}

process make_qnorm_ref {
  label 'qnorm'
  publishDir "${params.outdir}/qnorm_reference/"
  conda params.abc_code + "/abcenv.yml"

  when:
  params.mode =~ /(all)/

  input:
  tuple sn, path(enh) from MAKE_QN

  output:
  tuple sn, path(outname) into QN_REF_FILES

  script:
  outname = "${sn}.qnormReference.txt"
  """
  python ${params.abc_code}/src/makeQnormReference.py \
  --enhancers $enh \
  --outfile $outname \
  --cols ATAC.RPM
  """

}

QN_REF_FILES
  .combine(TO_NORM)
  .filter{ it[0] != it[3] } // filter out same QN_ref - sample pairs
  .set{ CATAC }

process atac_activity {
  label 'ccres'
  publishDir "${params.outdir}/neighborhoods/${snpdx}/"
  conda params.abc_code + "/abcenv.yml"

  memory { 2.GB * task.attempt }
  time { 1.hour * task.attempt }
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 5

  when:
  params.mode =~ /(all)/

  input:
  tuple qns, path(qnf), path(ccres), sn, path(bam), path(idx), path(g), path(e), path(s), path(sb), path(u) from CATAC

  output:
  tuple snpdx, path("EnhancerList.txt"), path("GeneList.txt") into PDX

  script:
  snpdx = "${sn}.${qns}QN"
  """
  python ${params.abc_code}/src/run.neighborhoods.py \
  --candidate_enhancer_regions $ccres \
  --genes $g \
  --ATAC $bam \
  --qnorm $qnf \
  --expression_table $e \
  --chrom_sizes $s \
  --ubiquitously_expressed_genes $u \
  --gene_name_annotations symbol \
  --primary_gene_identifier symbol \
  --cellType $snpdx \
  --outdir ./
  """

}

process atac_activity_no_qn {
  label 'ccres'
  publishDir "${params.outdir}/neighborhoods/${snpdx}/"
  conda params.abc_code + "/abcenv.yml"

  memory { 2.GB * task.attempt }
  time { 1.hour * task.attempt }
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 5

  when:
  params.mode =~ /(all)/

  input:
  tuple path(ccres), sn, path(bam), path(idx), path(g), path(e), path(s), path(sb), path(u) from CATAC_NO_QN

  output:
  tuple snpdx, path("EnhancerList.txt"), path("GeneList.txt") into PDX_NO_QN_FRESH

  script:
  snpdx = sn
  """
  python ${params.abc_code}/src/run.neighborhoods.py \
  --candidate_enhancer_regions $ccres \
  --genes $g \
  --ATAC $bam \
  --expression_table $e \
  --chrom_sizes $s \
  --ubiquitously_expressed_genes $u \
  --gene_name_annotations symbol \
  --primary_gene_identifier symbol \
  --cellType $snpdx \
  --outdir ./
  """

}

( PDX_NO_QN_NEW,
  PDX_NO_QN_OLD ) = ( params.pred_glob==""
                    ? [ PDX_NO_QN_FRESH,
                        Channel.empty() ]

                    : [ Channel.empty(),
                        Channel.fromPath(params.pred_glob)
                        .map{ it -> [ it.getParent().getName(), it ] }
                        .groupTuple()
                        .map{ it -> [ it[0], it[1][0], it[1][1] ] }
                      ]
                    )

PDX_NO_QN_NEW
  .mix(PDX_NO_QN_OLD)
  .set{ PDX_NO_QN }

// sample origin dependent on if HiChIP reps merged
def getSampleOrigin(x, hqn) {
  if (hqn.contains('rep')) {
    res = x.split("\\.")[0]
  } else {
    res = x.split("\\.")[0].split("_")[0]
  }
  return res
}

def getQnormOrigin(x, hqn) {
  if (hqn.contains('rep')) {
    res = x.split("\\.")[1] - ~/(QN)?$/
  } else {
    res = x.split("\\.")[1].split("_")[0]
  }
  return res
}

HICPRO_MATRIX = ( params.hicpro_matrix_dir==""
                  ? Channel.empty()
                  : Channel.fromPath( params.hicpro_matrix_dir )
                        .map{ it -> [ it.getName() - ~/(_5000)(_norm\.txt)?(\.matrix)?(_abs\.bed)?$/, it ] }
                  			.groupTuple()
                  			.dump(tag: 'bedpe_input')
                )

process make_bedpe {
  label 'bedpe'
  publishDir "${params.outdir}/${ctype}/"

  input:
  tuple sorig, path(matbed) from HICPRO_MATRIX

  output:
  tuple sorig, path(outname) into BEDPE

  script:
  if (params.hichip) {
    ctype = 'HiChIP'
  } else {
    ctype = 'HiC'
  }
  matbedsort = [ matbed[0].getName(), matbed[1].getName() ].sort()
  if (matbedsort[1].contains('_norm')) {
    mat = matbedsort[1]
    bed = matbedsort[0]
  } else {
    mat = matbedsort[0]
    bed = matbedsort[1]
  }
  outbase = mat.split("\\.")[0]
  outname = "${outbase}.${ctype}.bintotals.bedpe.gz"
  """
  Rscript ${params.bin}/hicpro_matrix_to_bedpe.R \
  --bin_totals \
  --matrix_dir . \
  --matrix $mat \
  --bed $bed \
  --name $outname \
  --outdir .
  """

}

if (params.hicpro_matrix_dir=="") {

  HIC = ( params.hic_dir==""
          ? Channel.empty()
          : Channel.fromPath(params.hic_dir)
              .map{ it -> ['HiC', it] }
        )

  PWRLAW = ( params.power_law
             ? Channel.of( ['powerlaw', ""] )
             : Channel.empty()
            )

  PDX
    .mix( PDXQN )
    .mix( PDX_NO_QN )
    .combine( Channel.fromPath(params.chrsize) )
    .combine( HIC.mix(PWRLAW) )
    .set{ PRED }

  process hichip_predict {
    label 'predict'
    publishDir "${params.outdir}/predictions/${ctype}/"
    conda params.abc_code + "/abcenv.yml"

    when:
    params.mode =~ /(all)/

    input:
    tuple sn, path(enh), path(gene), path(s), ctype, c from PRED

    output:
    tuple ctype, sn, path("${sn}/${sn}.EnhancerPredictionsAllPutative.txt.gz") into ANAL
    path("${sn}/*NonExpressedGenes.txt.gz")
    path("${sn}/*.bedpe")
    path("${sn}/*.txt")

    script:
    if (ctype=="HiC") {
      if (params.hichip) {
        ctype = ctype + 'hIP'
        """
        python ${params.abc_code}/src/predict.py \
        --enhancers $enh \
        --genes $gene \
        --HiCdir $c \
        --hic_type bedpe \
        --hichip \
        --chrom_sizes $s \
        --hic_resolution 5000 \
        --window ${params.window} \
        --score_column ABC.Score \
        --threshold .02 \
        --cellType $sn \
        --outdir ./ \
        --make_all_putative \
        --include_nonexpressed; \
        mkdir $sn; \
        for f in *txt* *bedpe* ; \
        do mv \$f ${sn}/${sn}.\$f ; \
        done
        """
      } else {
        """
        python ${params.abc_code}/src/predict.py \
        --enhancers $enh \
        --genes $gene \
        --HiCdir $c \
        --hic_type ${params.hic_type} \
        --chrom_sizes $s \
        --hic_resolution 5000 \
        --window ${params.window} \
        --scale_hic_using_powerlaw \
        --score_column ABC.Score \
        --threshold .02 \
        --cellType $sn \
        --outdir ./ \
        --make_all_putative \
        --include_nonexpressed; \
        mkdir $sn; \
        for f in *txt* *bedpe* ; \
        do mv \$f ${sn}/${sn}.\$f ; \
        done
        """
      }
    } else {
      """
      python ${params.abc_code}/src/predict.py \
      --enhancers $enh \
      --genes $gene \
      --chrom_sizes $s \
      --hic_resolution 5000 \
      --window ${params.window} \
      --scale_hic_using_powerlaw \
      --score_column powerlaw.Score \
      --threshold .02 \
      --cellType $sn \
      --outdir ./ \
      --make_all_putative \
      --include_nonexpressed; \
      mkdir $sn; \
      for f in *txt* *bedpe* ; \
      do mv \$f ${sn}/${sn}.\$f ; \
      done
      """
    }


  }

} else {

  if (params.qnorm_refs=="") {

    PDX_NO_QN
      .map{ it -> [ getSampleOrigin(it[0],params.qnorm_refs_hic[0]), it[0], it[1], it[2] ] }
      .combine( BEDPE, by: [0] )
      // sorig, sn, enh, gene, bedpe
      .combine( Channel.fromPath(params.chrsize) )
      .set{ PRED_NO_QN }

      process hichip_predict_no_qn {
        label 'predict'
        publishDir "${params.outdir}/predictions/${ctype}/"
        conda params.abc_code + "/abcenv.yml"

        when:
        params.mode =~ /(all)/

        input:
        tuple sorig, sn, path(enh), path(gene), path(bpe), path(s) from PRED_NO_QN

        output:
        tuple ctype, snanal, path("${snanal}/${snanal}.inclNonExpr.EnhancerPredictionsAllPutative.txt.gz") into ANAL_NO_QN
        path("${snanal}/*NonExpressedGenes.txt.gz")
        path("${snanal}/*.bedpe")
        path("${snanal}/*.txt")

        script:
        qnsuff = 'HiC_no_QN'
        if (!params.qnorm_refs_hic[0].contains('rep')) {
          qnsuff = 'CombRep' + qnsuff
        }
        if (bpe.getName().contains('_norm')) {
          qnsuff = 'norm' + qnsuff
        }
        snanal = "${sn}.${qnsuff}" // must be <pop>_<rep>.qnsuff for downstream processing
        bpe_uc = bpe.toString() - ~/(\.gz)$/
        if (params.hichip) {
          ctype = "HiChIP"
        } else {
          ctype = "HiC"
        }
        if (ctype=="HiC") {
          """
          gunzip -c $bpe \
          > $bpe_uc; \
          for chr in `cut -f 1 $bpe_uc | sort | uniq`; \
          do echo \$chr; \
          mkdir -p \$chr; \
          grep -w \$chr $bpe_uc \
          | gzip -c \
          > \${chr}/\${chr}.bedpe.gz; \
          done; \
          rm $bpe_uc; \
          python ${params.abc_code}/src/predict.py \
          --enhancers $enh \
          --genes $gene \
          --HiCdir . \
          --hic_type ${params.hic_type} \
          --chrom_sizes $s \
          --hic_resolution 5000 \
          --window ${params.window} \
          --scale_hic_using_powerlaw \
          --score_column ABC.Score \
          --threshold .015 \
          --cellType $snanal \
          --outdir ./ \
          --make_all_putative \
          --include_nonexpressed; \
          mkdir $snanal; \
          for f in *txt* *bedpe* ; \
          do mv \$f ${snanal}/${snanal}.inclNonExpr.\$f ; \
          done
          """
        } else {
          """
          gunzip -c $bpe \
          > $bpe_uc; \
          for chr in `cut -f 1 $bpe_uc | sort | uniq`; \
          do echo \$chr; \
          mkdir -p \$chr; \
          grep -w \$chr $bpe_uc \
          | gzip -c \
          > \${chr}/\${chr}.bedpe.gz; \
          done; \
          rm $bpe_uc; \
          python ${params.abc_code}/src/predict.py \
          --enhancers $enh \
          --genes $gene \
          --HiCdir . \
          --hic_type bedpe \
          --hichip \
          --chrom_sizes $s \
          --hic_resolution 5000 \
          --window ${params.window} \
          --score_column ABC.Score \
          --threshold .015 \
          --cellType $snanal \
          --outdir ./ \
          --make_all_putative \
          --include_nonexpressed; \
          mkdir $snanal; \
          for f in *txt* *bedpe* ; \
          do mv \$f ${snanal}/${snanal}.inclNonExpr.\$f ; \
          done
          """
        }

      }

      if (params.mean_qn) {

        ANAL_NO_QN
          .map{ it -> [ it[0], it[1].split("\\.")[1], it[2] ] }
          .groupTuple(by: [0,1])
          .set{ ANAL_TO_QN }

        process hichip_mean_qn {
          label 'mean_qnorm_hic'
          publishDir "${params.outdir}/predictions/${ctype}/meanQN/", pattern: "*.txt.gz"
          publishDir "${params.outdir}/predictions/${ctype}/meanQN/plots/", pattern: "*.png"

          when:
          params.mode =~ /(all)/

          input:
          tuple ctype, qnsuff, path(preds) from ANAL_TO_QN

          output:
          tuple ctype, qnsuffsplit, path("*.meanQN.decomp.allZerosFilt.EnhancerPredictionsAllPutative.txt.gz") into ANAL_GROUPED
          path("ATAC_RPM_QN.boxplot.png")
          path("ATAC_RPM.boxplot.png")

          script:
          qnsuffsplit = ".meanQN"
          """
          Rscript ${params.bin}/meanQN_ABCscore.R \
          --filt_zeros all \
          --decomp_scores \
          --plotdir ${params.plotdir} \
          --pred_suffix $qnsuff \
          --pred_dir . \
          --outdir .
          """

        }

        ANAL_GROUPED
          // flatten list of prediciton files into separate channels
          // will be reversed before differential_abc
          .transpose(by: 2) // reverse of groupTuple(by: [0,1])
          .set{ ANAL }

      } else {

        ANAL_NO_QN.set{ ANAL }

      }

  } else {

    BEDPE.into{ BEDPE_QN; BEDPE_TO_NORM }

    PDX
      .mix( PDXQN )
      // ${sn}.${sn}QN, path("EnhancerList.txt"), path("GeneList.txt")
      .tap{ PRED_TO_NORM_UNFILT }
      .combine( BEDPE_QN )
      .dump(tag: "refPreFilt")
      // ${sn}.${sn}QN, enh, gene, pop_rep, bedpe
      // OR
      // ${sn}.${sn}QN, enh, gene, pop, bedpe
      // use same population/sample origin for HiC/HiChIP QN as was used for ATAC QN
      // if using separate HiChIP reps to preserve ABC score independence and enable
      // reproducibility analysis, samples will be matched arbitrarily by replicate name
      .filter{ (getSampleOrigin(it[0],params.qnorm_refs_hic[0]) == it[3]) \
            && (getQnormOrigin(it[0],params.qnorm_refs_hic[0]) == it[3]) }
      .dump(tag: "refPostFilt")
      .combine( Channel.fromPath(params.chrsize) )
      .set{ PRED_QN_REFS }

    process hichip_predict_qnorm_ref {
      label 'predict'
      publishDir "${params.outdir}/predictions/${ctype}/"
      conda params.abc_code + "/abcenv.yml"

      when:
      params.mode =~ /(all)/

      input:
      tuple sn, path(enh), path(gene), sorig, path(bpe), path(s) from PRED_QN_REFS

      output:
      tuple ctype, snanal, path("${snanal}/${snanal}.EnhancerPredictionsAllPutative.txt.gz") into ANAL_REF
      tuple sn, path("${snanal}/${snanal}.EnhancerPredictionsAllPutative.txt.gz") into MAKE_QN_HIC_UNFILT
      path("${snanal}/*NonExpressedGenes.txt.gz")
      path("${snanal}/*.bedpe")
      path("${snanal}/*.txt")

      script:
      if (bpe.getName().contains('_norm')) {
        snanal = "${sn}_${sorig}_normHiCQN"
      } else {
        snanal = "${sn}_${sorig}_HiCQN"
      }
      bpe_uc = bpe.toString() - ~/(\.gz)$/
      if (params.hichip) {
        ctype = "HiChIP"
      } else {
        ctype = "HiC"
      }
      if (ctype=="HiC") {
        """
        gunzip -c $bpe \
        > $bpe_uc; \
        for chr in `cut -f 1 $bpe_uc | sort | uniq`; \
        do echo \$chr; \
        mkdir -p \$chr; \
        grep -w \$chr $bpe_uc \
        | gzip -c \
        > \${chr}/\${chr}.bedpe.gz; \
        done; \
        rm $bpe_uc; \
        python ${params.abc_code}/src/predict.py \
        --enhancers $enh \
        --genes $gene \
        --HiCdir . \
        --hic_type ${params.hic_type} \
        --chrom_sizes $s \
        --hic_resolution 5000 \
        --window ${params.window} \
        --scale_hic_using_powerlaw \
        --score_column ABC.Score \
        --threshold .015 \
        --cellType $snanal \
        --outdir ./ \
        --make_all_putative \
        --include_nonexpressed; \
        mkdir $snanal; \
        for f in *txt* *bedpe* ; \
        do mv \$f ${snanal}/${snanal}.\$f ; \
        done
        """
      } else {
        """
        gunzip -c $bpe \
        > $bpe_uc; \
        for chr in `cut -f 1 $bpe_uc | sort | uniq`; \
        do echo \$chr; \
        mkdir -p \$chr; \
        grep -w \$chr $bpe_uc \
        | gzip -c \
        > \${chr}/\${chr}.bedpe.gz; \
        done; \
        rm $bpe_uc; \
        python ${params.abc_code}/src/predict.py \
        --enhancers $enh \
        --genes $gene \
        --HiCdir . \
        --hic_type bedpe \
        --hichip \
        --chrom_sizes $s \
        --hic_resolution 5000 \
        --window ${params.window} \
        --score_column ABC.Score \
        --threshold .015 \
        --cellType $snanal \
        --outdir ./ \
        --make_all_putative \
        --include_nonexpressed; \
        mkdir $snanal; \
        for f in *txt* *bedpe* ; \
        do mv \$f ${snanal}/${snanal}.\$f ; \
        done
        """
      }

    }

    // filter to avoid making duplicate qnorm refs
    MAKE_QN_HIC_UNFILT
      .filter{ it[0].split("\\.")[0] in params.qnorm_refs }
      .set{ MAKE_QN_HIC }

    process make_qnorm_ref_hic {
      label 'qnorm_hic'
      publishDir "${params.outdir}/qnorm_reference/"
      conda params.abc_code + "/abcenv.yml"

      when:
      params.mode =~ /(all)/

      input:
      tuple sn, path(preds) from MAKE_QN_HIC

      output:
      tuple sorig, path(outname) into QN_REF_FILES_HIC

      script:
      sorig = getSampleOrigin(sn, params.qnorm_refs_hic[0])
      if (preds.getName().contains('_norm')) {
        outname = "${sorig}_norm.HiC.qnormReference.txt"
      } else {
        outname = "${sorig}.HiC.qnormReference.txt"
      }
      """
      python ${params.abc_code}/src/makeQnormReference.py \
      --enhancers $preds \
      --outfile $outname \
      --cols hic_contact
      """

    }

    PRED_TO_NORM_UNFILT
      // sorig, sn, path(enh), path(gene)
      .map{ it -> [ getSampleOrigin(it[0], params.qnorm_refs_hic[0]), it[0], it[1], it[2] ] }
      // sorig, path(hqn_ref)
      .combine( QN_REF_FILES_HIC )
      .dump(tag: "QNpreFilt")
      // only qnorm HiC/HiChIP if predictions are for non-reference population/sample origin
      .filter{ (it[0] != it[4]) && (getQnormOrigin(it[1], params.qnorm_refs_hic[0]) == it[4]) }
      .dump(tag: "QNpostFilt")
      // sorig, path(bpe)
      .combine( BEDPE_TO_NORM, by: [0] ) // should filter out bedpe sample origins not part of analysis
      .dump(tag: "QNbpeComb")
      .combine( Channel.fromPath(params.chrsize) )
      // sorig, sn, path(enh), path(gene), sorighqn, path(hqn_ref), path(bpe), path(s)
      .set{ PRED_TO_NORM }


    process hichip_predict_qnorm {
      label 'predict'
      publishDir "${params.outdir}/predictions/${ctype}/"
      conda params.abc_code + "/abcenv.yml"

      when:
      params.mode =~ /(all)/

      input:
      tuple sorig, sn, path(enh), path(gene), sorighqn, path(hqn_ref), path(bpe), path(s) from PRED_TO_NORM

      output:
      tuple ctype, snanal, path("${snanal}/${snanal}.EnhancerPredictionsAllPutative.txt.gz") into ANAL_NORMED
      path("${snanal}/*NonExpressedGenes.txt.gz")
      path("${snanal}/*.bedpe")
      path("${snanal}/*.txt")

      script:
      if (bpe.getName().contains('_norm')) {
        snanal = "${sn}_${sorighqn}_normHiCQN"
      } else {
        snanal = "${sn}_${sorighqn}_HiCQN"
      }
      bpe_uc = bpe.toString() - ~/(\.gz)$/
      if (params.hichip) {
        ctype = "HiChIP"
      } else {
        ctype = "HiC"
      }
      if (ctype=="HiC") {
        """
        gunzip -c $bpe \
        > $bpe_uc; \
        for chr in `cut -f 1 $bpe_uc | sort | uniq`; \
        do echo \$chr; \
        mkdir -p \$chr; \
        grep -w \$chr $bpe_uc \
        | gzip -c \
        > \${chr}/\${chr}.bedpe.gz; \
        done; \
        rm $bpe_uc; \
        python ${params.abc_code}/src/predict.py \
        --enhancers $enh \
        --genes $gene \
        --HiCdir . \
        --hic_type ${params.hic_type} \
        --qnorm $hqn_ref \
        --chrom_sizes $s \
        --hic_resolution 5000 \
        --window ${params.window} \
        --scale_hic_using_powerlaw \
        --score_column ABC.Score \
        --threshold .015 \
        --cellType $snanal \
        --outdir ./ \
        --make_all_putative \
        --include_nonexpressed; \
        mkdir $snanal; \
        for f in *txt* *bedpe* ; \
        do mv \$f ${snanal}/${snanal}.\$f ; \
        done
        """
      } else {
        """
        gunzip -c $bpe \
        > $bpe_uc; \
        for chr in `cut -f 1 $bpe_uc | sort | uniq`; \
        do echo \$chr; \
        mkdir -p \$chr; \
        grep -w \$chr $bpe_uc \
        | gzip -c \
        > \${chr}/\${chr}.bedpe.gz; \
        done; \
        rm $bpe_uc; \
        python ${params.abc_code}/src/predict.py \
        --enhancers $enh \
        --genes $gene \
        --HiCdir . \
        --hic_type bedpe \
        --hichip \
        --qnorm $hqn_ref \
        --chrom_sizes $s \
        --hic_resolution 5000 \
        --window ${params.window} \
        --score_column ABC.Score \
        --threshold .015 \
        --cellType $snanal \
        --outdir ./ \
        --make_all_putative \
        --include_nonexpressed; \
        mkdir $snanal; \
        for f in *txt* *bedpe* ; \
        do mv \$f ${snanal}/${snanal}.\$f ; \
        done
        """
      }

    }

    ANAL_REF
      .mix(ANAL_NORMED)
      .set{ ANAL }

  }

}

// process output of predict by qnorm reference
// ctype, sn, path("${sn}/EnhancerPredictionsAllPutative.txt.gz")
ANAL
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
  // path("${outbase}.${score_col}.scoresAllSamps.txt.gz") //
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

( OVERLAP_PRE_NEW,
  OVERLAP_PRE_OLD ) = ( params.diff_scores==""
                    ? [ OVERLAP_PRE_FRESH,
                        Channel.empty() ]

                    : [ Channel.empty(),
                        Channel.fromList(params.diff_scores)
                        .map{ it -> file(it) }
                        .map{ it -> [ it.getParent().getName(),
                                      it.getName().split("\\.")[1],
                                      "${(it.getName() - ~/(\.txt\.gz)?$/).split("\\.")[-2]}.${(it.getName() - ~/(\.txt\.gz)?$/).split("\\.")[-1]}",
                                      it ] }
                      ]
                    )

OVERLAP_PRE_NEW
  .mix(OVERLAP_PRE_OLD)
  .set{ OVERLAP_PRE }

// combine DE datasets for overlap
OVERLAP_PRE
  .tap{ CONCAT_PRE }
  .combine( Channel.fromList( params.de_datasets ) )
  .combine( Channel.fromPath( params.ihs ) )
  .set{ OVERLAP }

process overlap {
  label 'diff'
  publishDir "${params.outdir}/diffABC/${ctype}/overlap/allZerosFilt", pattern: "*.txt*"
  publishDir "${params.outdir}/diffABC/${ctype}/overlap/allZerosFilt/plots", pattern: "*.png"

  when:
  params.mode =~ /(all|overlap)/

  input:
  tuple ctype, qn, score_col, path(dabc), path(de), path(ihs) from OVERLAP

  output:
  tuple ctype, qn, score_col, path("*fisherEnrichments*.txt*") into AGG_ANAL_PRE
  path("*.txt") // raw overlap of E-G pairs with DE
  path("*.enhReinforcingDirection.gprofiler*.txt.gz")
  path("*.enhReinforcingDirection.allPairs.txt.gz")
  path("*.fisherEnrichments.txt")
  path("*.fisherEnrichmentsDirection.txt")
  path("*.fisherEnrichmentsDirectionTopN.txt")
  path("*.fisherEnrichments.png")
  path("*.fisherEnrichmentsDirection.png")
  path("*.enhReinforcingDirection.fisherEnrichments.png")
  path("*.effectSizeCorrelations.png")

  script:
  debase = de.getName() - ~/(\.txt)?(\.gz)?$/
  outbase = "allZerosFilt.${score_col}.${qn}.${params.min_nonzero_hic}.AFR_EUR.diff_overlap_${debase}"
  if (debase.contains('LeaEtAl2021')) {
    overlap_opts = "--pval_prefix pval --ES_prefix beta --enrichvar condition --enrichvarlab Condition"
  } else if (debase.contains('RandolphEtAl2021_scDecomp')) {
    overlap_opts = "--pval_prefix p --facetscales free --facetvar datatype --enrichgroup condition --enrichvar celltype --enrichvarlab Cell_type --ES_prefix beta"
  } else {
    overlap_opts = "--pval_prefix lfsr --ES_prefix beta --enrichgroup condition --enrichvar celltype --enrichvarlab Cell_type"
  }
  """
  Rscript ${params.bin}/overlap_diffABC_DE.R \
  $overlap_opts \
  --dist_thresh 10000 \
  --plotdir ${params.plotdir} \
  --DE $de \
  --diffABC_preds $dabc \
  --score_column $score_col \
  --iHS $ihs \
  --top_diffABC \
  --name $outbase \
  --outdir .
  """

}

AGG_ANAL_PRE
  // remove score_col then aggregate
  .map{ it -> [ it[0], it[1], it[3] ] }
  .groupTuple(by: [0,1] )
  // flatten grouped file lists to one path type list of files
  .map{ it -> [ it[0], it[1], it[2].flatten() ] }
  .set{ AGG_ANAL }

// for facetted MA plot of all scores
CONCAT_PRE
  // remove score_col then aggregate
  .map{ it -> [ it[0], it[1], it[3] ] }
  .groupTuple(by: [0,1] )
  .combine( AGG_ANAL, by: [0,1] )
  .combine( Channel.fromPath( params.ihs ) )
  .combine( Channel.fromPath( params.fst ) )
  .dump(tag: 'final')
  .set{ CONCAT }

process concat_scores {
  label 'diff' // shouldn't require tons of memory
  publishDir "${params.outdir}/diffABC/${ctype}/aggregate", pattern: "*.txt*"
  publishDir "${params.outdir}/diffABC/${ctype}/aggregate/plots", pattern: "*.png"

  when:
  params.mode =~ /(all|overlap)/

  input:
  tuple ctype, qn, path(dscores), path(fes), path(ihs), path(fst) from CONCAT

  output:
  tuple ctype, qn, path("${outbase}.txt.gz") into SPOT_CHECK, DE_MODEL_PRE
  path("*.fgseaHallmarkGOBP.txt.gz")
  path("*.selectionCandidatesBackground.txt.gz")
  path("*.DEoverlap.fisherEnrichments.txt")
  path("*.DEoverlap.fisherEnrichmentsDirection.txt")
  path("*.DEoverlap.fisherEnrichmentsDirectionTopN.txt")
  path("*.enhReinforcingDirection.fisherEnrichments.txt")
  path("*.selectionSNPs.fisherEnrichments.txt")
  path("*.fgseaTopPathways.png")
  path("*.componentScoreMAplot.png")
  path("*.fisherEnrichments.DEgene.LeaRandolph.png")
  path("*.fisherEnrichments.DEgene.scDecomp.png")
  path("*.fisherEnrichments.DEdirection.LeaRandolph.png")
  path("*.fisherEnrichments.DEdirection.scDecomp.png")
  path("*.fisherEnrichments.DEdirection.LeaRandolph.topDiffScore.png")
  path("*.fisherEnrichments.DEdirection.scDecomp.topDiffScore.png")
  path("*.fisherEnrichments.enhReinforcingDirection.overall.png")
  path("*.fisherEnrichments.enhReinforcingDirection.LeaRandolphDE.png")
  path("*.fisherEnrichments.enhReinforcingDirection.scDecompDE.png")
  path("*.selectionSNPs.fisherEnrichments.png")
  path("*.iHS.Wilcoxon.png")
  path("*.FST.Wilcoxon.png")

  script:
  outpre = dscores[0].getName() - ~/(\.hic)?(\.chip)?(\.atac)?(\.ABC)?(\.Score)?(\.txt)?(\.gz)?$/
  outbase = outpre + ".allComponents"
  """
  Rscript ${params.bin}/concatDecompABC.R \
  --dist_thresh 10000 \
  --plotdir ${params.plotdir} \
  --diff_pattern .Score.txt.gz \
  --reinf_pattern .enhReinforcingDirection.fisherEnrichments.txt.gz \
  --agg_dir . \
  --iHS $ihs \
  --FST $fst \
  --name $outbase \
  --outdir .
  """

}

Channel.fromList( params.de_model_datasets )
    .combine( Channel.fromList( params.de_model_conditions ), by: [0] )
    .transpose(by: 2)
    .combine( Channel.fromList( params.de_model_celltypes ), by: [0] )
    .transpose(by: 3)
    .combine( DE_MODEL_PRE )
    .combine( Channel.fromList( params.de_model_pred ) )
    .combine( Channel.fromList( params.de_model_response ) )
    .set{ DE_MODEL }


process de_model {
  label 'model' // shouldn't require tons of memory
  publishDir "${params.outdir}/diffABC/${ctype}/model_de", pattern: "*.txt"
  publishDir "${params.outdir}/diffABC/${ctype}/model_de/plots", pattern: "*.png"

  when:
  params.mode =~ /(model)/

  input:
  tuple auth, path(de), cdt, ct, ctype, qn, path(dscores), p, r from DE_MODEL

  output:
  tuple path("${outbase}*.txt") into AGG_MODELS
  path("*.png")

  script:
  outpre = dscores.getName() - ~/(\.activity)?(\.allComponents)?(\.txt)?(\.gz)?$/
  outbase = "${outpre}.${auth}_${ct}_${cdt}.${p}_enhancers.${r}_response"
  if (auth.contains('LeaEtAl2021')) {
    prefix_opts = "--pval_prefix pval --ES_prefix beta"
  } else {
    prefix_opts = "--pval_prefix lfsr --ES_prefix beta"
  }
  if (p=='all') {
    pred_opts = ""
  } else {
    pred_opts = "--top_enh ${p}"
  }
  """
  Rscript ${params.bin}/model_DE.R \
  --diff_all $dscores \
  --DE $de \
  --plotdir ${params.plotdir} \
  ${prefix_opts} \
  --condition $cdt \
  --celltype $ct \
  --response_var $r \
  ${pred_opts} \
  --name $outbase
  """

}
