Nextflow script to run 100x pySCENIC (modified version of https://github.com/aertslab/scenic-nf/blob/master/main.nf)

## Baselines

### Motif-based SCENIC w/o dropout masking

Nextflow script (main.nf):
```
#!/usr/bin/env nextflow


println( "\n***\nParameters in use:")
params.each { println "${it}" }

// channel for SCENIC databases resources:
motifDbs = Channel
    .fromPath( params.motif_dbs )
    .collect() // use all files together in the ctx command

n = Channel.fromPath(params.motif_dbs).count().get()
if( n==1 ) {
    println( "***\nWARNING: only using a single feather database:\n  ${motifDbs.get()[0]}.\nTo include all database files using pattern matching, make sure the value for the '--db' parameter is enclosed in quotes!\n***\n" )
} else {
    println( "***\nUsing $n feather databases:")
    motifDbs.get().each {
        println "  ${it}"
    }
    println( "***\n")
}

expr = file(params.expr)
tfs = file(params.TFs)
motifs = file(params.motif_tf_annotation)
// tracks = file(params.track_tf_annotation)
nbRuns = params.nb_runs

ctxMaskDropouts = params.ctx_mask_dropouts
ctxMaskDropoutsTag = ''
if(ctxMaskDropouts == 'no') {
    ctxMaskDropoutsTag = '_nodom'
}

// Resources limits
_maxForks = params.max_forks
maxCpus = params.global.threads

// UTILS
def runName = { it.getName().split('__')[0] }

if(params.parallel_framework == 'dask') {
    process dask_GRNinference {
        cache 'deep'

        clusterOptions "-l nodes=1:ppn=${params.global.threads} -l pmem=2gb -l walltime=${params.qsub_walltime_hours}:00:00 -A ${params.global.qsubaccount}"

        maxForks _maxForks
        cpus maxCpus

        input:
        each runId from 1..nbRuns
        file TFs from tfs
        file exprMat from expr

        output:
        file "run_${runId}__adj.tsv" into grn, grnSave

        """
        pyscenic grn \
            --num_workers ${params.global.threads} \
            -o "run_${runId}__adj.tsv" \
            --method ${params.grn} \
            --cell_id_attribute ${params.cell_id_attribute} \
            --gene_attribute ${params.gene_attribute} \
            ${exprMat} \
            ${TFs}
        """
    }
}

if(params.parallel_framework == 'multiprocessing_pool') {
    process multiprocessingPool_GRNinference {
        cache 'deep'

        clusterOptions "-l nodes=1:ppn=${params.global.threads} -l pmem=2gb -l walltime=${params.qsub_walltime_hours}:00:00 -A ${params.global.qsubaccount}"

        maxForks _maxForks
        cpus maxCpus

        input:
        each runId from 1..nbRuns
        file TFs from tfs
        file exprMat from expr

        output:
        file "run_${runId}__adj.tsv" into grn, grnSave

        """
        python /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/GitHub/scenic-nf/bin/grnboost2_without_dask.py \
            --output "run_${runId}__adj.tsv" \
            --num_workers ${params.global.threads} \
            --seed ${runId} \
            --cell_id_attribute ${params.cell_id_attribute} \
            --gene_attribute ${params.gene_attribute} \
            ${exprMat} \
            ${TFs}
        """
    }
}

process motif_cisTarget {

    maxForks _maxForks
    cpus maxCpus

    cache 'deep'

    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l pmem=2gb -l walltime=${params.qsub_walltime_hours}:00:00 -A ${params.global.qsubaccount}"

    input:
    file exprMat from expr
    file adj from grn
    file motifDb from motifDbs
    file motif from motifs

    output:
    file "${runName(adj)}__motif_reg${ctxMaskDropoutsTag}.csv" into motifRegulons, motifRegulonsSave

    """
    pyscenic ctx \
        ${adj} \
        ${motifDb} \
        --annotations_fname ${motif} \
        --expression_mtx_fname ${exprMat} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --mode "dask_multiprocessing" \
        --[mask_dropouts] "no" \
        --output "${runName(adj)}__motif_reg${ctxMaskDropoutsTag}.csv" \
        --num_workers ${params.global.threads} \
    """
}

process motif_AUCell {

    maxForks _maxForks
    cpus maxCpus

    cache 'deep'

    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l pmem=1gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
    file exprMat from expr
    file reg from motifRegulons

    output:
    file "${runName(reg)}__motif${ctxMaskDropoutsTag}_${params.output}" into motifAUCMatrix

    """
    pyscenic aucell \
        $exprMat \
        $reg \
        -o "${runName(reg)}__motif${ctxMaskDropoutsTag}_${params.output}" \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --num_workers ${params.global.threads}
    """
}

def save = {
    (full, run, filename, ext) = (it.getName() =~ /(.+)__(.+)\.(.+)/)[0]
    if( params.nb_runs==1 ) {
        outDir = file( params.global.outdir )
    } else if( params.nb_runs>1 ) {
        outDir = file( params.global.outdir+"/$run" )
    }
    result = outDir.mkdirs()
    println result ? "$run finished." : "Cannot create directory: $outDir"   
    Channel
        .fromPath(it)
        .collectFile(name: "${filename}.${ext}", storeDir: outDir)
}

grnSave.subscribe { 
    save(it) 
}
motifRegulonsSave.subscribe { 
    save(it) 
}
motifAUCMatrix.subscribe { 
    save(it) 
}
```

Run command
```
~/nextflow/19.07.0.5106/bin/nextflow run ~/main.nf \
    -profile singularity,qsub \
    --expr 10x_MEL_baselines_filtered_log_CPM.loom \
    --output pyscenic_output.loom \
    --TFs allTFs_hg38.txt \
    --threads 20 \
    --motif_dbs "cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-*.feather" \
    --motif_tf_annotation cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --nb_runs 100 \
    --parallel_framework multiprocessing_pool \
    --ctx_mask_dropouts no \
    --max_forks 10 \
    --qsubaccount lp_symbiosys \
    --qsub_walltime_hours 24 \
    -with-report report.html \
    -with-trace \
    -resume
```

### Track-based SCENIC w/o dropout masking

Nextflow script (main_tracks_wo_grn.nf):
```
#!/usr/bin/env nextflow


println( "\n***\nParameters in use:")
params.each { println "${it}" }

// channel for SCENIC databases resources:
trackDbs = Channel
    .fromPath( params.track_dbs )
    .collect() // use all files together in the ctx command

n = Channel.fromPath(params.track_dbs).count().get()
if( n==1 ) {
    println( "***\nWARNING: only using a single feather database:\n  ${trackDbs.get()[0]}.\nTo include all database files using pattern matching, make sure the value for the '--db' parameter is enclosed in quotes!\n***\n" )
} else {
    println( "***\nUsing $n feather databases:")
    trackDbs.get().each {
        println "  ${it}"
    }
    println( "***\n")
}

Channel
        .fromPath("${params.outdir}/**/adj.tsv")
        .set { grn }
        .println()


expr = file(params.expr)
tfs = file(params.TFs)
tracks = file(params.track_tf_annotation)
nbRuns = params.nb_runs

ctxMaskDropouts = params.ctx_mask_dropouts
ctxMaskDropoutsTag = ''
if(ctxMaskDropouts == 'no') {
    ctxMaskDropoutsTag = '_nodom'
}

// Resources limits
_maxForks = params.max_forks
maxCpus = params.threads

// UTILS
def runName = { it.toRealPath().getParent().getName() }

def runName2 = { it.getName().split('__')[0] }

process track_cisTarget {

    maxForks _maxForks
    cpus maxCpus

    cache 'deep'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l pmem=2gb -l walltime=${params.qsub_walltime_hours}:00:00 -A ${params.qsubaccount}"

    input:
    file exprMat from expr
    file adj from grn
    file trackDb from trackDbs
    file track from tracks

    output:
    file "${runName(adj)}__track_reg${ctxMaskDropoutsTag}.csv" into trackRegulons, trackRegulonsSave

    """
    pyscenic ctx \
        ${adj} \
        ${trackDb} \
        --annotations_fname ${track} \
        --expression_mtx_fname ${exprMat} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --mode "dask_multiprocessing" \
        --mask_dropouts "no" \
        --output "${runName(adj)}__track_reg${ctxMaskDropoutsTag}.csv" \
        --num_workers ${params.threads} \
    """
}

process track_AUCell {

    maxForks _maxForks
    cpus maxCpus

    cache 'deep'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l pmem=1gb -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
    file exprMat from expr
    file reg from trackRegulons

    output:
    file "${runName2(reg)}__track${ctxMaskDropoutsTag}_${params.output}" into trackAUCMatrix

    """
    pyscenic aucell \
        $exprMat \
        $reg \
        -o "${runName2(reg)}__track${ctxMaskDropoutsTag}_${params.output}" \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --num_workers ${params.threads}
    """
}

def save = {
    (full, run, filename, ext) = (it.getName() =~ /(.+)__(.+)\.(.+)/)[0]
    if( params.nb_runs==1 ) {
        outDir = file( params.outdir )
    } else if( params.nb_runs>1 ) {
        outDir = file( params.outdir+"/$run" )
    }
    result = outDir.mkdirs()
    println result ? "$run finished." : "Cannot create directory: $outDir"   
    Channel
        .fromPath(it)
        .collectFile(name: "${filename}.${ext}", storeDir: outDir)
}

trackRegulonsSave.subscribe { 
    save(it) 
}
trackAUCMatrix.subscribe { 
    save(it) 
}
```

Run command:
```
~/nextflow/19.07.0.5106/bin/nextflow run ~/main_tracks_wo_grn.nf \
    -profile singularity \
    --outdir Filtered_LogCPM_Matrix_MultiRuns/scenic \
    --expr 10x_MEL_baselines_filtered_log_CPM.loom \
    --output pyscenic_output.loom \
    --TFs allTFs_hg38.txt \
    --threads 10 \
    --track_dbs "cistarget/databases/homo_sapiens/hg19/refseq_r45/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg19-*.feather" \
    --track_tf_annotation icistarget-data/annotations/homo_sapiens/hg19/track_annotations/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg19.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv \
    --nb_runs 100 \
    --parallel_framework multiprocessing_pool \
    --ctx_mask_dropouts no \
    --max_forks 1 \
    --qsubaccount lp_symbiosys \
    --qsub_walltime_hours 24 \
    -with-report report.html \
    -with-trace
```

## M0087 SOX w/o TL

Nextflow config file:
```
manifest {
   name = 'aertslab/SingleCellTxBenchmark'
   description = 'A repository of pipelines for single-cell data in Nextflow DSL2'
   homePage = 'https://github.com/aertslab/SingleCellTxBenchmark'
   version = '0.3.1'
   mainScript = 'main.nf'
   defaultBranch = 'master'
   nextflowVersion = '!19.09.0-edge'
}

params {
   global {
      project_name = 'Human Melanoma Phenotype Switching'
      outdir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/10.M0087_SOX_wo_TL/pySCENIC_100x/out'
      tenx_folder = ''
      tracedir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/10.M0087_SOX_wo_TL/pySCENIC_100x/out/pipeline_reports'
      qsubaccount = 'lp_symbiosys'
   }
   sc {
      file_converter {
         iff = '10x_mtx'
         off = 'h5ad'
         useFilteredMatrix = true
      }
      file_annotator {
         iff = '10x_mtx'
         off = 'h5ad'
         type = 'sample'
         metaDataFilePath = ''
      }
      file_concatenator {
         join = 'outer'
         iff = '10x_mtx'
         off = 'h5ad'
      }
      star_concatenator {
         stranded = 'no'
         off = 'tsv'
      }
      scenic {
         container = 'docker://dweemx/sctx-pyscenic:0.9.19-4'
         scenicoutdir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/10.M0087_SOX_wo_TL/pySCENIC_100x/out'
         filteredloom = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/10.M0087_SOX_wo_TL/pySCENIC_100x/data/M0087_SOX_wo_TL_logCPM.loom'
         scenicOutputLoom = 'SCENIC_output.loom'
         scenicScopeOutputLoom = 'SCENIC_SCope_output.loom'
         numWorkers = '16'
         maxForks = 20
         mode = 'dask_multiprocessing'
         client_or_address = ''
         cell_id_attribute = 'CellID'
         gene_attribute = 'Gene'
         grn {
            TFs = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/tf/allTFs_hg38.txt'
            seed = ''
            pmem = '3750mb'
         }
         cistarget {
            adj = 'adj.tsv'
            mtfDB = '/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-*.feather'
            mtfANN = '/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
            trkDB = '/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg19/refseq_r45/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg19-*.feather'
            trkANN = '/ddn1/vol1/staging/leuven/stg_00002/icistarget-data/annotations/homo_sapiens/hg19/track_annotations/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg19.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv'
            type = ''
            output = 'reg.csv'
            rank_threshold = '5000'
            auc_threshold = '0.05'
            nes_threshold = '3.0'
            min_orthologous_identity = '0.0'
            max_similarity_fdr = '0.001'
            annotations_fname = ''
            thresholds = '0.75,0.90'
            top_n_targets = '50'
            top_n_regulators = '5,10,50'
            min_genes = '20'
            pmem = '7500mb'
         }
         aucell {
            output = 'aucell_output.loom'
            rank_threshold = '5000'
            auc_threshold = '0.05'
            nes_threshold = '3.0'
            min_genes_regulon = 5
            min_regulon_gene_occurrence = 5
            pmem = '7500mb'
         }
         aggregate_features {
            use_chunking = true
            output_format = 'csv'
            compression = 'gzip'
         }
         numRuns = 100
      }
      scope {
         genome = ''
         tree {
            level_1 = "Human_Melanoma"
            level_2 = "10x"
            level_3 = ""
         }
      }
   }
}

process {
   executor = 'local'
}

singularity {
   enabled = true
   autoMounts = true
   runOptions = '-B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ -B /staging/leuven/res_00001/:/staging/leuven/res_00001/ -B /ddn1/vol1/staging/leuven/res_00001/:/ddn1/vol1/staging/leuven/res_00001/'
}
```

Run command:
```
nextflow -C nextflow.config run \
   -r scenic_multi_runs_fixes \
   aertslab/SingleCellTxBenchmark \
      -entry scenic
```

## MM074 SOX

Nextflow config file,
```
manifest {
   name = 'aertslab/SingleCellTxBenchmark'
   description = 'A repository of pipelines for single-cell data in Nextflow DSL2'
   homePage = 'https://github.com/aertslab/SingleCellTxBenchmark'
   version = '0.3.1'
   mainScript = 'main.nf'
   defaultBranch = 'master'
   nextflowVersion = '!19.09.0-edge'
}

params {
   global {
      project_name = 'Human Melanoma Phenotype Switching'
      outdir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/5.MM074_SOX/pySCENIC_100x/out'
      tenx_folder = 'data/10x/1k_pbmc/1k_pbmc_*/outs/'
      tracedir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/5.MM074_SOX/pySCENIC_100x/out/pipeline_reports'
      qsubaccount = 'lp_symbiosys'
   }
   sc {
      file_converter {
         iff = '10x_mtx'
         off = 'h5ad'
         useFilteredMatrix = true
      }
      file_annotator {
         iff = '10x_mtx'
         off = 'h5ad'
         type = 'sample'
         metaDataFilePath = 'data/10x/1k_pbmc/metadata.tsv'
      }
      file_concatenator {
         join = 'outer'
         iff = '10x_mtx'
         off = 'h5ad'
      }
      star_concatenator {
         stranded = 'no'
         off = 'tsv'
      }
      scenic {
         container = 'docker://dweemx/sctx-pyscenic:0.9.19-4'
         scenicoutdir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/5.MM074_SOX/pySCENIC_100x/out/'
         filteredloom = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/5.MM074_SOX/pySCENIC_100x/data/MM074_SOX_logCPM.loom'
         scenicOutputLoom = 'SCENIC_output.loom'
         scenicScopeOutputLoom = 'SCENIC_SCope_output.loom'
         numWorkers = '16'
         maxForks = 20
         mode = 'dask_multiprocessing'
         client_or_address = ''
         cell_id_attribute = 'CellID'
         gene_attribute = 'Gene'
         grn {
            TFs = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/tf/allTFs_hg38.txt'
            seed = ''
            pmem = '3750mb'
         }
         cistarget {
            adj = 'adj.tsv'
            mtfDB = '/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-*.feather'
            mtfANN = '/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
            trkDB = '/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg19/refseq_r45/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg19-*.feather'
            trkANN = '/ddn1/vol1/staging/leuven/stg_00002/icistarget-data/annotations/homo_sapiens/hg19/track_annotations/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg19.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv'
            type = ''
            output = 'reg.csv'
            rank_threshold = '5000'
            auc_threshold = '0.05'
            nes_threshold = '3.0'
            min_orthologous_identity = '0.0'
            max_similarity_fdr = '0.001'
            annotations_fname = ''
            thresholds = '0.75,0.90'
            top_n_targets = '50'
            top_n_regulators = '5,10,50'
            min_genes = '20'
            pmem = '7500mb'
         }
         aucell {
            output = 'aucell_output.loom'
            rank_threshold = '5000'
            auc_threshold = '0.05'
            nes_threshold = '3.0'
            min_genes_regulon = 5
            min_regulon_gene_occurrence = 5
            pmem = '7500mb'
         }
         aggregate_features {
            use_chunking = true
            output_format = 'csv'
            compression = 'gzip'
         }
         numRuns = 100
      }
      scope {
         genome = ""
         tree {
            level_1 = "Human_Melanoma"
            level_2 = "10x"
            level_3 = ""
         }
      }
   }
}

process {
   executor = 'local'
}

singularity {
   enabled = true
   autoMounts = true
   runOptions = '-B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ -B /staging/leuven/res_00001/databases/:/staging/leuven/res_00001/databases/ -B /ddn1/vol1/staging/leuven/res_00001/databases/:/ddn1/vol1/staging/leuven/res_00001/databases/'
}
```

Run command,
```
nextflow -C nextflow.config run \
   -r scenic_multi_runs_fixes \
   aertslab/SingleCellTxBenchmark \
      -entry scenic
```

## MM057 SOX w/o TL

Nextflow config file:
```
manifest {
   name = 'aertslab/SingleCellTxBenchmark'
   description = 'A repository of pipelines for single-cell data in Nextflow DSL2'
   homePage = 'https://github.com/aertslab/SingleCellTxBenchmark'
   version = '0.3.1'
   mainScript = 'main.nf'
   defaultBranch = 'master'
   nextflowVersion = '!19.09.0-edge'
}

params {
   global {
      project_name = 'Human Melanoma Phenotype Switching'
      outdir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/9.MM057_SOX_wo_TL/pySCENIC_100x/out'
      tenx_folder = ''
      tracedir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/9.MM057_SOX_wo_TL/pySCENIC_100x/out/pipeline_reports'
      qsubaccount = 'lp_symbiosys'
   }
   sc {
      file_converter {
         iff = '10x_mtx'
         off = 'h5ad'
         useFilteredMatrix = true
      }
      file_annotator {
         iff = '10x_mtx'
         off = 'h5ad'
         type = 'sample'
         metaDataFilePath = ''
      }
      file_concatenator {
         join = 'outer'
         iff = '10x_mtx'
         off = 'h5ad'
      }
      star_concatenator {
         stranded = 'no'
         off = 'tsv'
      }
      scenic {
         container = 'docker://dweemx/sctx-pyscenic:0.9.19-4'
         scenicoutdir = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/9.MM057_SOX_wo_TL/pySCENIC_100x/out'
         filteredloom = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/aertslab/MEL/analysis_40k/20.analysis_with_sub_matrices/9.MM057_SOX_wo_TL/pySCENIC_100x/data/MM057_SOX_wo_TL_logCPM.loom'
         scenicOutputLoom = 'SCENIC_output.loom'
         scenicScopeOutputLoom = 'SCENIC_SCope_output.loom'
         numWorkers = '16'
         maxForks = 20
         mode = 'dask_multiprocessing'
         client_or_address = ''
         cell_id_attribute = 'CellID'
         gene_attribute = 'Gene'
         grn {
            TFs = '/ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/documents/resources/tf/allTFs_hg38.txt'
            seed = ''
            pmem = '3750mb'
         }
         cistarget {
            adj = 'adj.tsv'
            mtfDB = '/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-*.feather'
            mtfANN = '/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
            trkDB = '/staging/leuven/res_00001/databases/cistarget/databases/homo_sapiens/hg19/refseq_r45/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg19-*.feather'
            trkANN = '/ddn1/vol1/staging/leuven/stg_00002/icistarget-data/annotations/homo_sapiens/hg19/track_annotations/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg19.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv'
            type = ''
            output = 'reg.csv'
            rank_threshold = '5000'
            auc_threshold = '0.05'
            nes_threshold = '3.0'
            min_orthologous_identity = '0.0'
            max_similarity_fdr = '0.001'
            annotations_fname = ''
            thresholds = '0.75,0.90'
            top_n_targets = '50'
            top_n_regulators = '5,10,50'
            min_genes = '20'
            pmem = '7500mb'
         }
         aucell {
            output = 'aucell_output.loom'
            rank_threshold = '5000'
            auc_threshold = '0.05'
            nes_threshold = '3.0'
            pmem = '7500mb'
            min_genes_regulon = 5
            min_regulon_gene_occurrence = 5
         }
         aggregate_features {
            use_chunking = true
            output_format = 'csv'
            compression = 'gzip'
         }
         numRuns = 100
      }
      scope {
         genome = ''
         tree {
            level_1 = "Human_Melanoma"
            level_2 = "10x"
            level_3 = ""
         }
      }
   }
}

process {
   executor = 'local'
}

singularity {
   enabled = true
   autoMounts = true
   runOptions = '-B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ -B /staging/leuven/res_00001/:/staging/leuven/res_00001/ -B /ddn1/vol1/staging/leuven/res_00001/:/ddn1/vol1/staging/leuven/res_00001/'
}
```

Run command,
```
nextflow -C nextflow.config run \
   -r scenic_multi_runs_fixes \
   aertslab/SingleCellTxBenchmark
```
