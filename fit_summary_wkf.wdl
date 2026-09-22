version 1.0

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2026) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Nico Robine (nrobine@nygenome.org)
#    Will Liao (wliao@nygenome.org)
#    Valentina Grether
#    Zoe R. Goldstein (zgoldstein@nygenome.org)
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Timothy R. Chu (tchu@nygenome.org)
#    William F. Hooper (whooper@nygenome.org)
#    Heather Geiger (hgeigher @nygenome.org)
#    André Corvelo (acorvelo@nygenome.org)
#    Rachel Martini 
#    Melissa B. Davis
# 
#
# ================== /COPYRIGHT ===============================================

# Workflow from https://www.biorxiv.org/content/10.64898/2026.03.10.710815v1
# An explainable boosting machine model for identifying artifacts caused by formalin-fixed paraffin embedding
import "wdl/wdl_structs.wdl"
import "wdl/fifa.wdl" as fifaTasks



task DescribeMobsterFit {
    input {
        String sampleId
        File mobsterFitRds
        String fitPngPath = "~{sampleId}.mobster_fit.png"
        String fitCsvPath = "~{sampleId}.mobster_fit.csv"
        File fitPrintRscript = "/gpfs/commons/groups/compbio/projects/FFPE_filtering/repos/fifa/src/run_fit_print.R"
        # resources
        Int threads = 1
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 24
        Int diskSize = 20
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"

    }
    command <<<
        set -e -o pipefail
        sampleId="~{sampleId}"
        mobsterFitRds="~{mobsterFitRds}"
        fitPngPath="~{fitPngPath}"
        fitCsvPath="~{fitCsvPath}"
        

        # Extract features for the EBM model, reusing a precomputed MOBSTER fit instead of refitting per shard.
        Rscript ~{fitPrintRscript} \
            ~{sampleId} \
            ~{mobsterFitRds} \
            ~{fitCsvPath} \
            ~{fitPngPath}
    >>>

    output {
        File fitCsv = fitCsvPath
        File fitPng = fitPngPath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:98c2f925924537525d1f08224ca6a0cecb253e1defc6f769743310ac7da86bda"
        runtime_minutes: "6000"
        cpuPlatform : cpuPlatform
        partition: "cpu"
        qos: qos
    }
}


workflow FitPrintWkf {
    input {
        String projectId
        Array[String] sampleIds
        Array[IndexedVcf] vcfs
        # resources
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }

    scatter (i in range(length(sampleIds))){
        call fifaTasks.MobsterFit {
            input:
                sampleId = sampleIds[i],
                vcf = vcfs[i].vcf,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
        call DescribeMobsterFit {
            input:
                sampleId = sampleIds[i],
                mobsterFitRds =  MobsterFit.mobsterFitRds,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
    }

    call fifaTasks.ConcateTables {
        input:
            tables = DescribeMobsterFit.fitCsv,
            outputTablePath = "~{projectId}.mobster_fit.csv"
    }

    output {
        File fitCsv = ConcateTables.outputTable
        Array[File] fitPng = DescribeMobsterFit.fitPng
        Array[File] mobsterFitRds = MobsterFit.mobsterFitRds
    }
}
