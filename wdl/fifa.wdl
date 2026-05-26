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
import "wdl_structs.wdl"

task Extraction {
    input {
        Bam bam
        String sampleId
        String projectId
        IndexedVcf vcf
        IndexedReference referenceFa
        String extractedFeaturesPath = "~{sampleId}_extracted_features.csv"
        # resources
        Int threads = 4
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize
        String qos = "compbio"
        String partition = "cpu"

    }
    command <<<
        set -e -o pipefail
        sampleId="~{sampleId}"
        projectId="~{projectId}"
        vcf="~{vcf.vcf}"
        threads="~{threads}"
        referenceFa="~{referenceFa.fasta}"
        # Extract features for the EBM model. This script also extracts the reference sequence context for each variant, which is used in the mutational signature analysis.
        fifa \
            extract \
            -s ${sampleId} \
            -c ${projectId} \
            -v ${vcf} \
            -b ${bam} \
            -r ${referenceFa} \
            -o . \
            -n ${threads}
    >>>

    output {
        File extractedFeatures = extractedFeaturesPath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "6000"
        partition: "cpu"
        qos: qos
    }
}    

task Merge {
    input {
        Array[File] pickles
        String mergedModelId
        String modelPath = "~{mergedModelId}_extracted_features.csv"
        # resources
        Int threads = 1
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize = 20
        String qos = "compbio"
        String partition = "cpu"

    }
    command <<<
        set -e -o pipefail
        mergedModelId="~{mergedModelId}"
        modelPath="~{modelPath}"
        all_pickles=~{sep=' ' pickles}
        # Extract features for the EBM model. This script also extracts the reference sequence context for each variant, which is used in the mutational signature analysis.
        fifa \
            merge \
            -m ${all_pickles} \
            -o ${modelPath}
    >>>

    output {
        File model = modelPath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " LOCAL"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "300"
        partition: "cpu"
        qos: qos
    }
}    