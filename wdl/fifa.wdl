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
        Int threads = 2
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
        referenceFa="~{referenceFa.fasta}"
        # Extract features for the EBM model. This script also extracts the reference sequence context for each variant, which is used in the mutational signature analysis.
        python3 fifa \
            extract \
            -s ${sampleId} \
            -c ${projectId} \
            -v ${vcff} \
            -b ${bam} \
            -r ${referenceFa} \
            -o . \
            -n 2
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
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:363b97fcd9896b4c89ce1fe05df8f7a60ed60dbd61588e905333c432d4555ea2"
        runtime_minutes: "6000"
        partition: "cpu"
        qos: qos
    }
}    
