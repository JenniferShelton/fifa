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
# additional tasks for the ALCHEMIST pipeline

import "wdl_structs.wdl"



task FragCounter {
    input {
        Bram bram
        String bramPath = basename(bram.bram)
        String bramIndexPath = basename(bram.bramIndex)
        IndexedReference referenceFa
        String sampleId
        String rdsPath="~{sampleId}.cov.rds"
        String rawRdsPath="~{sampleId}.cov.raw.rds"
        String correctedBwPath="~{sampleId}.cov.corrected.bw"
        String oneKbRdsPath="~{sampleId}.1kb_cov.rds"
        Int diskSize = 20
        Int memoryGb = 16
    }

    command <<<
        set -e -o pipefail

        ln -s \
        "~{bram.bram}" \
        "~{bramPath}"
        ln -s \
        "~{bram.bramIndex}" \
        "~{bramIndexPath}"

        bram="~{bram.bram}"
        referenceFa="~{referenceFa.fasta}"
        sampleId="~{sampleId}"
        rdsPath="~{sampleId}.cov.rds"
        rawRdsPath="~{sampleId}.cov.raw.rds"
        correctedBwPath="~{sampleId}.cov.corrected.bw"
        oneKbRdsPath="~{sampleId}.1kb_cov.rds"

        frag \
        --bam "~{bramPath}" \
        --window 1000 \
        --gcmapdir /usr/local/lib/R/site-library/fragCounter/extdata/gcmap.38 \
        --paired TRUE \
        --reference ${referenceFa} \
        --outdir . \
        --libdir /usr/local/lib/R/site-library/fragCounter/extdata/

        find . -name "*.rds"

        mv \
        cov.rds ${rdsPath}
        mv \
        cov.raw.rds ${rawRdsPath}
        mv \
        cov.corrected.bw ${correctedBwPath} 
        # mv \
        # 1kb_cov.rds ${oneKbRdsPath}
    >>>

    output {
        Array[File] rds = glob("*.rds")
        # File rawRds = "~{rawRdsPath}"
        File correctedBw = "~{correctedBwPath}"
        # File oneKbRds = "~{oneKbRdsPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fragcounter@sha256:1f04990656cb7a6aecb763c0684caf5d885f66c7fc5f720014f3e68baf60806b"
        runtime_minutes: "90"
    }
}


task CountHetsSnps {
    input {
        File countHetsScript
        File hetSnpsScript
        File gnomADMaf
        Bram tumorBram
        Bram normalBram
        String tumorBramPath = basename(tumorBram.bram)
        String tumorBramIndexPath = basename(tumorBram.bramIndex)
        String normalBramPath = basename(normalBram.bram)
        String normalBramIndexPath = basename(normalBram.bramIndex)
        Int diskSize
        Int memoryGb
        String pairId
        String hetSnpsCountsTsvPath = "~{pairId}.het_snps_counts.tsv"
    }

    command <<<
        set -e -o pipefail

        hetSnpsScript=~{hetSnpsScript}
        gnomADMaf=~{gnomADMaf}
        tumorBram=~{tumorBram.bram}
        normalBram=~{normalBram.bram}
        hetSnpsCountsTsvPath=~{hetSnpsCountsTsvPath}

        ln -s \
        "~{tumorBram.bram}" \
        "~{tumorBramPath}"
        ln -s \
        "~{tumorBram.bramIndex}" \
        "~{tumorBramIndexPath}"
        ln -s \
        "~{normalBram.bram}" \
        "~{normalBramPath}"
        ln -s \
        "~{normalBram.bramIndex}" \
        "~{normalBramIndexPath}"

        python \
        ${hetSnpsScript} \
        --het-sites ${gnomADMaf} \
        --tumor-bam ~{tumorBramPath} \
        --normal-bam ~{normalBramPath} \
        --out ${hetSnpsCountsTsvPath}
    >>>

    output {
        File hetSnpsCountsTsv = "~{hetSnpsCountsTsvPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}



task HapaSegLocal {
    input {
        Bram tumorBram
        Bram normalBram
        String pairId
        File hetSnpsCountsTsv
        Int diskSize
        Int memoryGb = 24
        Int threads = 8
    }

    command <<<
        set -e -o pipefail

        hetSnpsCountsTsv=~{hetSnpsCountsTsv}
        tumorBram=~{tumorBram.bram}
        normalBram=~{normalBram.bram}

        hapaseg_local \
        --is-ffpe  \
        --ref-root-path $(pwd) \
        --ref-genome-build hg38 \
        --tumor-bam ${tumorBram} \
        --tumor-bai ~{tumorBram.bramIndex} \
        --normal-bam ${normalBram} \
        --normal-bai ~{normalBram.bramIndex} \
        --max-cpus ~{threads} \
        --max-memory ~{memoryGb} \
        --het-site-calls ${hetSnpsCountsTsv} \
        $(pwd) \
        ~{pairId}

        ## Clear out reference data
        rm -r hg38/ hg19/ || true

    >>>

    output {
        Array[File] hapaSegLocalFiles = glob("~{pairId}*")
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/hapaseg_local@sha256:29345c386af52082caf0621be73fbbac4b70a88435b8c64f8262b84d7716cd5f"
        runtime_minutes: "90"
    }
}


task FilterForPassSnps {
    input {
        String passVcfPath
        IndexedVcf unFilteredVcf

        Int memoryGb = 24
        Int diskSize = (ceil( size(unFilteredVcf.vcf, "GB") )  * 2 ) + 10
    }

    command {
        set -e -o pipefail

        bcftools view \
        -f PASS \
        ~{unFilteredVcf.vcf} \
        --types snps \
        | bgzip -c \
        > ~{passVcfPath}

        tabix \
        -p vcf \
        ~{passVcfPath}
    }

    output {
        IndexedVcf passVcf = object {
            vcf : "~{passVcfPath}",
            index : "~{passVcfPath}.tbi"
        }
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task ConcateTables {
    input {
        String outputTablePath
        Array[File] tables
        Int diskSize = 1
        Int memoryGb = 4
    }
    
    command {
       
       
       python3 \
        /concate_tables.py \
        --tables ~{sep=" " tables} \
        --output ~{outputTablePath}
    }
    
    output {
        File outputTable = "~{outputTablePath}"
    }
    
    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
    
    meta {
        internalOnly : "False, produces output for external or internal CSV files"
    }
}