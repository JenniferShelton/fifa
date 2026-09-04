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


workflow ExtractasticWkf {
    input {
        Bram bram
        String sampleId
        String projectId
        IndexedVcf vcf
        IndexedReference referenceFa
        # resources
        Int maxSplits = 40
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }

    call fifaTasks.MobsterFit {
        input:
            sampleId = sampleId,
            vcf = vcf.vcf,
            qos = qos,
            partition = partition,
            cpuPlatform = cpuPlatform
    }

    # gather split VCF filenames to avoid glob that can fail on prem
    Int vcfDiskSize = ceil(size(vcf.vcf, "MB")) 
    if ( vcfDiskSize > 10) {
        Int highMaxSplits = 100
    }
    Int count = select_first([highMaxSplits, maxSplits])
    scatter (i in range(count + 1)) {
        Int num = i + 1
        String suffixes = "${num}"
    }
    String prefix = "~{sampleId}.extractastic."
    String additionalSuffix = ".vcf"
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String splitVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
    }
    call fifaTasks.SplitVcf {
        input:
            vcf = vcf.vcf,
            prefix = prefix,
            diskSize = (ceil(size(vcf.vcf, "GB")) * 3) + 10,
            maxSplits = count,
            splitVcfPaths = splitVcfPaths
    }
    Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)

    Int diskSize = ceil(size(bram.bram, "GB")) * 3

    scatter (splitVcf in splitVcfs) {
        call fifaTasks.CompressIndexVcf {
            input:
                vcf = splitVcf
        }
        call fifaTasks.ExtractionWithMobsterFit {
            input:
                bram = bram,
                sampleId = sampleId,
                projectId = projectId,
                vcf = CompressIndexVcf.vcfCompressedIndexed,
                referenceFa = referenceFa,
                mobsterFitRds = MobsterFit.mobsterFitRds,
                diskSize = diskSize,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
    }

    call fifaTasks.ConcateTables {
        input:
            tables = ExtractionWithMobsterFit.extractedFeatures,
            outputTablePath = "~{sampleId}_extracted_features.csv"
    }

    output {
        File extractedFeatures = ConcateTables.outputTable
        File mobsterFitRds = MobsterFit.mobsterFitRds
    }
}
