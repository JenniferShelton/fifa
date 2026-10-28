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
import "wdl_structs.wdl"
import "fifa.wdl" as fifa


workflow PredictionWkf {
    input {
        File extractedFeatures
        String sampleId
        String projectId
        IndexedVcf vcf
        IndexedReference referenceFa
        File? optionalRnaFile
        Array[File] models
        Boolean mobsterFree = false
        # resources
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
        Int maxSplits = 10
        Int minSplits = 2
    }
    # gather split VCF filenames to avoid glob that can fail on prem
    Int count = 100
    scatter (i in range(count + 1)) {
        Int num = i + 1
        String suffixes = "${num}"
    }
    String pairId = sub(basename(vcf.vcf, ".gz"), "\\..*$", "")
    String tumorId = sub(pairId, "^.*--", "")
    String normalId = sub(pairId, "--.*$", "")
    String prefix = "~{pairId}.fifa."
    String additionalSuffix = ".vcf"
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String splitVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
    }
    call fifa.SplitVcf {
        input:
            vcf = vcf.vcf,
            prefix = prefix,
            diskSize = (ceil(size(vcf.vcf, "GB")) * 3) + 10,
            maxRows = 1000,
            minSplits = minSplits,
            maxSplits = maxSplits,
            splitVcfPaths = splitVcfPaths
    }
    Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)
    scatter(splitVcf in splitVcfs) {
        call fifa.ReorderVcfColumns {
            input:
                tumor = tumorId,
                normal = normalId,
                rawVcf = splitVcf,
                orderedVcfPath = "~{pairId}.renamedColumns.vcf",
                memoryGb = 4,
                diskSize = 10
        }

        call fifa.CompressIndexVcf {
            input:
                vcf = ReorderVcfColumns.orderedVcf,
                memoryGb = 1
        }
        if (defined(optionalRnaFile)) {
            call fifa.PredictionWithRna {
                input:
                    sampleId = sampleId,
                    vcf = CompressIndexVcf.vcfCompressedIndexed,
                    extractedFeatures = extractedFeatures,
                    fifaCsvPath = "~{pairId}.fifa.csv",
                    models = models,
                    optionalRnaFile = optionalRnaFile
            }
        }

        if (!defined(optionalRnaFile)) {
            call fifa.Prediction {
                input:
                    sampleId = sampleId,
                    vcf = CompressIndexVcf.vcfCompressedIndexed,
                    extractedFeatures = extractedFeatures,
                    fifaCsvPath = "~{pairId}.fifa.csv",
                    models = models
            }
        }
        File fifaVcfRun = select_first([PredictionWithRna.fifaVcf, Prediction.fifaVcf])
        File fifaCsvRun = select_first([PredictionWithRna.fifaCsv, Prediction.fifaCsv])
    }
    call fifa.Gatk4MergeSortVcf {
        input:
            tempVcfs = fifaVcfRun,
            referenceFa = referenceFa,
            sortedVcfPath = "~{pairId}.fifa.vcf"
    }
    output {
        File fifaVcf = Gatk4MergeSortVcf.sortedVcf.vcf
        File fifaCsv = fifaCsvRun[0]
    }
}
