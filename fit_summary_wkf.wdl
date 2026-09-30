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


workflow FitPrintWkf {
    input {
        String projectId
        Array[String] sampleIds
        Array[IndexedVcf] vcfs
        Array[File] extractedFeaturesFifas
        # resources
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }

    scatter (i in range(length(sampleIds))){
        call fifaTasks.MobsterFitCommpressed {
            input:
                sampleId = sampleIds[i],
                vcf = vcfs[i],
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
        if (length(extractedFeaturesFifas) > 0) {
            call fifaTasks.MergeMobsterFit {
                input:
                    sampleId = sampleIds[i],
                    extractedFeaturesFifa = extractedFeaturesFifas[i],
                    fitProbCsv = MobsterFitCommpressed.fitProbCsv,
                    fitProbCsvK2 = MobsterFitCommpressed.fitProbCsvK2
        } 
    }

    }

    call fifaTasks.ConcateTables {
        input:
            tables = MobsterFitCommpressed.fitProbCsv,
            outputTablePath = "~{projectId}.mobster_fit.prob.csv"
    }

    call fifaTasks.ConcateTables as concateTablesK2 {
        input:
            tables = MobsterFitCommpressed.fitProbCsvK2,
            outputTablePath = "~{projectId}.mobster_fit.k2.prob.csv"
    }
    if (length(extractedFeaturesFifas) > 0) {
        Array[File] extractedFeaturesTables = select_all(MergeMobsterFit.extractedFeatures)
        call fifaTasks.ConcateTables as concateTablesFeatures {
            input:
                tables = extractedFeaturesTables,
                outputTablePath = "~{projectId}_extracted_features.csv"
        }
    }
        

    output {
        File fitProbCsv = ConcateTables.outputTable
        File fitProbCsvK2 = concateTablesK2.outputTable
        File? extractedFeatures = concateTablesFeatures.outputTable
        Array[File] fitPng = MobsterFitCommpressed.fitPng
        Array[File] fitPngK2 = MobsterFitCommpressed.fitPngK2
        Array[File] mobsterFitRds = MobsterFitCommpressed.mobsterFitRds
    }
}
