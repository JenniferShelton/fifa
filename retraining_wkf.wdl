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
import "wdl/wdl_structs.wdl"
import "wdl/retraining_wkf.wdl" as retrainingExtract
import "wdl/fifa.wdl" as fifa


workflow RetrainingWkfs {
    input {
        Array[Bam] bams
        Array[String] sampleIds
        String projectId
        Array[IndexedVcf] vcfs
        IndexedReference referenceFa
        File labels
        String modelPath = "fifa_model.pkl"
        Boolean hyperParameter = false
        # resources
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    scatter (i in range(length(bams))) {
        call retrainingExtract.RetrainingExtractWkf {
            input:
                bam = bams[i],
                sampleId = sampleIds[i],
                projectId = projectId,
                vcf = vcfs[i],
                referenceFa = referenceFa,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
    }
    if (hyperParameter) {
        String hyperParameterFlagTrue = " --hyperparameter"
    }
    String hyperParameterFlag = select_first([hyperParameterFlagTrue, ""])
    call fifa.ReTraining {
        input:
            extractedFeatures = RetrainingExtractWkf.extractedFeatures,
            labels = labels,
            modelPath = modelPath,
            hyperParameterFlag = hyperParameterFlag
    }
    output {
        Array[File] extractedFeatures = RetrainingExtractWkf.extractedFeatures
        File model = ReTraining.model
    }
}
