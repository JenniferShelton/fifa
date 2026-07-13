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
import "wdl/prediction_wkf.wdl" as predictionWkf
import "wdl/fifa.wdl" as fifaTasks

workflow PredictionBucketedWkf {
    input {
        Boolean local = true
        String bamFilenamePath
        Int bamFileSize
        String bamDownloadUri
        String baiFilenamePath
        String baiDownloadUri
        String vcfFilenamePath
        String vcfDownloadUri
        String vcfIndexFilenamePath
        String vcfIndexDownloadUri
        File awsConfig
        File awsCredentials
        String endpointUrl
        String sampleId
        String tumorId
        String normalId
        String projectId
        IndexedReference referenceFa
        File? optionalRnaFile
        Array[File] models
        # resources
        String? gcpProject 
        File? serviceAccountKey
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    if (!local) {
        Int cloudMaxSplits = 5
        Int cloudMinSplits = 2
    }
    if (local) {
        Int hpcMaxSplits = 4
        Int hpcMinSplits = 2
    }

    call fifaTasks.Download as bamDownload {
            input:
                filenamePath = bamFilenamePath,
                downloadUri = bamDownloadUri,
                awsConfig = awsConfig,
                awsCredentials = awsCredentials,
                endpointUrl = endpointUrl,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform,
                diskSize = bamFileSize + 20
        }

    call fifaTasks.Download as baiDownload {
            input:
                filenamePath = baiFilenamePath,
                downloadUri = baiDownloadUri,
                awsConfig = awsConfig,
                awsCredentials = awsCredentials,
                endpointUrl = endpointUrl,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform,
                diskSize = 4
        }
    call fifaTasks.Download as vcfDownload {
            input:
                filenamePath = vcfFilenamePath,
                downloadUri = vcfDownloadUri,
                awsConfig = awsConfig,
                awsCredentials = awsCredentials,
                endpointUrl = endpointUrl,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
    
    call fifaTasks.Download as vcfIndexDownload {
            input:
                filenamePath = vcfIndexFilenamePath,
                downloadUri = vcfIndexDownloadUri,
                awsConfig = awsConfig,
                awsCredentials = awsCredentials,
                endpointUrl = endpointUrl,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }

    IndexedVcf vcf = object {
                vcf: vcfDownload.download,
                index: vcfIndexDownload.download
            }
    Bam bam = object {
                bam : bamDownload.download,
                bamIndex : baiDownload.download
            }
    call fifaTasks.MakeVariantBed {
        input:
            vcf = vcf.vcf,
            sampleId = sampleId,
            referenceFa = referenceFa,
            diskSize = 10
    }
    call fifaTasks.MakeVariantCram {
        input:
            finalBam = bam,
            gcpProject = gcpProject,
            serviceAccountKey = serviceAccountKey,
            features1000Bed = MakeVariantBed.features1000Bed,
            sampleId = sampleId,
            referenceFa = referenceFa,
            diskSize = 30
    }
    Bram bram = object {
        bram : MakeVariantCram.variantCram.cram,
        bramIndex : MakeVariantCram.variantCram.cramIndex
    }
    if (defined(optionalRnaFile)) {
        File rnaFile = select_first([optionalRnaFile])
        call predictionWkf.PredictionWkf as rnaPredictionWkf {
            input:
                bram = bram,
                sampleId = sampleId,
                projectId = projectId,
                vcf = vcf,
                optionalRnaFile = rnaFile,
                models = models,
                referenceFa = referenceFa,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
    }
    if (!defined(optionalRnaFile)) {
        call predictionWkf.PredictionWkf {
            input:
                bram = bram,
                sampleId = sampleId,
                projectId = projectId,
                vcf = vcf,
                models = models,
                referenceFa = referenceFa,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform
        }
    }
    File extractedFeaturesRun = select_first([rnaPredictionWkf.extractedFeatures, PredictionWkf.extractedFeatures])
    File fifaVcfRun = select_first([rnaPredictionWkf.fifaVcf, PredictionWkf.fifaVcf])
    # }
    # call fifaTasks.ConcateTables {
    #     input:
    #         tables = extractedFeaturesRun,
    #         outputTablePath =  "~{sampleId}.extracted_features.csv"
    # }
    # call fifaTasks.Gatk4MergeSortVcf {
    #     input:
    #         tempVcfs = fifaVcfRun,
    #         referenceFa = referenceFa,
    #         sortedVcfPath = "~{sampleId}.fifa.vcf"
    # }
    output {
        File extractedFeatures = extractedFeaturesRun
        File fifaVcf = fifaVcfRun
    }

}