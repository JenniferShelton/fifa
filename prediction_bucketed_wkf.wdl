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
import "wdl/tasks.wdl" as tasks

workflow PredictionBucketedWkf {
    input {
        File countHetsScript
        File hetSnpsScript
        Array[File] gnomADMaf
        Boolean local = false
        String bamFilenamePath
        Int bamFileSize
        Int normalBamFileSize
        String bamDownloadUri
        String baiFilenamePath
        String baiDownloadUri
        String normalBamFilenamePath
        String normalBamDownloadUri
        String normalBaiFilenamePath
        String normalBaiDownloadUri
        String vcfFilenamePath
        String vcfDownloadUri
        String vcfIndexFilenamePath
        String vcfIndexDownloadUri
        File awsConfig
        File awsCredentials
        String endpointUrl
        String pairId
        String tumorId
        String normalId
        String projectId
        IndexedReference referenceFa
        File? optionalRnaFile
        Array[File] models
        Boolean mobsterFree = true
        # resources
        String? gcpProject 
        File? serviceAccountKey
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    if (!local) {
        Int cloudMaxSplits = 30
        Int cloudMinSplits = 28
    }
    if (local) {
        Int hpcMaxSplits = 30
        Int hpcMinSplits = 28
    }
    # tumor
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
    # 
    call fifaTasks.Download as normalBamDownload {
            input:
                filenamePath = normalBamFilenamePath,
                downloadUri = normalBamDownloadUri,
                awsConfig = awsConfig,
                awsCredentials = awsCredentials,
                endpointUrl = endpointUrl,
                qos = qos,
                partition = partition,
                cpuPlatform = cpuPlatform,
                diskSize = normalBamFileSize + 20
        }

    call fifaTasks.Download as normalBaiDownload {
            input:
                filenamePath = normalBaiFilenamePath,
                downloadUri = normalBaiDownloadUri,
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
    Bram bram = object {
                bram : bamDownload.download,
                bramIndex : baiDownload.download
            }
    Bram normalBram = object {
                bram : normalBamDownload.download,
                bramIndex : normalBaiDownload.download
            }
    # gather split VCF filenames to avoid glob that can fail on prem
    Int count = 100
    scatter (i in range(count)) {
        Int num = i + 1
        String suffixes = "${num}"
    }
    String prefix = "~{pairId}.fifa."
    String additionalSuffix = ".vcf"
    scatter (index in range(length(suffixes))) {
        String suffix = suffixes[index]
        String splitVcfPaths = "~{prefix}.~{suffix}~{additionalSuffix}"
    }
    Int passDiskSize = (ceil( size(vcf.vcf, "GB") )  * 2 ) + 10
    call tasks.FilterForPassSnps {
        input:
            passVcfPath = "~{pairId}.pass.snps.vcf.gz",
            unFilteredVcf = vcf,
            diskSize = passDiskSize
    }
    call fifaTasks.SplitVcf {
        input:
            vcf = FilterForPassSnps.passVcf.vcf,
            prefix = prefix,
            diskSize = (ceil(size(FilterForPassSnps.passVcf.vcf, "GB")) * 3) + 10,
            maxRows = 1000,
            minSplits = select_first([cloudMinSplits, hpcMinSplits]),
            maxSplits = select_first([cloudMaxSplits, hpcMaxSplits]),
            splitVcfPaths = splitVcfPaths
    }
    Array[File] splitVcfs = select_all(SplitVcf.splitVcfs)
    scatter(splitVcf in splitVcfs) {
        call fifaTasks.ReorderVcfColumns {
            input:
                tumor = tumorId,
                normal = normalId,
                rawVcf = splitVcf,
                orderedVcfPath = "~{pairId}.renamedColumns.vcf",
                memoryGb = 4,
                diskSize = 10
        }

        call fifaTasks.CompressIndexVcf {
            input:
                vcf = ReorderVcfColumns.orderedVcf,
                memoryGb = 1
        }
        call fifaTasks.MakeVariantBed {
            input:
                vcf = splitVcf,
                sampleId = pairId,
                referenceFa = referenceFa,
                diskSize = 10
        }
        call fifaTasks.MakeVariantCram {
            input:
                finalBram = bram,
                gcpProject = gcpProject,
                serviceAccountKey = serviceAccountKey,
                features1000Bed = MakeVariantBed.features1000Bed,
                sampleId = pairId,
                referenceFa = referenceFa,
                diskSize = 30
        }
        Bram tumorVariantBram = object {
            bram : MakeVariantCram.variantCram.cram,
            bramIndex : MakeVariantCram.variantCram.cramIndex
        }
        File tumorVariantBrams = MakeVariantCram.variantCram.cram
        if (defined(optionalRnaFile)) {
            File rnaFile = select_first([optionalRnaFile])
            call predictionWkf.PredictionWkf as rnaPredictionWkf {
                input:
                    bram = tumorVariantBram,
                    sampleId = tumorId,
                    projectId = projectId,
                    vcf = CompressIndexVcf.vcfCompressedIndexed,
                    optionalRnaFile = rnaFile,
                    models = models,
                    mobsterFree = mobsterFree,
                    referenceFa = referenceFa,
                    qos = qos,
                    partition = partition,
                    cpuPlatform = cpuPlatform
            }
        }
        if (!defined(optionalRnaFile)) {
            call predictionWkf.PredictionWkf {
                input:
                    bram = tumorVariantBram,
                    sampleId = tumorId,
                    projectId = projectId,
                    mobsterFree = mobsterFree,
                    vcf = CompressIndexVcf.vcfCompressedIndexed,
                    models = models,
                    referenceFa = referenceFa,
                    qos = qos,
                    partition = partition,
                    cpuPlatform = cpuPlatform
            }
        }
        File extractedFeaturesRun = select_first([rnaPredictionWkf.extractedFeatures, PredictionWkf.extractedFeatures])
        File fifaVcfRun = select_first([rnaPredictionWkf.fifaVcf, PredictionWkf.fifaVcf])
    }
    call fifaTasks.ConcateTables {
        input:
            tables = extractedFeaturesRun,
            outputTablePath =  "~{pairId}.extracted_features.csv"
    }
    call fifaTasks.Gatk4MergeSortVcf {
        input:
            tempVcfs = fifaVcfRun,
            referenceFa = referenceFa,
            sortedVcfPath = "~{pairId}.fifa.vcf"
    }
    # merge tumor variant BAMs
    Int tumorVariantBramsSize = ceil(size(tumorVariantBrams, "GB"))
    call fifaTasks.MergeSortAlignments {
        input:
            inputs = tumorVariantBrams,
            outputPath = "~{pairId}.tumor.variant.bam",
            outputIndexPath = "~{pairId}.tumor.variant.bai",
            threads = 4,
            diskSize = (3 * tumorVariantBramsSize),
            referenceFa = referenceFa
    }

    call fifaTasks.MakeVariantBed as normalMakeVariantBed {
            input:
                vcf = vcf.vcf,
                sampleId = normalId,
                referenceFa = referenceFa,
                diskSize = 10
        }
    call fifaTasks.MakeVariantCram as normalMakeVariantCram {
        input:
            finalBram = normalBram,
            gcpProject = gcpProject,
            serviceAccountKey = serviceAccountKey,
            features1000Bed = normalMakeVariantBed.features1000Bed,
            sampleId = normalId,
            referenceFa = referenceFa,
            diskSize = 30
    }
    call tasks.FragCounter as normalFragCounter {
        input:
            bram = normalBram,
            referenceFa = referenceFa,
            sampleId = normalId,
            diskSize = normalBamFileSize + 10
    }
    call tasks.FragCounter as tumorFragCounter {
        input:
            bram = bram,
            referenceFa = referenceFa,
            sampleId = tumorId,
            diskSize = bamFileSize + 10
    }

    scatter (gnomADMaf in gnomADMafs) {
        call fifaTasks.MakeVariantCram as tumor2MakeVariantCram {
            input:
                finalBram = bram,
                gcpProject = gcpProject,
                serviceAccountKey = serviceAccountKey,
                features1000Bed = gnomADMaf,
                sampleId = tumorId,
                referenceFa = referenceFa,
                diskSize = 30
        }
        call fifaTasks.MakeVariantCram as normal2MakeVariantCram {
            input:
                finalBram = normalBram,
                gcpProject = gcpProject,
                serviceAccountKey = serviceAccountKey,
                features1000Bed = gnomADMaf,
                sampleId = normalId,
                referenceFa = referenceFa,
                diskSize = 30
        }
        Int tumor2VariantBramsSize = ceil(size(tumor2MakeVariantCram.variantCram.cram, "GB"))
        Int normal2VariantBramsSize = ceil(size(normal2MakeVariantCram.variantCram.cram, "GB"))
        
        Bram tumor2VariantBram = object {
            bram : tumor2MakeVariantCram.variantCram.cram,
            bramIndex : tumor2MakeVariantCram.variantCram.cramIndex
        }
        Bram normal2VariantBram = object {
            bram : normal2MakeVariantCram.variantCram.cram,
            bramIndex : normal2MakeVariantCram.variantCram.cramIndex
        }
        call tasks.CountHetsSnps as countHetsSnps {
            input:
                countHetsScript = countHetsScript,
                hetSnpsScript = hetSnpsScript,
                gnomADMaf = gnomADMaf,
                tumorBram = tumor2VariantBram,
                normalBram = normal2VariantBram,
                pairId = pairId,
                diskSize = tumor2VariantBramsSize + normal2VariantBramsSize + 20,
                memoryGb = 16
        }
    }

    call tasks.ConcateTables as ConcateTablesHs {
        input :
            outputTablePath = "hetSnpsCounts.tsv",
            tables = countHetsSnps.hetSnpsCountsTsv,
            diskSize = 20
    }
    # call tasks.CountHetsSnps as countHetsSnps {
    #     input:
    #         countHetsScript = countHetsScript,
    #         hetSnpsScript = hetSnpsScript,
    #         gnomADMaf = gnomADMaf,
    #         tumorBram = bram,
    #         normalBram = normalBram,
    #         pairId = pairId,
    #         diskSize = normalBamFileSize + bamFileSize + 10,
    #         memoryGb = 16
    # }

    call tasks.HapaSegLocal {
        input:
            tumorBram = bram,
            normalBram = normalBram,
            pairId = pairId,
            hetSnpsCountsTsv = ConcateTablesHs.outputTable,
            diskSize = normalBamFileSize + bamFileSize + 30,
            memoryGb = 24,
            threads = 8
    }
    output {
        Array[File] hapaSegLocalFiles = HapaSegLocal.hapaSegLocalFiles
        File extractedFeatures = ConcateTables.outputTable
        File fifaVcf = Gatk4MergeSortVcf.sortedVcf.vcf
        Bram tumorVariantCram = MergeSortAlignments.mergedBram
        Cram normalVariantCram = normalMakeVariantCram.variantCram
        Array[File] normalRds = normalFragCounter.rds 
        Array[File] tumorRds = tumorFragCounter.rds
        # File normalRds = normalFragCounter.rds
        # File normalRawRds = normalFragCounter.rawRds
        File normalCorrectedBw = normalFragCounter.correctedBw
        # File normalOneKbRds = normalFragCounter.oneKbRds
        # File tumorRds = tumorFragCounter.rds
        # File tumorRawRds = tumorFragCounter.rawRds
        File tumorCorrectedBw = tumorFragCounter.correctedBw
        # File tumorOneKbRds = tumorFragCounter.oneKbRds
    }
}