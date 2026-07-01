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


task ReorderVcfColumns {
    input {
        Int diskSize
        Int memoryGb
        String normal
        String tumor
        File rawVcf
        String orderedVcfPath
    }

    command {
        set -e -o pipefail

        if [[ ~{rawVcf} == *.vcf.gz ]]; then
            gunzip \
            -c \
            ~{rawVcf} \
            > ~{rawVcf}.vcf
            runVcf=~{rawVcf}.vcf
        elif [[ ~{rawVcf} == *.vcf.bgz ]]; then
            gunzip \
            -c \
            ~{rawVcf} \
            > ~{rawVcf}.vcf
            runVcf=~{rawVcf}.vcf
        else
            runVcf=~{rawVcf}
        fi

        python \
        /reorder_vcf.py \
        $runVcf \
        ~{orderedVcfPath} \
        ~{normal} ~{tumor}
    }

    output {
        File orderedVcf = "~{orderedVcfPath}"
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b"
        runtime_minutes: "90"
    }
}

task MakeVariantBed {
    input {
        File vcf
        String sampleId
        String featuresBedPath = "~{sampleId}.features.bed"
        IndexedReference referenceFa
        Int slopSizeBp = 1000
        String features1000BedPath = "~{sampleId}.features.1000.bed"
        # resources
        Int diskSize
        Int threads = 1
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
    }

    command <<<
        set -e -o pipefail

        cut -f1,2 \
        ~{referenceFa.index} \
        > sizes_genome.tsv

        python3 \
        /vcf_to_bed.py \
        ~{vcf} \
        | cut \
        -f 1,2,3 \
        | sort \
        -k 1,1 \
        -k2,2n \
        | bedtools \
        merge \
        > ~{featuresBedPath}

        bedtools \
        slop \
        -i ~{featuresBedPath} \
        -g sizes_genome.tsv \
        -b ~{slopSizeBp} \
        > ~{features1000BedPath}
    >>>
    output {
        File features1000Bed = "~{features1000BedPath}"
    }
    runtime {
        docker: "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        disks: "local-disk " + diskSize + " HDD"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu : threads
        cpus: runRequestThreads
        runtime_minutes: "600"
        maxRetries : 4
    }
}

task MakeVariantCram {
    input {
        Bam finalBam
        File features1000Bed
        String sampleId
        String variantCramPath = "~{sampleId}.variantRegions.cram"
        String variantCraiPath = "~{sampleId}.variantRegions.crai"
        IndexedReference referenceFa
        Int slopSizeBp = 1000
        # resources
        String? gcpProject 
        File? serviceAccountKey
        Int diskSize
        Int threads = 8
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 16
    }
    command <<<
        set -e -o pipefail
        serviceAccountKey=~{serviceAccountKey}
        if [ -f "$serviceAccountKey" ]; then
            export GCS_REQUESTER_PAYS_PROJECT=~{gcpProject}
            export GOOGLE_APPLICATION_CREDENTIALS=~{serviceAccountKey}
            export GOOGLE_CLOUD_PROJECT=~{gcpProject}
            # expires in 60 min
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        fi

        samtools view \
        --threads ~{threads} \
        -h \
        --cram \
        --reference ~{referenceFa.fasta} \
        -L ~{features1000Bed} \
        ~{finalBam.bam} \
        > ~{variantCramPath}

        samtools \
        index \
        -@ ~{threads} \
        ~{variantCramPath} \
        ~{variantCraiPath}
    >>>

    output {
        Cram variantCram = object {
                cram : "~{variantCramPath}",
                cramIndex : "~{variantCraiPath}"
             }
    }
    runtime {
        docker:  "us.gcr.io/nygc-comp-s-fd4e/samtoolsgcs@sha256:e44efa1effd03df21f14ce1d7ca586276bee053d64e3420897ae9d4fdf11d1a9"
        disks: "local-disk " + diskSize + " HDD"
        memory: memoryGb + "GB"
        mem: memoryGb + "G"
        cpu : threads
        cpus: runRequestThreads
        runtime_minutes: "600"
        maxRetries : 2
    }
    parameter_meta {
        finalBam: {
            description: "Input BAM filename",
            category: "required",
            localization_optional: true
        }
    }
}

task SplitVcf {
    input {
        File vcf
        String vcfPath = "~{vcf}"
        String prefix
        Array[String] splitVcfPaths
        Int maxRows = 1000
        Int minSplits = 2
        Int maxSplits = 10
        Int diskSize
        Int memoryGb = 1
    }

    command <<<
        set -e -o pipefail
        prefix=~{prefix}
        maxRows=~{maxRows}
        maxSplits=~{maxSplits}
        minSplits=~{minSplits}
        vcfPath=~{vcfPath}
        vcf=~{vcf}

        rm -f ${prefix}*.vcf

        filename=$( basename  ${vcfPath} )
        extension="${filename##*.}"
        filename="${filename%.*}"

        if [[ $extension == gz ]]; then
            input_path=${filename}

            gunzip -c \
            ${vcf} \
            > $input_path
        else
            input_path=${vcf}
        fi

        python /split_vcf.py \
        --vcf ${input_path} \
        --output-prefix ${prefix} \
        --max-rows ${maxRows} \
        --max-splits ${maxSplits} \
        --min-splits ${minSplits}
    >>>

    output {
        Array[File?] splitVcfs = splitVcfPaths
    }

    runtime {
        mem: memoryGb + "G"
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }
}

task CompressIndexVcf {
    input {
        File vcf
        String vcfCompressedPath = sub(basename(vcf), ".vcf$", ".vcf.gz")
        String vcfIndexedPath = sub(basename(vcf), ".vcf$", ".vcf.gz.tbi")
        Int memoryGb = 8
        Int diskSize = (ceil( size(vcf, "GB") )  * 2 ) + 4
    }

    command {
        set -e -o pipefail

        bgzip -c \
        ~{vcf} \
        > ~{vcfCompressedPath}

        tabix \
        -p vcf \
        ~{vcfCompressedPath}
    }

    output {
        IndexedVcf vcfCompressedIndexed = object {
            vcf : "~{vcfCompressedPath}",
            index : "~{vcfIndexedPath}"
        }
    }

    runtime {
        mem: memoryGb + "G"
        memory : memoryGb + "GB"
        disks: "local-disk " + diskSize + " HDD"
        docker: "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e"
        runtime_minutes: "1440"
    }

    meta {
        internalOnly : "False, produces output for external or internal BED files"
    }
}

task Download {
    input {
        String filenamePath
        String downloadUri
        File awsConfig
        File awsCredentials
        String endpointUrl
        # resources
        Int threads = 4
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize = 8
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    command <<<
        set -e -o pipefail
        awsConfig="~{awsConfig}"
        awsCredentials="~{awsCredentials}"
        filenamePath="~{filenamePath}"
        endpointUrl="~{endpointUrl}"
        downloadUri="~{downloadUri}"
        mkdir ~/.aws/
        cp \
        ${awsConfig} \
        ${awsCredentials} \
        ~/.aws/

        aws s3 \
        --endpoint-url ${endpointUrl} \
        cp \
        ${downloadUri} \
        .
    >>>

    output {
        File download = filenamePath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/awscli@sha256:67dc9052c4c286cafbfd9b93d6189e4f530a645a83f3f5fd96bd87563213127b"
        runtime_minutes: "6000"
        cpuPlatform : cpuPlatform
        partition: "cpu"
        qos: qos
    }
}

task ConcateTables {
    input {
        String outputTablePath
        Array[File] tables
        Int diskSize = 1
        Int memoryGb = 1
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

task Gatk4MergeSortVcf {
    input {
        String sortedVcfPath
        Array[File] tempVcfs
        IndexedReference referenceFa
        Boolean gzipped = false
        String suffix = if gzipped then ".tbi" else ".idx"
        Int threads = 1
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize = 10
    }

    Int jvmHeap = memoryGb * 750  # Heap size in Megabytes. mem is in GB. (75% of mem)

    command {
        export LANG=en_US.UTF-8

        gatk \
        SortVcf \
        --java-options "-Xmx~{jvmHeap}m -XX:ParallelGCThreads=4" \
        -SD ~{referenceFa.dict} \
        -I ~{sep=" -I " tempVcfs} \
        -O ~{sortedVcfPath}
    }

    output {
        IndexedVcf sortedVcf = object {
            vcf : "~{sortedVcfPath}",
            index : "~{sortedVcfPath}~{suffix}"
        }
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "gcr.io/nygc-public/broadinstitute/gatk4@sha256:b3bde7bc74ab00ddce342bd511a9797007aaf3d22b9cfd7b52f416c893c3774c"
        runtime_minutes: "500"
    }
}



task Extraction {
    input {
        Bram bram
        String sampleId
        String projectId
        IndexedVcf vcf
        IndexedReference referenceFa
        String extractedFeaturesPath = "~{sampleId}_extracted_features.csv"
        # resources
        Int threads = 1
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"

    }
    command <<<
        set -e -o pipefail
        sampleId="~{sampleId}"
        bram="~{bram.bram}"
        projectId="~{projectId}"
        vcf="~{vcf.vcf}"
        threads="~{threads}"
        referenceFa="~{referenceFa.fasta}"
        # Extract features for the EBM model. This script also extracts the reference sequence context for each variant, which is used in the mutational signature analysis.
        python3 /opt/fifa/src/cli.py \
            extract \
            -s ${sampleId} \
            -c ${projectId} \
            -v ${vcf} \
            -b ${bram} \
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "6000"
        cpuPlatform : cpuPlatform
        partition: "cpu"
        qos: qos
    }
}

task Prediction {
    input {
        String sampleId
        IndexedVcf vcf
        File extractedFeatures
        String fifaVcfPath = sub(basename(vcf.vcf, ".gz"), ".vcf$", ".fifa.vcf")
        Array[File] models
        # resources
        Int threads = 4
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize = 20
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }

    command <<<
        set -e -o pipefail
        sampleId="~{sampleId}"
        vcf="~{vcf.vcf}"
        extractedFeatures="~{extractedFeatures}"
        all_models="~{sep=' ' models}"
        python3 /opt/fifa/src/cli.py \
            predict \
            -s ${sampleId} \
            -v ${vcf} \
            -f ${extractedFeatures} \
            -o . \
            -m ${all_models}
    >>>

    output {
        File fifaVcf= fifaVcfPath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "300"
        partition: "cpu"
        cpuPlatform : cpuPlatform
        qos: qos
    }
}

task PredictionWithRna {
    input {
        String sampleId
        IndexedVcf vcf
        File extractedFeatures
        String fifaVcfPath = sub(basename(vcf.vcf, ".gz"), ".vcf$", ".fifa.vcf")
        Array[File] models
        File? optionalRnaFile
        # resources
        Int threads = 4
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize = 20
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }

    command <<<
        set -e -o pipefail
        sampleId="~{sampleId}"
        vcf="~{vcf.vcf}"
        extractedFeatures="~{extractedFeatures}"
        all_models="~{sep=' ' models}"
        all_rna="~{sep=' ' optionalRnaFile}"
        python3 /opt/fifa/src/cli.py \
            predict \
            -s ${sampleId} \
            -v ${vcf} \
            -f ${extractedFeatures} \
            -o . \
            -m ${all_models} \
            -r ${all_rna}
    >>>

    output {
        File fifaVcf= fifaVcfPath
    }

    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "6000"
        cpuPlatform : cpuPlatform
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
        String cpuPlatform = "Intel Cascade Lake"

    }
    command <<<
        set -e -o pipefail
        mergedModelId="~{mergedModelId}"
        modelPath="~{modelPath}"
        all_pickles="~{sep=' ' pickles}"
        python3 /opt/fifa/src/cli.py \
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
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "300"
        partition: "cpu"
        cpuPlatform : cpuPlatform
        qos: qos
    }
}

task ReTraining {
    input {
        Array[File] extractedFeatures
        File labels
        String modelPath = "fifa_model.pkl"
        String hyperParameterFlag
        # resources
        Int threads = 1
        Int runRequestThreads =  ceil(threads / 2.0)
        Int memoryGb = 4
        Int diskSize = 20
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    
    command <<<
        set -e -o pipefail
        extractedFeatures="~{sep=' ' extractedFeatures}"
        modelPath="~{modelPath}"
        labels="~{labels}"
        hyperParameterFlag="~{hyperParameterFlag}"

        mkdir extractedFeatures/
        ln -s \
        ${extractedFeatures} \
        extractedFeatures/
        python3 /opt/fifa/src/cli.py \
            retrain \
            -d extractedFeatures/ \
            -o ${modelPath} \
            --labels_path ${labels} ${hyperParameterFlag}
    >>>

    output {
        File model = modelPath
    }
    runtime {
        mem: memoryGb + "G"
        cpus: runRequestThreads
        cpu : threads
        disks: "local-disk " + diskSize + " HDD"
        memory : memoryGb + "GB"
        docker : "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:cc013f61e319445203d727bf7efa9d1ddc0f658eb65e7e74f335db45e79619de"
        runtime_minutes: "300"
        partition: "cpu"
        cpuPlatform : cpuPlatform
        qos: qos
    }
}