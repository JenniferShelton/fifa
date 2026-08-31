from pathlib import Path

configfile: "extractastic_config.yaml"

SAMPLE = config["sample_id"]
PROJECT = config["project_id"]
BAM = config["bam"]
BAM_INDEX = config["bam_index"]
VCF = config["vcf"]
REFERENCE = config["reference"]
REFERENCE_INDEX = config["reference_index"]
THREADS = config.get("threads", 1)
MAX_SPLITS = 20
FIFA_CONTAINER = config.get(
    "fifa_container",
    "us.gcr.io/nygc-comp-s-fd4e/fifa@sha256:c48512a22d04097750e29572e77101f48a5435a2e5a944f8f01a8bb2da1797e3",
)
HTSLIB_CONTAINER = config.get(
    "htslib_container",
    "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e",
)
SPLIT_CONTAINER = config.get(
    "split_container",
    "gcr.io/nygc-comp-s-fd4e/gcs_htslib_suite@sha256:47eba58683641905a58f31c05605c953dc1d888466e8d434a6ceff47b76df03e",
)
CONCAT_CONTAINER = config.get(
    "concat_container",
    "gcr.io/nygc-public/somatic_dna_tools@sha256:f281e73ddf515f3a5db6e766ef12fb2331dafa2b27c88e426764dcd8b4a7e19b",
)


def extracted_feature_files(wildcards):
    split_directory = checkpoints.split_vcf.get(**wildcards).output.directory
    split_ids = glob_wildcards(
        str(Path(split_directory) / "{split}.vcf")
    ).split
    return expand(
        f"work/extraction/{{split}}/{SAMPLE}_extracted_features.csv",
        split=split_ids,
    )


rule all:
    input:
        f"results/{SAMPLE}_extracted_features.csv"


rule mobster_fit:
    input:
        vcf=VCF
    output:
        rds=f"work/{SAMPLE}.mobster_fit.rds"
    container:
        FIFA_CONTAINER
    shell:
        """
        Rscript /opt/fifa/src/run_mobster_fit.R {SAMPLE} {input.vcf} {output.rds}
        """


checkpoint split_vcf:
    input:
        vcf=VCF
    output:
        directory="work/splits"
    params:
        prefix="work/splits/",
        split_vcf_script=config.get("split_vcf_script", "/split_vcf.py")
    container:
        SPLIT_CONTAINER
    shell:
        """
        set -euo pipefail
        mkdir -p {output.directory}
        input_vcf={input.vcf}
        if [[ "$input_vcf" == *.vcf.gz || "$input_vcf" == *.vcf.bgz ]]; then
            input_vcf=work/{SAMPLE}.input.vcf
            gunzip -c {input.vcf} > "$input_vcf"
        fi
        python {params.split_vcf_script} \
            --vcf "$input_vcf" \
            --output-prefix {params.prefix} \
            --max-rows 1000 \
            --max-splits {MAX_SPLITS} \
            --min-splits 2
        """


rule compress_index_vcf:
    input:
        vcf="work/splits/{split}.vcf"
    output:
        vcf="work/splits/{split}.vcf.gz",
        index="work/splits/{split}.vcf.gz.tbi"
    container:
        HTSLIB_CONTAINER
    shell:
        """
        set -euo pipefail
        bgzip -c {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_features:
    input:
        bam=BAM,
        bam_index=BAM_INDEX,
        vcf="work/splits/{split}.vcf.gz",
        vcf_index="work/splits/{split}.vcf.gz.tbi",
        reference=REFERENCE,
        reference_index=REFERENCE_INDEX,
        mobster_fit=f"work/{SAMPLE}.mobster_fit.rds"
    output:
        features=f"work/extraction/{{split}}/{SAMPLE}_extracted_features.csv"
    threads:
        THREADS
    params:
        vcf_link=lambda wildcards: f"{wildcards.split}.vcf.gz",
        reference_link=lambda wildcards: Path(REFERENCE).name,
        bam_link=lambda wildcards: Path(BAM).name,
        output_directory=lambda wildcards: f"work/extraction/{wildcards.split}"
    container:
        FIFA_CONTAINER
    shell:
        """
        set -euo pipefail
        ln -sf {input.vcf} {params.vcf_link}
        ln -sf {input.vcf_index} {params.vcf_link}.tbi
        ln -sf {input.bam} {params.bam_link}
        ln -sf {input.bam_index} {params.bam_link}.bai
        ln -sf {input.reference} {params.reference_link}
        ln -sf {input.reference_index} {params.reference_link}.fai
        mkdir -p {params.output_directory}
        python3 /opt/fifa/src/cli.py extract \
            -s {SAMPLE} \
            -c {PROJECT} \
            -v {params.vcf_link} \
            -b {params.bam_link} \
            -r {params.reference_link} \
            -o {params.output_directory} \
            -n {threads} \
            --mobster-fit-rds {input.mobster_fit}
        """


rule concatenate_features:
    input:
        extracted_feature_files
    output:
        f"results/{SAMPLE}_extracted_features.csv"
    params:
        concatenate_script=config.get("concatenate_script", "/concate_tables.py")
    container:
        CONCAT_CONTAINER
    shell:
        """
        mkdir -p results
        python3 {params.concatenate_script} --tables {input} --output {output}
        """
