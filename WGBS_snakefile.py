configfile: "config.yaml"

rule all:
        input:
            'genome_folder/Bisulfite_Genome',
            'lambda_genome_folder/Bisulfite_Genome',
            expand('lambda.{sample}_trim_bismark_bt2.bam', sample=config["samples"]),
            expand('{sample}.CpG_report.merged_CpG_evidence.cov', sample=config["samples"])

rule trimming:
    input:
        "{sample}.fastq"
    output:
        "{sample}_trim.fastq"
    shell:
        """
        module load trimmomatic/0.36
        trimmomatic SE -threads 1 {input} {output} \
        ILLUMINACLIP:illumina_adapters.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:20 HEADCROP:5
        """

rule genome_prep:
    input:
        "genome_folder"
    output:
        "genome_folder/Bisulfite_Genome"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
        /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation {input}
        """

rule lambda_genome_prep:
    input:
        "lambda_genome_folder"
    output:
        "lambda_genome_folder/Bisulfite_Genome"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
        /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_genome_preparation {input}
        """

rule alignment:
    input:
        fastq="{sample}_trim.fastq",
        genome_dir="genome_folder/Bisulfite_Genome" # This is required to make sure previous rules have run
    output:
        "{sample}_trim_bismark_bt2.bam"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
        /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
        --multicore 3 --genome genome_folder --single_end {input.fastq}
        """

rule lambda_alignment:
    input:
        fastq="{sample}_trim.fastq",
        genome_dir="lambda_genome_folder/Bisulfite_Genome" # This is required to make sure previous rules have run
    output:
        "lambda.{sample}_trim_bismark_bt2.bam"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
        /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark \
        --multicore 3 --genome lambda_genome_folder --prefix lambda --single_end {input.fastq}
        """

rule deduplication:
    input:
        "{sample}_trim_bismark_bt2.bam"
    output:
        "{sample}_trim_bismark_bt2.deduplicated.bam"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
       /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/deduplicate_bismark {input}
        """

rule meth_extraction:
    input:
        "{sample}_trim_bismark_bt2.deduplicated.bam"
    output:
        "{sample}_trim_bismark_bt2.deduplicated.bismark.cov.gz"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
        /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/bismark_methylation_extractor \
        --single-end --no_overlap --comprehensive --bedgraph --report --cytosine_report \
        --genome_folder genome_folder {input}
        """

rule merge_cpgs:
    input:
        "{sample}_trim_bismark_bt2.deduplicated.bismark.cov.gz"
    output:
        "{sample}.CpG_report.merged_CpG_evidence.cov"
    params:
        "{sample}"
    shell:
        """
        module load bowtie2/2.3.5.1
        module load samtools/1.9 
        /scratch/monoallelic/hjm32/bin/Bismark-0.22.3/coverage2cytosine \
        -o {params} --merge_CpGs \
        --genome_folder genome_folder {input}
        """