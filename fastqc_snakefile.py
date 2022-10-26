configfile: "config.yaml"

rule all:
        input:
            expand('{sample}_fastqc.html', sample=config["samples"])

rule fastqc:
    input:
        "{sample}.fastq"
    output:
        html="{sample}_fastqc.html",
        zip="{sample}_fastqc.zip"
    shell: 
        """
        module load fastqc/0.11.5
        fastqc -t 15 {input}
        """