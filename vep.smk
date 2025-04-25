configfile: "configfile.yaml"
input_dir = config["input_dir"]
species = config["species"]
reference = config["ref"]
vcf = config["vcf"]

rule all:
    input:
        input_dir + f"results/{species}_annotated.vcf"


rule vep_annotate:
    message:
        """Annotating input vcf file"""
    input:
        vcf = input_dir + "/" + vcf,
        fasta = input_dir + "/" + fasta
        gtf = input_dir + "/" + gtf
    output:
         input_dir + f"/results/{species}_annotated.vcf"
    params:
        species = f"{species}",
        singularity = input_dir + "/vep.sif",
    resources:
        mem="50G",
        time="00:45:00",
        cpus=30 
    threads: 30
    shell:
        """
        apptainer run -B {input_dir}:{input_dir} \
            {params.singularity} vep \
            --input_file {input.vcf} \
            --fasta {input.fasta} \
            --gtf {input.gtf} \
            --output_file {output}  \
            --species {params.species} \
            --vcf \
            --fork {threads} \
            --force_overwrite
        """

