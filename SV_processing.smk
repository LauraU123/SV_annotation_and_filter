# this workflow is to test the inputs and outputs of the variant calling file in sheep
configfile: "config_SV.yaml"
OUTPUT_DIR = config["outputdir"]
INPUT_DIR = config["inputdir"]
vcf_name = config["vcf"]
mode = config["mode"]
chrs = config["chr"]
species = config["species"]

rule all:
    input:
        in_ = expand(OUTPUT_DIR + f"{species}_{chrs}_unique_case_SVs.vcf.gz", species=species, chrs=chrs)

rule extract_chromosome:
    message:
        """Extracting chromosome {chrs}"""
    input:
        INPUT_DIR + vcf_name
    output:
        OUTPUT_DIR + f"chromosome_{chrs}.sv.vcf.gz"
    resources:
        mem="50G",
        time="00:45:00",
        cpus=30 
    params:
        chromosome = f"{chrs}",
    threads: 30
    shell:
        """
        module load BCFtools
        bcftools view -i '(CHROM="{params.chromosome}")'  {input} | grep -v  '^##contig' | bgzip >  {output}
        tabix -p vcf {output}
        """   

rule filter_svs:
    message:
        """Filtering the structural variants in chromosome {chrs} based on quality..."""
    input: rules.extract_chromosome.output
    output: OUTPUT_DIR + "{chrs}_filtered.sv.vcf.gz"
    resources:
        mem="50G",
        time="00:30:00",
        cpus=30 
    threads: 30
    shell:
        """
        module load BCFtools
        bcftools filter -i 'INFO/SVLEN>1000 && QUAL>30' {input} > {output}
        """
        
rule get_samples:
    input:
        vcf = rules.filter_svs.output,
        controls = "controls/{species}_chr{chrs}.txt"
    output:
        cases = "temp/{species}_{chrs}/cases.txt",
        ctrls = "temp/{species}_{chrs}/controls_verified.txt"
    resources:
        mem="20G",
        time="00:30:00",
        cpus=30 
    params:
        all_samples = "temp/{species}_{chrs}/all_samples.txt"
    shell:
        """
        module load BCFtools
        bcftools query -l {input.vcf} > {params.all_samples}
        
        # Verify controls exist in VCF
        grep -wFf {input.controls} {params.all_samples} > {output.ctrls}
        
        # Cases = All samples not in controls
        grep -v -wFf {output.ctrls} {params.all_samples} > {output.cases}
        """

# Step 2: Filter variants (1/1 in all cases, not 1/1 in any control)
rule filter_variants_cases:
    message:
        """Filtering based on case genotypes..."""
    input:
        vcf = rules.filter_svs.output,
        cases = rules.get_samples.output.cases,
    resources:
        mem="10G",
        time="00:30:00",
        cpus=30 
    output:
        out_ = OUTPUT_DIR + "{species}_{chrs}_cases.vcf.gz",
        keepsites  = "temp/{species}_{chrs}/keep_sites.txt"
    params:
        mode = mode
    shell:
        """
        module load BCFtools
        # if you are filtering strictly, only sites where cases are either 1/1 or 1|1 are kept.
        if [ "{params.mode}" = "recessive" ]; then
            bcftools query \
                -f '%CHROM\\t%POS[\\t%GT]\\n' \
                -s $(paste -sd, {input.cases}) \
                {input.vcf} | \
            awk '
                {{
                    pass=1;
                    for(i=3; i<=NF; i++) {{
                        if($i != "1/1" && $i != "1|1") {{
                            pass=0;
                            break;
                        }}
                    }}
                    if(pass) print $1"\\t"$2;
                }}
            ' > {output.keepsites}
        # if you're not using strict mode, both homozygous and heterozygous alternative genotypes will be kept in cases. 
        else
            bcftools query \
                -f '%CHROM\\t%POS[\\t%GT]\\n' \
                -s $(paste -sd, {input.cases}) \
                {input.vcf} | \
            awk '{{ 
                pass=1;
                for(i=3; i<=NF; i++) {{
                    if($i != "1/1" && $i != "1|1" && $i != "0/1" && $i != "0|1" && $i != "1/0" && $i != "1|0") {{
                        pass=0;
                        break;
                    }}
                }}
                if(pass) print $1"\\t"$2;
            }}' > {output.keepsites}
        fi
        #only locations of interest are kept. 
        bcftools view \
            -T {output.keepsites} \
            {input.vcf} -o {output.out_}

        tabix -p vcf {output.out_}
        """
rule filter_variants_controls:
    message:
        """Filtering based on control genotypes..."""
    input:
        filteredinput =  rules.filter_variants_cases.output.out_,
        controls = rules.get_samples.output.ctrls
    resources:
        mem="10G",
        time="00:30:00",
        cpus=30 
    output:
        out_ = OUTPUT_DIR + "{species}_{chrs}_unique_case_SVs.vcf.gz",
    params:
        excludesites  =  "temp/{species}_{chrs}/exclude_sites.txt",
        mode = mode
    shell:
        """
        module load BCFtools
        # Getting variants to exclude where any control is 1/1 or 1|1
        if [ "{params.mode}" = "recessive" ]; then
            bcftools query \
                -f '%CHROM\\t%POS[\\t%GT]\\n' \
                -s $(paste -sd, {input.controls}) \
                {input.filteredinput} | \
            awk '{{ 
                for(i=3; i<=NF; i++) {{
                    if($i == "1/1" || $i == "1|1") {{
                        print $1"\\t"$2;
                        break;
                    }}
                }}
            }}' > {params.excludesites}
        # if the mode is set for dominant genotypes, 0/1 are also excluded from the controls
        else
            bcftools query \
                -f '%CHROM\\t%POS[\\t%GT]\\n' \
                -s $(paste -sd, {input.controls}) \
                {input.filteredinput} | \
            awk '{{ 
                for(i=3; i<=NF; i++) {{
                    if($i == "1/1" || $i == "1|1" || $i == "0/1" || $i == "0|1" || $i == "1/0" || $i == "1|0") {{
                        print $1"\\t"$2;
                        break;
                    }}
                }}
            }}' > {params.excludesites}
        fi
        # these are excluded here. 
        bcftools view \
            -T ^{params.excludesites} \
            {input.filteredinput} -o {output.out_}
        tabix -p vcf {output.out_}
        """