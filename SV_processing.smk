
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
    output: 
        filtered_vcf = OUTPUT_DIR + f"{chrs}_filtered.sv.vcf.gz"
    resources:
        mem="50G",
        time="00:30:00",
        cpus=30 
    threads: 30
    params:
        exclude_samples= ",".join(config["exclude_samples"])
    shell:
        """
        module load BCFtools
        if [ -n "{params.exclude_samples}" ]; then
            echo "{params.exclude_samples}" | tr ',' '\\n' > temp/exclude_samples.txt
            echo "Excluding samples: $(cat temp/exclude_samples.txt)"
            bcftools view -S ^temp/exclude_samples.txt {input} > {output.filtered_vcf}
        else
            cp {input} {output.filtered_vcf}
        fi
        """

"""
        module load BCFtools
        if [ -n "{params.exclude_samples}" ]; then
            echo "{params.exclude_samples}" | tr ',' '\\n' > temp/exclude_samples.txt
            echo "Excluding samples: $(cat temp/exclude_samples.txt)"
            bcftools filter -i 'INFO/SVLEN>1000 && QUAL>30' {input} | bcftools view -S ^temp/exclude_samples.txt -o {output.filtered_vcf}
        else
            bcftools filter -i 'INFO/SVLEN>1000 && QUAL>30' {input} > {output.filtered_vcf}
        fi
"""
        
rule get_samples:
    message:
        """Finding the cases based on the controls..."""
    input:
        vcf = rules.filter_svs.output,
    output:
        cases = "temp/{species}/cases.txt",
        ctrls = "temp/{species}/controls_verified.txt"
    resources:
        mem="20G",
        time="00:30:00",
        cpus=30 
    params:
        all_samples = "temp/{species}/all_samples.txt",
        case_samples= ",".join(config["case_samples"]),
    shell:
        """
        module load BCFtools
        bcftools query -l {input.vcf} > {params.all_samples}
        
        # Verify cases exist in VCF
        echo "{params.case_samples}" | tr ',' '\\n' > temp/case_samples.txt
        grep -wFf temp/case_samples.txt {params.all_samples} > {output.cases}
        
        grep -v -wFf {output.cases} {params.all_samples} > {output.ctrls}
        """

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