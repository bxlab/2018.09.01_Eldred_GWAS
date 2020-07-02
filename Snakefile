import glob

import pandas

# Read sample sheet csv. Looks like:
# ```
# Project,FCID,Lane,Index,SM_Tag,File_Name
# keldred1_GWAS,CCEPGANXX,1,GAGATTCC~TAATCTTA,12-1091_1,CCEPGANXX_1_GAGATTCC~TAATCTTA
# ```
sample_dict = {}
sample_sheet = pandas.read_csv( "sample_barcodes.csv", header=0 )
for index, row in sample_sheet.iterrows():
    sample_id = row['SM_Tag'].split("_")[0]
    fname_base = row['File_Name']
    sample_dict[sample_id] = fname_base

SAMPLES = list( sample_dict.keys() )
assert len(SAMPLES) == 738, "Expected 738 samples"
print( "Sample base names:", SAMPLES[:10], "..." )

# Targeted regions
REGIONS = [ line.strip().split()[:3] for line in open("regions.bed") ]
REGIONS_CONCAT = [ ( r[0] + ":" + r[1] + "-" + r[2] ) for r in REGIONS ]
print( "Regions:", REGIONS_CONCAT )
DATASETS = [ "gatk", "freebayes" ]

singularity: "docker://continuumio/miniconda3:4.4.10"

rule all:
    input:
        expand( "plots/freebayes${ethnicity}_coneratios_manhattan.pdf", ethnicity=['', 'AA', 'Asian', 'Cau']),
        "plots/freebayes_coneratios_genotypes.pdf",
        "plots/freebayes_coneratios_ethnicity.pdf",
        expand("stats/{sample}.isize.txt", sample=SAMPLES),
        expand("stats/{sample}.alignment_summary.txt", sample=SAMPLES),
        expand( "qc/qualimap/{sample}.bamqc", sample=SAMPLES )

rule multiqc_summary:
    input:
        "qc/multiqc.html"


# ---- Mapping ----

rule get_genome:
    output:
        "genome/{genome}.fa"
    params:
        genome=lambda wildcards: {wildcards.genome}
    shell:
        "curl http://hgdownload.soe.ucsc.edu/goldenPath/{params.genome}/bigZips/{params.genome}.fa.gz | gzip -d > {output}"

rule bwa_index:
    input:
        "genome/{genome}.fa"
    output:
        "genome/{genome}.amb",
        "genome/{genome}.ann",
        "genome/{genome}.bwt",
        "genome/{genome}.pac",
        "genome/{genome}.sa"
    log:
        "logs/bwa_index/{genome}.log"
    params:
        prefix="{genome}",
        algorithm="bwtsw"
    wrapper:
        "0.27.1/bio/bwa/index"

rule create_gatk_dict:
    input:
        "genome/{genome}.fa"
    output:
        "genome/{genome}.dict"
    log:
        "logs/picard/{genome}.dict.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/picard/createsequencedictionary"

rule create_fasta_index:
    input:
        "genome/{genome}.fa"
    output:
        "genome/{genome}.fa.fai"
    log:
        "logs/samtools/{genome}.fasta_index.log"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input}"

rule fastqc:
    input:
        "fastq/{sample}_S{lane}_L001_R1_001.fastq.gz"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}.zip"
    params:
        ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"

rule bwa_mem:
    input:
        # rawdata/keldred1_GWAS/FASTQ/CCEPGANXX_5_TCCGCGAA~CAGGACGT_1.fastq.gz
        reads=lambda wildcards: [ "fastq/" + sample_dict[wildcards.sample] + "_1.fastq.gz",
                                  "fastq/" + sample_dict[wildcards.sample] + "_2.fastq.gz" ]
    output:
        "mapped/{sample}.bam" 
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="genome/hg38",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        sort="samtools",
        sort_order="coordinate",
        sample="{sample}"
    threads:
        16
    wrapper:
        "0.27.1/bio/bwa/mem"

rule mark_duplicates:
    input:
        "mapped/{sample}.bam"
    output:
        bam="markdup/{sample}.bam",
        metrics="markdup/{sample}.metrics.txt"
    log:
        "logs/picard_markdup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/picard/markduplicates"

rule samtools_index:
    input:
        "markdup/{sample}.bam"
    output:
        "markdup/{sample}.bam.bai"
    log:
        "logs/samtools/{sample}_index.log"
    wrapper:
        "0.27.1/bio/samtools/index"


# ---- Quality Control ----

rule insert_size:
    input:
        "mapped/{sample}.bam"
    output:
        txt="stats/{sample}.isize.txt",
        pdf="stats/{sample}.isize.pdf"
    log:
        "logs/picard/insert_size/{sample}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/picard/collectinsertsizemetrics"

rule alignment_summary:
    input:
        ref="genome/hg38.fa",
        bam="markdup/{sample}.bam"
    output:
        "stats/{sample}.alignment_summary.txt"
    log:
        "logs/picard/alignment-summary/{sample}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/picard/collectalignmentsummarymetrics"

rule qualimap:
    input:
        bam="markdup/{sample}.bam",
        regions="regions.bed"
    output:
        directory("qc/qualimap/{sample}.bamqc")
    log:
        "logs/qualimap/{sample}.log"
    conda:
        "envs/qualimap.yaml"
    threads: 16
    shell:
        "qualimap bamqc -bam {input.bam} --feature-file {input.regions} --collect-overlap-pairs -nt {threads} --outside-stats -outdir {output}"

rule multiqc:
    input:
        "qc/",
        "stats/",
        "markdup/"
    output:
        "qc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/multiqc"


# ---- dbSNP ----

rule get_feature_file:
    output:
        "Homo_sapiens_assembly38.dbsnp138.vcf"
    log:
        "logs/download_feature_file.log"
    conda:
        "envs/gsutil.yaml"
    shell:
        "gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/{output} ."

rule get_variation_set:
    output:
        "All_20180418.vcf.gz"
    log:
        "logs/download_variation_sets.log"
    shell:
        """
        wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/{output}; \
        wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/{output}.tbi
        """

rule create_dbsnp_feature_index:
    input:
        "Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        "Homo_sapiens_assembly38.dbsnp138.vcf.idx"
    log:
        "logs/gatk/index_feature_file.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/indexfeaturefile"


# ---- GATK ----

rule gatk_bqsr:
    input:
        bam="markdup/{sample}.bam",
        ref="genome/hg38.fa",
        known="Homo_sapiens_assembly38.dbsnp138.vcf",
        refdict="genome/hg38.fa.fai",
        knownidx="Homo_sapiens_assembly38.dbsnp138.vcf.idx"
    output:
        bam="recal/{sample}.bam"
    log:
        "logs/gatk/bqsr/{sample}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/baserecalibrator"

rule gatk_call_variants:
    input:
        bam="recal/{sample}.bam",
        ref="genome/hg38.fa",
        known="Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        gvcf="gatk_called/{sample}.{region}.g.vcf.gz"
        #gvcf=temp("gatk_called/{sample}.{region}.g.vcf.gz")
    params:
        extra=lambda wildcards: f"--intervals {wildcards.region} --ploidy {'1' if 'chrX' in wildcards.region else '2'}"
    log:
        "logs/gatk/haplotypecaller/{sample}.{region}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/haplotypecaller"

rule gatk_combine_calls:
    input:
        ref="genome/hg38.fa",
        gvcfs=expand("gatk_called/{sample}.{{region}}.g.vcf.gz", sample=SAMPLES),
        index=expand("gatk_called/{sample}.{{region}}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        gvcf="gatk_called_combined/all.{region}.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.{region}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/combinegvcfs"

rule create_feature_index:
    input:
        "{folder}/{feature}.vcf.gz"
    output:
        "{folder}/{feature}.vcf.gz.tbi"
    log:
        "logs/gatk/{folder}_index_{feature}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/indexfeaturefile"

rule gatk_genotype_variants:
    input:
        ref="genome/hg38.fa",
        gvcf="gatk_called_combined/all.{region}.g.vcf.gz",
        index="gatk_called_combined/all.{region}.g.vcf.gz.tbi"
    output:
        vcf="gatk_genotyped/all.{region}.vcf.gz"
    params:
        extra=lambda wildcards: f"--use-new-qual-calculator --intervals {wildcards.region} {'--ploidy 1' if 'chrX' in wildcards.region else '--ploidy 2'}"
    log:
        "logs/gatk/genotypegvcfs.{region}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/genotypegvcfs"


rule gatk_merge_variants:
    input:
        vcfs=expand("gatk_genotyped/all.{region}.vcf.gz", region=REGIONS_CONCAT),
        index=expand("gatk_genotyped/all.{region}.vcf.gz.tbi", region=REGIONS_CONCAT)
    output:
        vcf="gatk_genotyped/all.vcf.gz"
    log:
        "logs/gatk/picard-merge-genotyped.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/picard/mergevcfs"

# # hard filtering as outlined in GATK docs
# # (https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)
# snvs:
#   "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
# indels:
#   "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")

rule gatk_select_calls:
    input:
        ref="genome/hg38.fa",
        vcf="gatk_genotyped/all.vcf.gz"
    output:
        vcf="gatk_filtered/all.{vartype,\w+}.vcf.gz"
        #vcf=temp("gatk_filtered/all.{vartype,\w+}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{vartype}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/selectvariants"

def get_gatk_filter(wildcards):
    if wildcards.vartype == "snvs":
        return { "snv-hard-filter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" }
    elif wildcards.vartype == "indels":
        return { "snv-hard-filter": "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" }
    else:
        assert False, "Unknown variant type: " + wildcards.vartype

rule gatk_hard_filter_calls:
    input:
        ref="genome/hg38.fa",
        vcf="gatk_filtered/all.{vartype}.vcf.gz"
    output:
        vcf="gatk_filtered/all.{vartype}.hardfiltered.vcf.gz"
        #vcf=temp("gatk_filtered/all.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_gatk_filter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/gatk/variantfiltration"

rule gatk_merge_calls:
    input:
        vcfs=expand("gatk_filtered/all.{vartype}.{filtertype}.vcf.gz",
                    vartype=["snvs", "indels"],
                    filtertype="hardfiltered"),
        indices=expand("gatk_filtered/all.{vartype}.{filtertype}.vcf.gz.tbi",
                       vartype=["snvs", "indels"],
                       filtertype="hardfiltered"),
    output:
        vcf="gatk_filtered/all.vcf.gz"
    log:
        "logs/picard/merge-filtered.log"
    wrapper:
        "https://raw.githubusercontent.com/bxlab/snakemake-wrappers/0.60.1b/bio/picard/mergevcfs"

rule gatk_add_ids:
    input: 
        vcf="gatk_filtered/all.vcf.gz",
        dbsnp="All_20180418.vcf.gz"
    output:
        "gatk_annotated/all.vcf.gz"
    log:
        "logs/add_ids.log"
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        '''SnpSift annotate -a {input.dbsnp} {input.vcf} | bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" | bgzip > {output}'''


# ---- FreeBayes ----

rule freebayes_windows:
    input:
        "regions.bed"
    output:
        "freebayes/{region}.windows"
    log:
        "logs/freebayes/{region}.windows.log"
    conda:
        "envs/bedtools.yaml"
    shell:
        '''echo "{wildcards.region}" | tr ":-" "\t\t" | bedtools makewindows -b /dev/stdin -w 10000 | sed "s/\t/:/" | sed "s/\t/-/" > {output}'''

rule freebayes:
    input:
        ref="genome/hg38.fa",
        windows="freebayes/{region}.windows",
        samples=expand("markdup/{sample}.bam", sample=SAMPLES),
        indices=expand("markdup/{sample}.bam.bai", sample=SAMPLES)
    output:
        "freebayes/{region}.vcf"  
    log:
        "logs/freebayes/{region}.vcf.log"
    params:
        ploidy=lambda wc: "-p 1" if "chrX" in wc.region else ""
    conda:
        "envs/freebayes.yaml"
    threads:
        5
    shell:
       "freebayes-parallel {input.windows} {threads} {params.ploidy} -f {input.ref} --use-best-n-alleles 4 {input.samples} > {output}"

rule freebayes_cleanup_vcf:
    input:
        vcf="freebayes/{region}.vcf",
        dbsnp="All_20180418.vcf.gz"
    output:
        "freebayes_annotated/{region}.vcf.gz",
    params:
        temp="freebayes/{region}_temp"
    log:
        "logs/freebayes/annotate_{region}.vcf.log"
    wildcard_constraints:
        region='chr.+'
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        """
        cat {input.vcf} \
        | vcffilter -s -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
        | vcfallelicprimitives -m --keep-geno --keep-info > {params.temp}; \
        vcfsort {params.temp} \
        | SnpSift annotate -a {input.dbsnp} \
        | bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" \
        | sw/vcfremovedups \
        | vcfuniq \
        | bgzip > {output}; \
        rm {params.temp}
        """

rule freebayes_subset:
    input:
        vcf="freebayes_annotated/{region}.vcf.gz",
        pheno="ethnicity.csv"
    output:
        "freebayes{ethnicity}_annotated/{region}.vcf.gz"
    log:
        "logs/freebayes/{ethnicity}_{region}.vcf.log"
    wildcard_constraints:
        region='chr.+',
        ethnicity='.+'
    params:
        ethnicity="{ethnicity}"
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        """
        sw/filter_subset {input.pheno} {input.vcf} {params.ethnicity} \
        | bgzip > {output}
        """

rule freebayes_all_vcf:
    input:
        lambda wc: expand( "{folder}_annotated/{region}.vcf.gz", folder=wc.folder, region=REGIONS_CONCAT )
    output:
        "{folder}_annotated/all.vcf.gz"
    params:
        temp="{folder}_annotated/temp"
    log:
        "logs/{folder}/all.vcf.log"
    wildcard_constraints:
        folder="freebayes.*"
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        """
        zcat {input} | vcffirstheader > {params.temp}; \
        vcfsort {params.temp} | bgzip > {output};
        rm {params.temp}
        """


# ---- Plink ----

rule vcf_to_plink:
    input:
        vcf="{dataset}_annotated/all.vcf.gz",
        ph="coneratios.pheno"
    output:
        "plink/{dataset}.bed"
    log:
        "logs/plink/{dataset}.bed.log"
    params:
        base="{dataset}"
    conda:
        "envs/plink.yaml"
    shell:
        """
        awk '{{ print $1, $1, 1 }}' {input.ph} > plink/{params.base}.update_sex
        plink --vcf {input.vcf} --vcf-filter --make-bed --update-sex plink/{params.base}.update_sex --out plink/{params.base}
        """

rule plink_recode:
    input:
        "plink/{dataset}.bed"
    output:
        "plink/{dataset}_recoded.ped"
    log:
        "logs/plink/{dataset}_recoded.ped.log"
    params:
        base="{dataset}"
    conda:
        "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --recode 12 --maf 0.035 --geno 0.01 --out plink/{params.base}_recoded"

rule plink_ped_to_bim:
    input:
        "plink/{dataset}.ped"
    output:
        "plink/{dataset}.bim"
    log:
        "logs/plink/ped_to_bed.{dataset}.log"
    params:
        base="{dataset}"
    conda:
        "envs/plink.yaml"
    shell:
        "plink --file plink/{params.base} --make-bed --out plink/{params.base}"

rule plink_pca:
    input:
        "plink/{dataset}.bed"
    output:
        "plink/{dataset}.eigenval", "plink/{dataset}.eigenvec"
    log:
        "logs/plink/{dataset}.pca.log"
    params:
        base="{dataset}"
    conda:
        "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --pca --out plink/{params.base}"

rule plink_cluster:
    input:
        "plink/{dataset}.bed"
    output:
        "plink/{dataset}.mds"
    log:
        "logs/plink/{dataset}.mds.log"
    params:
        base="{dataset}"
    conda:
        "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --allow-extra-chr --cluster --mds-plot 3 --out plink/{params.base}"

rule plink_assoc:
    input:
        "plink/{dataset}.bed",
        ph="{pheno}.pheno"
    output:
        "plink/{dataset}_{pheno}.qassoc"
    log:
        "logs/plink/{dataset}_{pheno}.qassoc.log"
    params:
        base="{dataset}", pheno="{pheno}"
    conda:
        "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --allow-extra-chr --maf 0.01 --pheno {input.ph} --assoc --allow-no-sex --out plink/{params.base}_{params.pheno}"

rule plink_linear:
    input:
        "plink/{dataset}.bed",
        ph="{pheno}.pheno",
        eigenvec="plink/{dataset}.eigenvec"
    output:
        #"plink/{dataset}_{pheno}_eigen.assoc.linear"
        "plink/{dataset}_{pheno}.assoc.linear"
    log:
        "logs/plink/{dataset}_{pheno}.qassoc.log"
    params:
        base="{dataset}", pheno="{pheno}"
    threads:
        24
    conda:
        "envs/plink.yaml"
    shell:
        #"plink --bfile plink/{params.base} --allow-extra-chr --maf 0.05 --geno 0.1 --pheno {input.ph} --threads {threads} --linear mperm=1000 --covar {input.eigenvec} --out plink/{params.base}_{params.pheno}_eigen"
        "plink --bfile plink/{params.base} --allow-extra-chr --maf 0.05 --geno 0.1 --pheno {input.ph} --threads {threads} --linear mperm=1000 --out plink/{params.base}_{params.pheno}"

rule extract_pheno:
    input:
        "coneratios.xlsx"
    output:
        "coneratios.pheno"
    log:
        "logs/coneratios.pheno.log"
    conda:
        "envs/csv_stuff.yaml"
    shell:
        """in2csv {input} | csvcut -x -c "DNA ID4","DNA ID4","Cone Ratio" | csvformat -T | tail -n +2 > {output}"""


# ---- LDAK ----

rule get_ldak:
    output:
        "sw/ldak"
    log:
        "logs/ldak_download.log"
    shell:
        """
        wget http://dougspeed.com/wp-content/uploads/ldak5.1.linux_.fast_.zip
        gunzip -c ldak5.1.linux_.fast_.zip > {output}
        chmod a+rx {output}
        rm ldak5.1.linux_.fast_.zip
        """

rule ldak_prune:
    input:
        ldak="sw/ldak",
        bim="plink/{dataset}_recoded.bim"
    output:
        "ldak/{dataset}.in"
    log:
        "logs/ldak/{dataset}.prune.log"
    params:
        base="{dataset}"
    shell:
        "{input.ldak} --bfile plink/{params.base}_recoded --thin ldak/{params.base} --window-prune 0.75 --window-kb 100"

rule ldak_get_kinship:
    input:
        ldak="sw/ldak",
        bim="plink/{dataset}_recoded.bim",
        extract="ldak/{dataset}.in"
    output:
        "ldak/{dataset}.grm.bin"
    log:
        "logs/ldak/{dataset}.kinship.log"
    params:
        base="{dataset}"
    shell:
        "{input.ldak} --bfile plink/{params.base}_recoded --ignore-weights YES --power -1 --extract {input.extract} --calc-kins-direct ldak/{params.base}"

rule ldak_pca:
    input:
        ldak="sw/ldak",
        bim="plink/{dataset}_recoded.bim",
        extract="ldak/{dataset}.in",
        grm="ldak/{dataset}.grm.bin"
    output:
        "ldak/{dataset}.vect"
    log:
        "logs/ldak/{dataset}.pca.log"
    params:
        base="{dataset}"
    shell:
        "{input.ldak} --bfile plink/{params.base}_recoded --pca ldak/{params.base} --axes 5 --grm ldak/{params.base} --extract {input.extract}"

rule ldak_linear:
    input:
        ldak="sw/ldak",
        bim="plink/{dataset}_recoded.bim",
        extract="ldak/{dataset}.in",
        grm="ldak/{dataset}.grm.bin",
        pca="ldak/{dataset}.vect",
        pheno="{pheno}.pheno"
    output:
        "ldak/{dataset}_{pheno}.pvalues"
    log:
        "logs/ldak/{dataset}_{pheno}.linear.log"
    params:
        base="{dataset}",
        pheno="{pheno}"
    shell:
        "{input.ldak} --bfile plink/{params.base}_recoded --linear ldak/{params.base}_{params.pheno} --pheno {input.pheno} --grm ldak/{params.base} --covar {input.pca} --extract {input.extract}"

rule append_ldak_pvalues:
    input:
        bim="plink/{dataset}_recoded.bim",
        pvalues="ldak/{dataset}_{pheno}.pvalues"
    output:
        "ldak/{dataset}_{pheno}.pvalues.annotated"
    log:
        "log/ldak/{dataset}_{pheno}.annotate.log"
    shell:
        """
        tail -n +2 {input.pvalues} | awk 'BEGIN{{zero="0"; OFS="\t";}} {{print $1,zero}}' | sort -k1,1 > all_snps.txt
        cat {input.bim}  | grep -F -f all_snps.txt | sort -k2,2 > A
        tail -n +2 {input.pvalues} | sort -k1,1 | awk '{{print $2}}' > B
        echo -e "CHR\tSNP\tSEX\tBP\tA1\tA2\tP" > {output}
        paste A B >> {output}
        rm A B all_snps.txt
        """


# ---- Plots ----

rule manhattan_plot:
    input:
        region="regions.bed",
        pvalues="ldak/{dataset}_{pheno}.pvalues.annotated"
    output:
        pdf="plots/{dataset}_{pheno}_manhattan.pdf",
        snps="ldak/{dataset}_{pheno}.sig_snps"
    log:
        "log/matplotlib/{dataset}_{pheno}_manhattan.log"
    conda:
        "envs/matplotlib.yaml"
    shell:
        "sw/plot_manhattan {input.region} {input.pvalues} {output.snps} {output.pdf} 0.05"

rule genotype_plot:
    input:
        pheno="{pheno}.pheno",
        region="regions.bed",
        vcf="{dataset}_annotated/all.vcf.gz",
        snps="ldak/{dataset}_{pheno}.sig_snps"
    output:
        "plots/{dataset}_{pheno}_genotypes.pdf"
    log:
        "log/matplotlib/{dataset}_{pheno}_genotypes.log"
    conda:
        "envs/matplotlib.yaml"
    shell:
        "sw/plot_genotypes {input.pheno} {input.region} {input.vcf} {input.snps} {output}"

rule ethnicity_plot:
    input:
        region="regions.bed",
        vcf="{dataset}_annotated/all.vcf.gz",
        pheno="{pheno}.pheno",
        ethnicity="ethnicity.csv",
        snps="ldak/{dataset}_{pheno}.sig_snps"
    output:
        "plots/{dataset}_{pheno}_ethnicity.pdf"
    log:
        "log/matplotlib/{dataset}_{pheno}_ethnicity.log"
    conda:
        "envs/matplotlib.yaml"
    shell:
        "sw/plot_ethnic_genotypes {input.pheno} {input.ethnicity} {input.region} {input.vcf} {input.snps} {output}"



