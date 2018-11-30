import glob

import pandas

# Read sample sheet csv. Looks like:
# ```
# Project,FCID,Lane,Index,SM_Tag,File_Name
# keldred1_GWAS,CCEPGANXX,1,GAGATTCC~TAATCTTA,12-1091_1,CCEPGANXX_1_GAGATTCC~TAATCTTA
# ```
sample_dict = {}
sample_sheet = pandas.read_csv( "rawdata/keldred1_GWAS/keldred1_GWAS.csv", header=0 )
for index, row in sample_sheet.iterrows():
    sample_id = row['SM_Tag'].split("_")[0]
    fname_base = row['File_Name']
    sample_dict[sample_id] = fname_base

SAMPLES = list( sample_dict.keys() )
assert len(SAMPLES) == 748, "Expected 748 samples"
print( "Sample base names:", SAMPLES[:10], "..." )

# Targetted regions
REGIONS = [ line.strip().split()[:3] for line in open("regions.bed") ]
REGIONS_CONCAT = [ ( r[0] + ":" + r[1] + "-" + r[2] ) for r in REGIONS ]
print( "Regions:", REGIONS_CONCAT )

singularity: "docker://continuumio/miniconda3:4.4.10"

rule all:
    input:
        expand( "recal/{sample}.bam", sample=SAMPLES ),
        expand( "stats/{sample}.isize.txt", sample=SAMPLES ),
        expand( "qc/qualimap/{sample}.bamqc", sample=SAMPLES ),
        expand( "gatk_called/{sample}.{region}.g.vcf.gz", sample=SAMPLES, region=REGIONS_CONCAT )

rule get_genome:
    output:
        "genome/hg38.fa"
    shell:
        "curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gzip -d > {output}"

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
        "genome/hg38.fa"
    output:
        "genome/hg38.dict"
    log:
        "logs/picard/hg38.dict.log"
    wrapper:
        "0.27.1/bio/picard/createsequencedictionary"

rule fastqc:
    input:
        "rawdata/{sample}_S{lane}_L001_R1_001.fastq.gz"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"

rule bwa_mem:
    input:
        # rawdata/keldred1_GWAS/FASTQ/CCEPGANXX_5_TCCGCGAA~CAGGACGT_1.fastq.gz
        reads=lambda wildcards: [ "rawdata/keldred1_GWAS/FASTQ/" + sample_dict[wildcards.sample] + "_1.fastq.gz",
                                  "rawdata/keldred1_GWAS/FASTQ/" + sample_dict[wildcards.sample] + "_2.fastq.gz" ]
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="/scratch/groups/jtayl139/cache/bwa_indices/hg38/hg38",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        sort="samtools",
        sort_order="coordinate"
    threads: 16
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
        "0.27.1/bio/picard/markduplicates"

rule gatk_bqsr:
    input:
        bam="markdup/{sample}.bam",
        ref="genome/hg38.fa",
        known="Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        bam="recal/{sample}.bam"
    log:
        "logs/gatk/bqsr/{sample}.log"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"

rule samtools_index:
    input: "markdup/{sample}.bam"
    output: "markdup/{sample}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"

rule insert_size:
    input:
        "mapped/{sample}.bam"
    output:
        txt="stats/{sample}.isize.txt",
        pdf="stats/{sample}.isize.pdf"
    log:
        "logs/picard/insert_size/{sample}.log"
    wrapper:
        "0.27.1/bio/picard/collectinsertsizemetrics"

rule alignment_summary:
    input:
        ref="/scratch/groups/jtayl139/cache/Fasta/hg38/hg38.fa",
        bam="markdup/{sample}.bam"
    output:
        "stats/{sample}.alignment_summary.txt"
    log:
        "logs/picard/alignment-summary/{sample}.log"
    wrapper:
        "0.27.1/bio/picard/collectalignmentsummarymetrics"

rule qualimap:
    input:
        bam="markdup/{sample}.bam",
        regions="regions.bed"
    output:
        directory("qc/qualimap/{sample}.bamqc")
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
        "https://bitbucket.org/james_taylor/snakemake-wrappers/raw/master/bio/multiqc"

rule gatk_call_variants:
    input:
        bam="recal/{sample}.bam",
        ref="genome/hg38.fa",
        known="Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        gvcf=temp("gatk_called/{sample}.{region}.g.vcf.gz")
    params:
        extra=lambda wildcards: f"--intervals {wildcards.region} --ploidy {'1' if 'chrX' in wildcards.region else '2'}"
    log:
        "logs/gatk/haplotypecaller/{sample}.{region}.log"
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"

rule gatk_combine_calls:
    input:
        ref="genome/hg38.fa",
        gvcfs=expand("gatk_called/{sample}.{{region}}.g.vcf.gz", sample=SAMPLES)
    output:
        gvcf="gatk_called_combined/all.{region}.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.{region}.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"

rule gatk_genotype_variants:
    input:
        ref="genome/hg38.fa",
        gvcf="gatk_called_combined/all.{region}.g.vcf.gz"
    output:
        vcf="gatk_genotyped/all.{region}.vcf.gz"
    params:
        extra=lambda wildcards: f"--use-new-qual-calculator --intervals {wildcards.region} {'--ploidy 1' if 'chrX' in wildcards.region else '--ploidy 2'}"
    log:
        "logs/gatk/genotypegvcfs.{region}.log"
    wrapper:
        "https://bitbucket.org/james_taylor/snakemake-wrappers/raw/master/bio/gatk/genotypegvcfs"


rule gatk_merge_variants:
    input:
        vcf=expand("gatk_genotyped/all.{region}.vcf.gz", region=REGIONS_CONCAT)
    output:
        vcf="gatk_genotyped/all.vcf.gz"
    log:
        "logs/gatk/picard-merge-genotyped.log"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"

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
        vcf=temp("gatk_filtered/all.{vartype,\w+}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/selectvariants"

def get_filter(wildcards):
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
        vcf=temp("gatk_filtered/all.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"

rule gatk_merge_calls:
    input:
        vcf=expand("gatk_filtered/all.{vartype}.{filtertype}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="hardfiltered")
    output:
        vcf="gatk_filtered/all.vcf.gz"
    log:
        "logs/picard/merge-filtered.log"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"

rule gatk_add_ids:
    input: "gatk_filtered/all.vcf.gz"
    output: "gatk_annotated/all.vcf.gz"
    log: "logs/add_ids.log"
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        '''SnpSift annotate -a All_20180418.vcf.gz {input} | bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" | bgzip > {output}'''

# ---- FreeBayes ----

rule freebayes_windows:
    input:
        "regions.bed"
    output:
        "freebayes/{region}.windows"
    conda:
        "envs/bedtools.yaml"
    shell:
        '''echo "{wildcards.region}" | tr ":-" "\t\t" | bedtools makewindows -b /dev/stdin -w 10000 | sed "s/\t/:/" | sed "s/\t/-/" > {output}'''


rule freebayes:
    input:
        ref="genome/hg38.fa",
        windows="freebayes/{region}.windows",
        samples=expand("markdup/{sample}.bam", sample=SAMPLES ),
        indices=expand("markdup/{sample}.bam.bai", sample=SAMPLES )
    output:
        "freebayes/{region}.vcf"  
    log:
        "logs/freebayes/{region}.vcf.log"
    params:
        ploidy=lambda wc: "-p 1" if "chrX" in wc.region else ""
    conda:
        "envs/freebayes.yaml"
    threads: 32
    shell:
       "freebayes-parallel {input.windows} {threads} {params.ploidy} -f {input.ref} --use-best-n-alleles 4 {input.samples} > {output}"

rule freebayes_all_vcf:
    input:
        expand( "freebayes/{region}.vcf", region=REGIONS_CONCAT )
    output:
        "freebayes/all.vcf"
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        "cat {input} | vcffirstheader | vcfstreamsort | vcfuniq > {output}"

rule freebayes_cleanup_vcf:
    input: "freebayes/all.vcf"
    output: "freebayes_annotated/all.vcf.gz"
    conda:
        "envs/vcf_stuff.yaml"
    shell:
        """
        cat {input} \
        | vcffilter -f "QUAL > 40" \
        | vcfallelicprimitives -m --keep-geno --keep-info \
        | SnpSift annotate -a All_20180418.vcf.gz \
        | bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" \
        | bgzip > {output}
        """

        # | vt normalize -r genome/hg38.fa - \

# ---- Plink ----

rule vcf_to_plink:
    input: "calls/{dataset}.vcf.gz"
    output: "plink/{dataset}.bed"
    params: base="{dataset}"
    conda: "envs/plink.yaml"
    shell:
        """
        awk '{{ print $1, $1, 1 }}' coneratios.pheno > plink/{params.base}.update_sex
        plink --vcf {input} --vcf-filter --make-bed --update-sex plink/{params.base}.update_sex --out plink/{params.base}
        """

rule pca:
    input: "plink/{dataset}.bed"
    output: "plink/{dataset}.eigenval", "plink/{dataset}.eigenvec"
    params: base="{dataset}"
    conda: "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --pca --out plink/{params.base}"

rule cluster:
    input: "plink/{dataset}.bed"
    output: "plink/{dataset}.mds"
    params: base="{dataset}"
    conda: "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --allow-extra-chr --cluster --mds-plot 3 --out plink/{params.base}"

rule assoc:
    input: "plink/{dataset}.bed", ph="{pheno}.pheno"
    output: "plink/{dataset}_{pheno}.qassoc"
    params: base="{dataset}", pheno="{pheno}"
    conda: "envs/plink.yaml"
    shell:
        "plink --bfile plink/{params.base} --allow-extra-chr --maf 0.01 --pheno {input.ph} --assoc --allow-no-sex --out plink/{params.base}_{params.pheno}"

rule linear:
    input: "plink/{dataset}.bed", ph="{pheno}.pheno", eigenvec="plink/{dataset}.eigenvec"
    output: "plink/{dataset}_{pheno}.assoc.linear"
    params: base="{dataset}", pheno="{pheno}"
    threads: 24
    conda: "envs/plink.yaml"
    shell:
        # "plink --bfile plink/{params.base} --allow-extra-chr --maf 0.05 --geno 0.1 --pheno {input.ph} --threads {threads} --linear mperm=1000 --covar {input.eigenvec} --out plink/{params.base}_{params.pheno}"
        "plink --bfile plink/{params.base} --allow-extra-chr --maf 0.05 --geno 0.1 --pheno {input.ph} --threads {threads} --linear mperm=1000 --out plink/{params.base}_{params.pheno}"

rule extract_pheno:
    input: "coneratios.xlsx"
    output: "coneratios.pheno"
    shell:
        """in2csv {input} | csvcut -x -c "DNA ID4","DNA ID4","Cone Ratio" | csvformat -T | tail -n +2 > {output}"""
