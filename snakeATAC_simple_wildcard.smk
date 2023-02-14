import os
import socket
import re

genome = config["genome"]
chrom_sizes = "~/genomes/" + genome + "/" + genome + ".chrom.sizes"
SUFFIX = ".fastq.gz"
FASTQ = "fastq"
TRIM_DIR = "trim"
BOWTIE2_DIR = "bowtie2"
METHYLQA_DIR = "methylQA"
BW_DIR = "bw"

threads_trim = 6
threads_bowtie2 = 6
sambamba_options = " --hash-table-size=3000000 --overflow-list-size=3000000 "

adapter_1="CTGTCTCTTATACACATCT"
adapter_2="CTGTCTCTTATACACATCT"

bowtie2_unal = "--no-unal"
if "unal" in config.keys():
    bowtie2_unal = ""

hostname = socket.gethostname()
if len(re.findall("compute1", hostname)) > 0:
    server = "ris"
elif len(re.findall("^n\\d+$", hostname)) == 0:
    server = "wanglab"
else:
    server = "htcf"

if server == "ris":
    Rscript = "/usr/bin/Rscript"
    cutadapt = "cutadapt"
    bowtie2 = "bowtie2"
    samtools = "samtools"
    methylQA = "~/software/methylQA_0.2.1/methylQA"
    bedtools = "bedtools"
    bedGraphToBigWig = "~/software/kent/bin/bedGraphToBigWig"
    sambamba = "~/software/sambamba/sambamba"
    bamCoverage = "~/anaconda3/envs/common/bin/bamCoverage"
    tabix = "tabix"
    bgzip = "bgzip"
else:
    Rscript = "/opt/apps/R/3.6.1/bin/Rscript"
    cutadapt = "/opt/apps/python3/bin/cutadapt"
    bowtie2 = "/opt/apps/bowtie2/2.3.4.1/bowtie2"
    samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools"
    methylQA = "/bar/cfan/software/methylQA_0.2.1/methylQA"
    bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools"
    bedGraphToBigWig = "/bar/cfan/software/kent/bin/bedGraphToBigWig"
    sambamba = "~/software/sambamba/sambamba"
    bamCoverage = "/bar/cfan/anaconda2/envs/jupyter/bin/bamCoverage"
    tabix = "/bar/cfan/anaconda2/envs/jupyter/bin/tabix"
    bgzip = "/bar/cfan/anaconda2/envs/jupyter/bin/bgzip"

# my scripts:
bam_filter_by_reads = "~/R_for_bash/bam_filter_by_reads.R"
bw_zero = "~/R_for_bash/bw_zero.R"

help="""
Alignment pipeline for both ATAC and ChIP. 
For ATAC data: bowtie2 alignment -> methylqa-based dedup and reads shifting 
For ChIP data: bowtie2 alignment -> sambamba-based dedup.
For bowtie2 alignment: -a/-k not specified. Please don't change this because in the 
step-wise processing pipeline there is no filtering for secondary alignments.
"""

print(help)

def bt2_input(wildcards):
    pe = FASTQ + "/" + wildcards.sample + "_R1" + SUFFIX
    if os.path.exists(pe):
        if "noCutadapt" in config.keys():
            return expand(FASTQ + "/{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
        else:
            return expand(TRIM_DIR + "/{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
    else:
        return FASTQ + "/{sample}" + SUFFIX

rule full:
    input:
        BW_DIR + "/{sample}_mqa_cov_RPKM_75.bw",
        BW_DIR + "/{sample}_mqa_cov_None_1.bw",
        BW_DIR + "/{sample}_mqa_cov_zero.bed.gz",
        BOWTIE2_DIR + "/{sample}.bam",
        BW_DIR + "/{sample}.open.bw"
    output:
        "snake_res/{sample}_full.txt"
    shell:
        "touch {output}"

rule chip:
    input:
        BW_DIR + "/{sample}_cov_RPKM_75.bw",
        BW_DIR + "/{sample}_cov_None_1.bw",
        BOWTIE2_DIR + "/{sample}.bam"
    output:
        "snake_res/{sample}_chip.txt"
    shell:
        "touch {output}"
rule simple:
    input:
        BW_DIR + "/{sample}_cov_RPKM_75.bw"
    output:
        "snake_res/{sample}_simple.txt"
    shell:
        "touch {output}"

rule cutadapt:
    input: 
        expand(FASTQ + "/{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
    output:
        expand(TRIM_DIR + "/{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
    log:
        TRIM_DIR + "/{sample}.cutadapt.log"
    threads:
        threads_trim
    run:
        shell("""{cutadapt}  -a {adapter_1} -A {adapter_2} \
        --quality-cutoff=15,10 --minimum-length=36 \
        -j {threads} \
        -o {output[0]} -p {output[1]} \
        {input[0]} {input[1]} > {log}""")
rule bowtie2:
    input:
        bt2_input
    output:
        BOWTIE2_DIR + "/{sample}.bam"
    log:
        BOWTIE2_DIR + "/{sample}.log"
    threads:
        threads_bowtie2
    run:
        if (len(input) > 1):
            fastqs = " -1 " + input[0] + " -2 " + input[1]
        else:
            fastqs = " -U " + input[0]
        
        shell("""
            {bowtie2} -x ~/genomes/{genome}/bowtie2/{genome} {fastqs} \
            --very-sensitive --xeq --dovetail -p {threads} --mm {bowtie2_unal} 2>{log} | \
            {samtools} sort - -O bam -m 2G -o {output} -@ {threads}
        """)
        shell("""
            {samtools} index {output}
        """)

rule AS0:
    input:
        BOWTIE2_DIR + "/{sample}.bam"
    output:
        BOWTIE2_DIR + "/{sample}_AS0.bam"
    threads:
        threads_bowtie2
    log:
        BOWTIE2_DIR + "/{sample}_AS0.log"
    run:
        shell("""
            {sambamba} view -f bam -F '[AS] == 0' -t {threads} {input} 1>{output} 2>{log}
        """)
        shell("""
            {samtools} index {output}
        """)

rule mapq:
    input:
        BOWTIE2_DIR + "/{sample}.bam"
    output:
        BOWTIE2_DIR + "/{sample}_mapq8.bam"
    threads:
        threads_bowtie2
    log:
        BOWTIE2_DIR + "/{sample}_mapq8.log"
    run:
        shell("""
            {sambamba} view -f bam -F 'mapping_quality >= 8 ' -t {threads} {input} 1>{output} 2>{log}
        """)
        shell("""
            {samtools} index {output}
        """)

rule sambamba_mkdup:
    input:
        BOWTIE2_DIR + "/{sample}.bam"
    output:
        BOWTIE2_DIR + "/{sample}_mkdup.bam"
    threads:
        threads_bowtie2
    log:
        BOWTIE2_DIR + "/{sample}_mkdup.log"
    run:
        shell("""
            {sambamba} markdup -t {threads} {sambamba_options} {input} {output} 1>{log} 2>&1
        """)
        shell("""
            {samtools} index {output}
        """)

rule sambamba_rmdup:
    input:
        BOWTIE2_DIR + "/{sample}_mkdup.bam"
    output:
        BOWTIE2_DIR + "/{sample}_rmdup.bam"
    threads:
        threads_bowtie2
    log:
        BOWTIE2_DIR + "/{sample}_rmdup.log"
    run:
        shell("""
            {sambamba} view -f bam -F 'not duplicate' -t {threads} {input} 1>{output} 2>{log}
        """)
        shell("""
            {samtools} index {output}
        """)

rule methylQA:
    input:
        BOWTIE2_DIR + "/{sample}.bam"
    output:
        METHYLQA_DIR + "/{sample}.open.tmp.bed"
    log:
        METHYLQA_DIR + "/{sample}.log"
    params:
        rootname = METHYLQA_DIR + "/{sample}"
    run:
        # {singularity_} exec -B /scratch:/scratch ~/software/simg/ATAC_mm10_target_181103.simg \
        shell("""
            {methylQA} atac -X 38 -o {params.rootname} -Q 1 \
            {chrom_sizes} {input} 1>{log} 2>&1
        """)
        shell("""
            mv {params.rootname}.open.bed {output}
        """)

rule methylQA_filter:
    input:
        METHYLQA_DIR + "/{sample}.open.tmp.bed"
    output:
         METHYLQA_DIR + "/{sample}.open.bed"
    run:
        shell("""
            awk -F "\t" 'BEGIN {{OFS = FS}} $2 < 3000000000 && $3 < 3000000000 {{print $0}}' {input} > {output} && \
            rm {input}
        """)

rule cutsite_bw:
    input: 
        METHYLQA_DIR + "/{sample}.open.bed"
    output: 
        BW_DIR + "/{sample}.open.bw"
    params:
        bdg = BW_DIR + "/{sample}.open.bedGraph"
    log:
        BW_DIR + "/{sample}.open.bwgen.log"
    run:
        shell("""
            awk '{{mid=int(($3+$2)/2); if($6=="+") {{print $1"\t"mid"\t"mid+1"\t"1}} else {{print $1"\t"mid-1"\t"mid"\t"1}}}}' \
            {input} | LC_COLLATE=C sort -k1,1 -k2,2n | uniq -c | awk -F " " '{{print $2"\t"$3"\t"$4"\t"$1}}' > {params.bdg}
        """)
        # note: LC_COLLATE=C is added by FANC: this leads to a case sensitive sort on chromosome names. Otherwise bedGraphToBigWig wouldn't work
        shell("""
            {bedGraphToBigWig} {params.bdg} {chrom_sizes} {output} 1> {log} 2>&1
        """)

rule mqa_reads:
    input:
        METHYLQA_DIR + "/{sample}.open.bed"
    output:
        METHYLQA_DIR + "/{sample}.open.reads.txt"
    shell:
        "cut -f4 {input} > {output}"

rule bam_filter_methylQA:
    input:
        bam = BOWTIE2_DIR + "/{sample}.bam",
        reads = METHYLQA_DIR + "/{sample}.open.reads.txt"
    output:
        BOWTIE2_DIR + "/{sample}_mqa.bam"
    log:
        BOWTIE2_DIR + "/{sample}_mqa.log"
    run:
        shell("""
            {Rscript} --vanilla {bam_filter_by_reads} -b {input.bam} -r {input.reads} 1>{log} 2>&1 
        """)

rule bam_coverage:
    input:
        BOWTIE2_DIR + "/{sample}.bam"
    output:
        BW_DIR + "/{sample}_cov_{norm}_{res}.bw"
    threads:
        threads_bowtie2
    log:
        BW_DIR + "/{sample}_cov_{norm}_{res}.log"
    shell:
        "{bamCoverage} -b {input} -o {output} -p {threads} --normalizeUsing {wildcards.norm} -bs {wildcards.res} 1>{log} 2>&1"

rule zero_cov_regions:
    input:
        BW_DIR + "/{sample}_mqa_cov_None_1.bw"
    output:
        BW_DIR + "/{sample}_mqa_cov_zero.bed.gz"
    params:
        bed = BW_DIR + "/{sample}_mqa_cov_zero.bed"
    log:
        BW_DIR + "/{sample}_mqa_cov_zero.log"
    shell:
        "{Rscript} --vanilla {bw_zero} -b {input} -o {params.bed} -t {tabix} -z {bgzip} 1>{log} 2>&1"
