# Manual
A snakemake pipeline for ATAC-seq or ChIP/Cut&Run. It's designed to run on either Wang Lab server or WashU RIS. The most up-to-date version is at `/bar/cfan/software/snakeATAC/`. Each smk file is a version of the pipeline. Currently the most up-to-date pipeline is in `snakeATAC_simple_wildcard.smk`.
## Dependencies
* snakemake (tested on version 6.12.3)
* To run it without being Changxu Fan, you perhaps need to modify the paths to these other dependencies:
```
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
```
* someother R scripts that you need to copy from `/bar/cfan/R_for_bash/`:
```
bam_filter_by_reads = "~/R_for_bash/bam_filter_by_reads.R"
bw_zero = "~/R_for_bash/bw_zero.R"
```
* in addition, your bowtie2 index needs to be located in: `~/genomes/${genome}/bowtie2/${genome}`. Take a look at `/bar/cfan/genomes/mm10/bowtie2/mm10*` for reference.

## Running the pipeline:
Note the logic of snakemake: the pipeline is dictated by what you specify as output. Therefore, the snakeATAC pipeline can be run in many different ways, depending on how you specify the output. The output themselves are just empty txt files. But the names of these files tell the pipeline what to do. 

### input:
* `fastq/${sample}_R1.fastq.gz`; `fastq/${sample}_R2.fastq.gz`
### possible outputs:
* full output (ATAC only): `snake_res/${sample}_AS0_mapq8_full.txt`
    * AS0: only retain reads that matches the references perfectly.
    * mapq8: only retain reads with mapq > 8.
    * Note that `full` indicates that `methylQA` will be run, which removes PCR duplicates.
* partial output can be specified to skip some steps. For example, when aligning human data you probably want to specify `snake_res/${sample}_mapq8_full.txt` and skip the AS0 filtering since humans have SNPs.
* ChIP/Cut&Run vs ATAC: 
    * the pipeline is by default for ATAC data. After alignment, methylQA (by Daofeng Li) is used for Tn5 cutsite identification and PCR duplicate removal.
    * for ChIP/Cut&Run, specify `chip` such as `snake_res/${sample}_mapq8_rmdup_chip.txt`. This will skip the methylQA step. Deduplication is instead performed by `sambamba`.

### running the pipeline:
Most common example:
```
snakemake snake_res/${sample}_mapq8_rmdup_chip.txt -p --snakefile snakeATAC_simple_wildcard.smk -c 12 --config genome=hg38
```
* before actually running the pipeline, add `-n` for a dry-run first!

