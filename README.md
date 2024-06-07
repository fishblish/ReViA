# search-for-regulatory-variants
Preparing a tool for searching regulatory variants in the genome.

All key functions are located in directory *src* in the file *firstmodule.py*, and their execution takes place in the file *main.py*.

In directory *data* are located example files, which can be used while running scripts. This directory does not contain following files: VCF with variants *Symfonia_all.hard_filtered.selected.vcf.gz*, promoter activity signal *h3k27ac_promoters_normalized_cov_with_noise.csv* and gene expression *transcrtpts_rnaseq_quantile_normalized.csv*, which was used to obtain results presented in *output* directory.

Having GATK and Annovar programs is necessary to run scripts. These programs were used in the following versions: GATK version 4.4.0.0, Annovar version: Date: 2020-06-07.
To use motif search it is also obligatory to have installed following R packages: motifbreakR (2.12.3), BSgenome.Hsapiens.UCSC.hg38 (1.4.5), MotifDb (1.40.0).

Directory *output* contains example results. The example results were obtained by running code with following inputs and default parameters:
- input_vcf = Symfonia_all.hard_filtered.selected.vcf.gz,
- active_promoters = active_promoters_20pts_1kbup_200bpdown_not_merged.bed
- enhancer_regions = brain_enhancers_active.bed
- enhancer_activity = h3k27ac_coverage_quantile_normalized.csv
- promoter_activity = h3k27ac_promoters_normalized_cov_with_noise.csv
- gene_expression = transcrtpts_rnaseq_quantile_normalized.csv
- chromatin_loops = looplist_outside_promoters.txt.gz
- genes_info = hg38_genes_transcripts_info.tsv

Example program execution can look like that: 
```
/bin/python3 {base_path}/src/main.py --freq_filter_cutoff 0.01 --population ALL NFE --reference_population ALL --input_vcf {input_vcf} --enhancer_activity {enhacer_activity} --promoter_activity {promoter_activity} --gene_expression {gene_expression}
```
where names in brackets stand for paths to corresponding files.
Computing time for the VCF file 'Symfonia_all.hard_filtered.selected.vcf.gz' is approximately 27 minutes with motif search and 6 minutes without it.  This file consists of around 290,000 variants and 25 samples.

The successive steps in the analysis are as follows:
- Selection of variants that are located in regulatory regions and are biallelic.
- Assignment of variant frequency in the population and filtering based on this value. This is done using ANNOVAR.
- Performing a binomial test - selecting variants that are significantly more often (or rarely) present in the studied group than in the population. Benjamin-Hochberg correction is applied.
- Assignment of motifs which are influenced by variants.
- Assignment of transcripts to promoters (overlaping and those in chromatin contact).
- Assignment of transcripts to enhancers: overlaping, the closest ones and those in chromatin contact.
- Determination of correlation between enhancer activity (h3k27ac signal) and gene expression.
- Determination of correlation between promoter activity and gene expression.
- Determination of correlation between genotype and gene expression (at the transcript level).

It is possible to conduct a limited analysis using only the VCF file, without including expression and acetylation data.






