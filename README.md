# search-for-regulatory-variants
Preparing a tool for searching regulatory variants in the genome.

All key functions are located in directory *src* in the file *firstmodule.py*, and their execution takes place in the file *main.py*.

In directory *data* are located example files, which can be used while running scripts. This directory does not contain .vcf file with variants *Symfonia_all.hard_filtered.selected.vcf.gz*, which was used to obtain results presented in *output* directory.

Having GATK and Annovar programs is necessary to run scripts. These programs were used in following versions: GATK version 4.4.0.0, Annovar version: Date: 2020-06-07.


Directory *output* contains example results. They was obtained by running code with following inputs and default parameters:
- input_vcf = Symfonia_all.hard_filtered.selected.vcf.gz,
- promoter_regions = brain_promoters_active.bed
- enhancer_regions = brain_enhancers_active.bed
- enhancer_activity = h3k27ac_coverage_quantile_normalized.csv
- gene_expression = transcrtpts_rnaseq_quantile_normalized.csv
- chromatin_contacts = predicted_contacts.bed
- genes_info = hg38_full.genes.gtf


The successive steps in the analysis are as follows:
- Selection of variants that are located in regulatory regions and are biallelic.
- Assignment of variant frequency in the population and filtering based on this value. This is done using ANNOVAR.
- Performing a binomial test - selecting variants that are significantly more often (or rarely) present in the studied group than in the population. Benjamin-Hochberg correction is applied.
- Assignment of genes to promoters and intronic enhancers.
- Assignment of genes to enhancers: the closest ones and those in chromatin contact.
- Determination of correlation between enhancer activity (h3k27ac signal) and gene expression. Selection of the best gene.
- Determination of correlation between genotype and gene expression (at the transcript level, followed by selecting the best transcript).








