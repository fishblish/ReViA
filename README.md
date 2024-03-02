# search-for-regulatory-variants
Preparing a tool for searching regulatory variants in the genome.

All key functions are located in the file firstmodule.py, and their execution takes place in the file main.py.

The successive steps in the analysis are as follows:
- Selection of variants that are located in regulatory regions and are biallelic.
- Assignment of variant frequency in the population and filtering based on this value. This is done using ANNOVAR.
- Performing a binomial test - selecting variants that are significantly more often (or rarely) present in the studied group than in the population. Benjamin-Hochberg correction is applied.
- Assignment of genes to promoters and intronic enhancers.
- Assignment of genes to enhancers: the closest ones and those in chromatin contact.
- Determination of correlation between enhancer activity (h3k27ac signal) and gene expression. Selection of the best gene.
- Determination of correlation between genotype and gene expression (at the transcript level, followed by selecting the best transcript).








