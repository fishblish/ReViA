import argparse
import os
import firstmodule as fm
import pandas as pd

base_path = '/home/julia/Desktop/uni/projekt_warianty/all/brain_reg_var-main/'
programs_path = '/home/julia/Desktop/uni/projekt_warianty/programs'
data_path = base_path + '/data'


parser = argparse.ArgumentParser()
parser.add_argument("--gatk_path", type=str, nargs='?', default= programs_path + \
                    "/gatk/gatk-4.4.0.0/gatk", help="Path to directory containing GATK program.")
parser.add_argument("--annovar_path", type=str, nargs='?', default= programs_path + \
                    "/annovar/", help="Path to directory containing ANNOVAR program.")
parser.add_argument("--input_vcf", type=str, nargs='?', default= data_path + "/Symfonia_all.hard_filtered.selected.vcf.gz", \
                    help="Path to directory containing input VCF file. It can be gzipped.")
parser.add_argument("--promoter_regions", type=str, nargs='?', default= data_path + "/brain_promoters_active.bed", \
                    help="Path to directory containing bed file with promoter regions. Last column should contain gene names, comma separated if promoters of several genes overlap.")
parser.add_argument("--enhancer_regions", type=str, nargs='?', default= data_path + "/brain_enhancers_active.bed", \
                    help="Path to directory containing bed file with enhancer regions. It should have 4 columns: chr, start, end, gene.\
                        Last column should contain gene names if enhancer is located inside gene, \
                        comma separated, '.' for intergenic enhancers = no gene overlaps.")
parser.add_argument("--output", type=str, nargs='?', default= "output", help="Name of directory to save outputs in.")
parser.add_argument("--enhancer_activity", type=str, nargs='?', default= data_path + "/h3k27ac_coverage_quantile_normalized.csv", \
                    help="Path to directory containing csv file with enhancer activity. It should have columns 'chr', 'start', 'end' and sample names.")
parser.add_argument("--gene_expression", type=str, nargs='?', default= data_path + "/transcrtpts_rnaseq_quantile_normalized.csv", \
                    help="Path to directory containing csv file with gene expression. It should have column 'Transcript' with data in form transcipt/gene name.")
parser.add_argument("--chromatin_contacts", type=str, nargs='?', default= data_path + "/predicted_contacts.bed", \
                    help="Path to directory containing bed file with chromatin contacts.")
parser.add_argument("--genes_info", type=str, nargs='?', default= data_path + '/hg38_full.genes.gtf', \
                    help="Path to directory containing file with information about genes.")
parser.add_argument("--default_enhancers_genes_assignment", type=str, default='output/assigned_genes_to_enhancers.csv', nargs='?', \
                    help="You can provide path to already assigned genes to enhancers. It should have columns: enh_start, enh_end, Gene, genomic element, H3K27ac-expression correlation p-values, relation.\
                        If you won't provide path, genes will be assigned to enhancers during analysis.")
# 

parser.add_argument("--freq_filter_target", type=str, choices=['r', 'c'], default='r', nargs='?', \
                    help="Choose 'r' (rare) or 'c' (common). It expresses if you want to select rare or common variants considering frequency in population")
parser.add_argument("--freq_filter_cutoff", type=float, default=0.01, choices=range(0,1), nargs='?', \
                    help="Set cut-off point for frequency filter. It expresses how rare or how common variant you what to select (rare/common depends on freq_filter_target argument)")
parser.add_argument("--population", type=str, nargs='*', choices=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL'], default=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL'],\
                    help="Choose population which input data will be copared with from a list ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL']. You can pass multiple arguments. If you won't pass any argument, whole population frequency will be used.")
parser.add_argument("--freq_filter_missing", type=str, choices=['r', 'c'], default='c', nargs='?', \
                    help="Choose 'r' (rare) or 'c' (common). It expresses if you want to treat variant with missing frequency data as rare or common variant.")
parser.add_argument("--reference_population", type=str, choices=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL'], default='ALL', nargs='?', \
                    help="Choose reference population for binomial test. It should be in population list: ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL']. If you won't pass any argument, frequency from whole population will be used.")
parser.add_argument("--bh_alpha", type=float, default=0.01, choices=range(0,1), nargs='?', \
                    help="Set cut-off point for Benjamini-Hochberg correction.")



args = parser.parse_args()

ANNOVAR = args.annovar_path
GATK = args.gatk_path
INPUT_VCF = args.input_vcf
PROMOTER_REGIONS = args.promoter_regions
ENHANCER_REGIONS = args.enhancer_regions
OUTPUT = args.output
ENHANCER_ACTIVITY = args.enhancer_activity
GENE_EXPRESSION = args.gene_expression
CHROMATIN_CONTACTS = args.chromatin_contacts
GENES_INFO = args.genes_info
DEFAULT_ENHANCERS_GENES_ASSIGNMENT = args.default_enhancers_genes_assignment
freq_filter_target = args.freq_filter_target
freq_filter_cutoff = args.freq_filter_cutoff
population = args.population
freq_filter_missing = args.freq_filter_missing
reference_population = args.reference_population
bh_alpha = args.bh_alpha

if reference_population != 'ALL':
    assert reference_population in population, "Reference population must be in population list"


try:
    os.mkdir(OUTPUT)
except:

    print('Directory ', OUTPUT, ' already exist. Output files will be saved in this localization.')

'''
#check if ANNOVAR and GATK are installed in provided localization
flags = fm.check_programs(GATK, ANNOVAR) 
if flags!= [0,0]:
    from itertools import compress
    raise TypeError('You should provide correct path to program(s)', list(compress(['GATK', 'ANNOVAR'], flags)))

if fm.prepare_annovar_database(ANNOVAR)==0:
    raise TypeError(f'ANNOVAR database cound not be installed. You can install it manually from annovar website. It should be placed in {ANNOVAR}/humandb directory and unpacked.')

#check if input .vcf file is indexed - it has to be to use some GATK functions
fm.index_input_vcf(INPUT_VCF, GATK)

# print how many variants are in input file
fm.count_variants(GATK, INPUT_VCF) 

#select biallelic variants located inside regulatory regions
biallelic_variants_in_regulatory_regions_path = fm.select_biallelic_inside_regulatory_regions(GATK, INPUT_VCF, PROMOTER_REGIONS, ENHANCER_REGIONS, OUTPUT)

#found variants that are enriched in analyzed cohort
annotated_variants_path = fm.annotate_freq(biallelic_variants_in_regulatory_regions_path, ANNOVAR)

#filter variants by frequency
filtered_by_freq_variants_file = fm.select_snps_by_freq(annotated_variants_path, GATK, target=freq_filter_target, cutoff=freq_filter_cutoff, population=population, treating_gaps=freq_filter_missing)

#select enriched snps
enriched_promoter_snps_df, enriched_enhancer_snps_df = fm.select_enriched_snps(filtered_by_freq_variants_file, GATK, bh_alpha=bh_alpha, target=freq_filter_target, reference_population=reference_population)

#assigning genes to promoters and enhancers
enriched_promoter_snps_gene = fm.assign_genes_to_promoter_snps(enriched_promoter_snps_df, PROMOTER_REGIONS)

enriched_enhancer_snps_gene = fm.assign_genes_intronic_enhancer_snps(enriched_enhancer_snps_df, ENHANCER_REGIONS)

DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available = 0
if DEFAULT_ENHANCERS_GENES_ASSIGNMENT:
    try:
        df = pd.read_csv(DEFAULT_ENHANCERS_GENES_ASSIGNMENT)
        assert all(element in df.columns for element in ['enh_start', 'enh_end', 'Gene', 'genomic element', 'H3K27ac-expression correlation p-values', 'relation'])
        DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available=1
        print('Default enhancer-gene assignment will be used.')
    except:
        pass

if DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available==0:
    genes_info_prepared = fm.prepare_genes_info(GENES_INFO)
    enriched_enhancer_snps_gene_closest = \
        fm.assign_closest_gene_to_enhancers(enriched_enhancer_snps_gene, genes_info_prepared)

    enriched_enhancer_snps_gene_closest_contacts = \
        fm.assign_chromatin_contacting_gene_to_enhancer(enriched_enhancer_snps_gene_closest, genes_info_prepared, CHROMATIN_CONTACTS)
    enriched_enhancer_snps_genes_collected = \
        fm.reformat_target_genes_enh(enriched_enhancer_snps_gene_closest_contacts)



    #calculate correlation between enhancer activity and gene expression
    enriched_enhancer_snps_genes_collected_corelations = \
        fm.check_signal_gene_expression_correlation_enhancer(enriched_enhancer_snps_genes_collected, ENHANCER_ACTIVITY, GENE_EXPRESSION)

    #change format to store one gene per row
    enhancer_snps = fm.change_table_format_enh(enriched_enhancer_snps_genes_collected_corelations)

    assigned_genes_to_enhancers = enhancer_snps[['enh_start','enh_end','Gene','genomic element','H3K27ac-expression correlation p-values','relation']].drop_duplicates()
    assigned_genes_to_enhancers.to_csv(OUTPUT + '/assigned_genes_to_enhancers.csv')

else:
    assigned_genes_to_enhaners = pd.read_csv(DEFAULT_ENHANCERS_GENES_ASSIGNMENT)
'''
#saved intermediate results to testing
# enriched_enhancer_snps_gene.to_csv(OUTPUT + '/middle_result_enriched_enhancer_snps.csv')
# enriched_promoter_snps_gene.to_csv(OUTPUT + '/middle_result_enriched_promoter_snps.csv')

assigned_genes_to_enhaners = pd.read_csv(DEFAULT_ENHANCERS_GENES_ASSIGNMENT)
enriched_enhancer_snps_gene = pd.read_csv(OUTPUT + '/middle_result_enriched_enhancer_snps.csv')
enriched_promoter_snps_gene = pd.read_csv(OUTPUT + '/middle_result_enriched_promoter_snps.csv')
    


#enriched_enhancer_snps_df join with assigned_genes_to_enhancers
enriched_enhancer_snps = enriched_enhancer_snps_gene.drop(columns=['Gene'])
enriched_enhancer_snps_genes_collected_corelations = enriched_enhancer_snps.merge(assigned_genes_to_enhaners, left_on=['enh_start', 'enh_end', 'genomic element'], right_on=['enh_start', 'enh_end','genomic element'], how='right')

enriched_promoter_snps_gene = fm.change_table_format_promoter(enriched_promoter_snps_gene)
#calculate correlation between genotype and gene expression
promoter_snps_sample_lvl, enhancer_snps_sample_lvl = fm.import_vcf_sample_level(INPUT_VCF, OUTPUT, GATK, enriched_promoter_snps_gene, enriched_enhancer_snps_genes_collected_corelations)
promoter_snps, enhancer_snps = fm.check_gene_genotype_correlation(GENE_EXPRESSION, promoter_snps_sample_lvl, enhancer_snps_sample_lvl)

enhancer_snps = fm.check_genotype_enhancer_correlation(enhancer_snps,ENHANCER_ACTIVITY)

#remove unnecessary files created during analysis
fm.remove_unnecessary_files(OUTPUT)

#save and visualize results
fm.visualize_results(promoter_snps, enhancer_snps, GENE_EXPRESSION, OUTPUT,ENHANCER_ACTIVITY)

