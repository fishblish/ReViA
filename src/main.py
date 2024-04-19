import argparse
import os
import firstmodule as fm
import pandas as pd
import sys
from itertools import compress
import time
start_time = time.ctime()

base_path = os.path.dirname(os.path.abspath(__file__)).replace('/src', '')
programs_path = '/home/julia/Desktop/uni/projekt_warianty/programs'
data_path = base_path + '/data'
print('Default input data will be taken from: ', data_path)


parser = argparse.ArgumentParser()
parser.add_argument("--gatk_path", type=str, nargs='?', default= programs_path + \
                    "/gatk/gatk-4.4.0.0/gatk", help="Path to directory containing GATK program.")
parser.add_argument("--annovar_path", type=str, nargs='?', default= programs_path + \
                    "/annovar/", help="Path to directory containing ANNOVAR program.")
parser.add_argument("--input_vcf", type=str, nargs='?', default= data_path + "/Symfonia_all.hard_filtered.selected.vcf.gz", \
                    help="Path to directory containing input VCF file. It can be gzipped.")
parser.add_argument("--promoter_regions", type=str, nargs='?', default= data_path + "/promoter_transcripts_hg38_full.bed", \
                    help="Path to directory containing bed file with promoter regions. Last column should contain gene names, comma separated if promoters of several genes overlap.")
parser.add_argument("--enhancer_regions", type=str, nargs='?', default= data_path + "/brain_enhancers_active.bed", \
                    help="Path to directory containing bed file with enhancer regions. It should contain this 3 columns: chr, start, end.\
                        Last column should contain gene names if enhancer is located inside gene, \
                        comma separated, '.' for intergenic enhancers = no gene overlaps.")
parser.add_argument("--transcripts_regions", type=str, nargs='?', default= data_path + "/transcripts.bed", \
                    help="Path to directory containing bed file with transcripts regions. It should have columns 'chr', 'start', 'end' and 'Transcript'.")
parser.add_argument("--output", type=str, nargs='?', default= "output", help="Name of directory to save outputs in.")
parser.add_argument("--enhancer_activity", type=str, nargs='?', default= data_path + "/h3k27ac_coverage_quantile_normalized.csv", \
                    help="Path to directory containing csv file with signal h3k27ac for enhancers. It should have columns 'chr', 'start', 'end' and sample names.")
parser.add_argument("--promoter_activity", type=str, nargs='?', default= data_path + "/h3k27ac_promoters_normalized_cov_with_noise.csv", \
                    help="Path to directory containing csv file with signal h3k27ac for promoters. It should have column 'Transcript' and sample names.")

parser.add_argument("--gene_expression", type=str, nargs='?', default= data_path + "/transcrtpts_rnaseq_quantile_normalized.csv", \
                    help="Path to directory containing csv file with gene expression. It should have column 'Transcript' with data in form transcipt/gene name.")
parser.add_argument("--chromatin_loops", type=str, nargs='?', default= data_path + '/loop_lists/GSE63525_NHEK_HiCCUPS_looplist.txt.gz',help="Path to directory containing file with chromatin loops. It should have tab as separator.")
parser.add_argument("--genes_info", type=str, nargs='?', default= data_path + '/hg38_full.transcripts.protein_coding.gtf', \
                    help="Path to directory containing file with information about genes.")
parser.add_argument("--default_enhancers_genes_assignment",  type=str, nargs='?', \
                    help="You can provide path to already assigned genes to enhancers. It should have columns: enh_start, enh_end, Gene, genomic element, H3K27ac-expression correlation p-values, relation.\
                        If you won't provide path, genes will be assigned to enhancers during analysis.")
# default='output/assigned_genes_to_enhancers.csv',

parser.add_argument("--freq_filter_target", type=str, choices=['r', 'c'], default='r', nargs='?', \
                    help="Choose 'r' (rare) or 'c' (common). It expresses if you want to select rare or common variants considering frequency in population")
parser.add_argument("--freq_filter_cutoff", type=float, default=0.01, nargs='?', \
                    help="Set cut-off point for frequency filter. It expresses how rare or how common variant you what to select (rare/common depends on freq_filter_target argument)")
parser.add_argument("--population", type=str, nargs='*', choices=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL'], default=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL'],\
                    help="Choose population which input data will be copared with from a list ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL']. You can pass multiple arguments. If you won't pass any argument, whole population frequency will be used.")
parser.add_argument("--freq_filter_missing", type=str, choices=['r', 'c'], default='c', nargs='?', \
                    help="Choose 'r' (rare) or 'c' (common). It expresses if you want to treat variant with missing frequency data as rare or common variant.")
parser.add_argument("--reference_population", type=str, choices=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL'], default='ALL', nargs='?', \
                    help="Choose reference population for binomial test. It should be in population list: ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'ALL']. If you won't pass any argument, frequency from whole population will be used.")
parser.add_argument("--bh_alpha", type=float, default=0.01, nargs='?', \
                    help="Set cut-off point for Benjamini-Hochberg correction.")
parser.add_argument("-limited_analysis", action='store_true', help="If you have only input vcf file, you can use this flag to run analysis without enhancer activity, gene expression.")

parser.add_argument("--motifs", action='store_false', help="If you want to search for motifs in enhancers and promoters.")
parser.add_argument("--r_packages_path", type=str, nargs='?', default="/home/julia/miniconda3/lib/R/library/", help="Path to directory containing R packages: motifbreakR, BSgenome.Hsapiens.UCSC.hg38, MotifDb. Providing this path is necessary to search for motifs.")

args = parser.parse_args()

ANNOVAR = args.annovar_path
GATK = args.gatk_path
INPUT_VCF = args.input_vcf
PROMOTER_REGIONS = args.promoter_regions
ENHANCER_REGIONS = args.enhancer_regions
TRANSCRIPTS_REGIONS = args.transcripts_regions
OUTPUT = args.output
ENHANCER_ACTIVITY = args.enhancer_activity
PROMOTER_ACTIVITY = args.promoter_activity
GENE_EXPRESSION = args.gene_expression
CHROMATIN_LOOPS = args.chromatin_loops
GENES_INFO = args.genes_info
DEFAULT_ENHANCERS_GENES_ASSIGNMENT = args.default_enhancers_genes_assignment
freq_filter_target = args.freq_filter_target
freq_filter_cutoff = args.freq_filter_cutoff
population = args.population
freq_filter_missing = args.freq_filter_missing
reference_population = args.reference_population
bh_alpha = args.bh_alpha
limited_analysis = args.limited_analysis
motifs = args.motifs
R_PACKAGES_PATH = args.r_packages_path

if reference_population != 'ALL':
    assert reference_population in population, "Reference population must be in population list"

assert freq_filter_cutoff >= 0 and freq_filter_cutoff <= 1, "Frequency filter cut-off must be between 0 and 1."
assert bh_alpha >= 0 and bh_alpha <= 1, "Benjamini-Hochberg alpha must be between 0 and 1."

try:
    os.mkdir(OUTPUT)
except:

    print('Directory ', OUTPUT, ' already exist. Output files will be saved in this localization.')

#check if input files are in provided localizations
fm.check_input_files([INPUT_VCF, PROMOTER_REGIONS, ENHANCER_REGIONS, ENHANCER_ACTIVITY, PROMOTER_ACTIVITY, GENE_EXPRESSION, CHROMATIN_LOOPS, GENES_INFO], DEFAULT_ENHANCERS_GENES_ASSIGNMENT)    
genes_info_prepared = fm.prepare_genes_info(GENES_INFO)

DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available = 0
if DEFAULT_ENHANCERS_GENES_ASSIGNMENT:
    try:
        df = pd.read_csv(DEFAULT_ENHANCERS_GENES_ASSIGNMENT)
        assert all(element in df.columns for element in ['enh_start', 'enh_end', 'Gene','Transcript', 'genomic element', 'H3K27ac-expression correlation p-values', 'relation'])
        DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available=1
        print('Default enhancer-gene assignment will be used.')
    except:
        pass

#check if ANNOVAR and GATK are installed in provided localization
flags = fm.check_programs(GATK, ANNOVAR) 
if flags!= [0,0]:
    raise TypeError('You should provide correct path to program(s)', list(compress(['GATK', 'ANNOVAR'], flags)))

#check if R packages are installed
if motifs:
    flags = fm.check_R_packages(R_PACKAGES_PATH)
    if flags!= [0,0,0]:
        raise TypeError('You should provide correct path to R packages:', list(compress(['motifbreakR', 'BSgenome.Hsapiens.UCSC.hg38', 'MotifDb'], flags)))

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
'''
enriched_enhancer_snps_df.to_csv(OUTPUT + '/enriched_enhancer_snps.csv')
enriched_promoter_snps_df.to_csv(OUTPUT + '/enriched_promoter_snps.csv')

enriched_enhancer_snps_df = pd.read_csv(OUTPUT + '/enriched_enhancer_snps.csv')
enriched_promoter_snps_df = pd.read_csv(OUTPUT + '/enriched_promoter_snps.csv')
'''
if motifs:
    #select SNPs in motifs
    enriched_promoter_snps_df, enriched_enhancer_snps_df = fm.find_motifs(enriched_promoter_snps_df, enriched_enhancer_snps_df, OUTPUT, R_PACKAGES_PATH)

#assigning genes to promoters and enhancers
enriched_promoter_snps_gene = fm.assign_genes_to_promoter_snps(enriched_promoter_snps_df, PROMOTER_REGIONS)
enriched_promoter_snps_gene = fm.assign_chromatin_contacting_gene_to_promoter_with_loops(enriched_promoter_snps_gene, TRANSCRIPTS_REGIONS, CHROMATIN_LOOPS)


enriched_enhancer_snps_gene = fm.assign_genes_intronic_enhancer_snps(enriched_enhancer_snps_df, ENHANCER_REGIONS, genes_info_prepared)
enriched_promoter_snps_gene = fm.change_table_format_promoter(enriched_promoter_snps_gene, genes_info_prepared)
if limited_analysis:
    enriched_enhancer_snps_gene_closest = \
        fm.assign_closest_gene_to_enhancers(enriched_enhancer_snps_gene, genes_info_prepared)

    enriched_enhancer_snps_gene_closest_contacts = \
        fm.assign_chromatin_contacting_gene_to_enhancer_with_loops(enriched_enhancer_snps_gene_closest, TRANSCRIPTS_REGIONS, CHROMATIN_LOOPS)
    enriched_enhancer_snps_genes_collected = \
        fm.reformat_target_genes_enh(enriched_enhancer_snps_gene_closest_contacts, genes_info_prepared)

    fm.remove_unnecessary_files(OUTPUT)
    if DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available==1:
        assigned_genes_to_enhancers = pd.read_csv(DEFAULT_ENHANCERS_GENES_ASSIGNMENT)
        enriched_enhancer_snps_genes_collected = enriched_enhancer_snps_genes_collected.merge(assigned_genes_to_enhancers[['enh_start','enh_end','H3K27ac-expression correlation p-values']], on=['enh_start', 'enh_end'], how='left').fillna('.')
    fm.save_limited_results(enriched_promoter_snps_gene, enriched_enhancer_snps_genes_collected, OUTPUT)
    sys.exit("Limited analysis finished.")



if DEFAULT_ENHANCERS_GENES_ASSIGNMENT_is_available==0: # czy skorzystać w tym kroku jeszcze z korelacji między ekspresją a aktywnością enhancerów wyznaczoną na innych danych?
    
    enriched_enhancer_snps_gene_closest = \
        fm.assign_closest_gene_to_enhancers(enriched_enhancer_snps_gene, genes_info_prepared)

    enriched_enhancer_snps_gene_closest_contacts = \
        fm.assign_chromatin_contacting_gene_to_enhancer_with_loops(enriched_enhancer_snps_gene_closest, TRANSCRIPTS_REGIONS, CHROMATIN_LOOPS)
    enriched_enhancer_snps_genes_collected = \
        fm.reformat_target_genes_enh(enriched_enhancer_snps_gene_closest_contacts, genes_info_prepared)


    #calculate correlation between enhancer activity and gene expression
    enriched_enhancer_snps_genes_collected_corelations = \
        fm.check_signal_gene_expression_correlation_enhancer(enriched_enhancer_snps_genes_collected, ENHANCER_ACTIVITY, GENE_EXPRESSION)

    assigned_genes_to_enhancers = enriched_enhancer_snps_genes_collected_corelations[['enh_start','enh_end','Transcript','Gene','Gene_ID','genomic element','H3K27ac-expression correlation p-values','relation']].drop_duplicates()
    assigned_genes_to_enhancers.to_csv(OUTPUT + '/assigned_genes_to_enhancers.csv')
    DEFAULT_ENHANCERS_GENES_ASSIGNMENT = OUTPUT + '/assigned_genes_to_enhancers.csv'

else:
    assigned_genes_to_enhaners = pd.read_csv(DEFAULT_ENHANCERS_GENES_ASSIGNMENT)


#enriched_enhancer_snps_df join with assigned_genes_to_enhancers
enriched_enhancer_snps = enriched_enhancer_snps_gene.drop(columns=['Transcript'])
enriched_enhancer_snps_genes_collected_corelations = enriched_enhancer_snps.merge(assigned_genes_to_enhancers, left_on=['enh_start', 'enh_end', 'genomic element'], right_on=['enh_start', 'enh_end','genomic element'], how='right')


#calculate correlation between genotype and gene expression
promoter_snps_sample_lvl, enhancer_snps_sample_lvl = fm.import_vcf_sample_level(INPUT_VCF, OUTPUT, GATK, enriched_promoter_snps_gene, enriched_enhancer_snps_genes_collected_corelations)
promoter_snps, enhancer_snps = fm.check_gene_genotype_correlation(GENE_EXPRESSION, promoter_snps_sample_lvl, enhancer_snps_sample_lvl)

enhancer_snps, promoter_snps = fm.check_genotype_signal_correlation(enhancer_snps,promoter_snps,ENHANCER_ACTIVITY, PROMOTER_ACTIVITY)

#remove unnecessary files created during analysis
fm.remove_unnecessary_files(OUTPUT)

#save and visualize results
fm.visualize_results(promoter_snps, enhancer_snps, GENE_EXPRESSION, OUTPUT,ENHANCER_ACTIVITY, PROMOTER_ACTIVITY)
print('Time of analysis: start - ',start_time,' end - ', time.ctime(), )

