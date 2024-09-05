import pandas as pd
import numpy as np
from scipy.stats import binomtest, spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools as pbt
import os
import subprocess
from pybedtools import BedTool
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def check_programs(GATK, ANNOVAR):
    cmd = GATK + ' --help'
    flags = [1, 1]
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
        print('GATK is installed correctly.')
        flags[0] = 0
    except subprocess.CalledProcessError as e:
        print('GATK is not installed correctly in provided localization.')
        print(f"Command '{''.join(cmd)}' failed with error code {e.returncode}:\n{e.output}")
    cmd = 'perl ' + ANNOVAR + 'annotate_variation.pl -h'
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
        print('ANNOVAR is installed correctly.')
        flags[1] = 0
    except subprocess.CalledProcessError as e:
        if e.returncode == 1:
            print('ANNOVAR is installed correctly.')
            flags[1] = 0

        else:
            print('ANNOVAR is not installed correctly in indicated localization.')
            print(f"Command '{''.join(cmd)}' failed with error code {e.returncode}:\n{e.output}")

    return flags

def check_input_files(files_to_check_list , DEFAULT_ENHANCERS_GENES_ASSIGNMENT):
    for path in files_to_check_list:
        assert os.path.isfile(path), f'File {path} does not exist.'
        
    if DEFAULT_ENHANCERS_GENES_ASSIGNMENT:
        assert os.path.isfile(DEFAULT_ENHANCERS_GENES_ASSIGNMENT), f'File {DEFAULT_ENHANCERS_GENES_ASSIGNMENT} does not exist.'


def prepare_genes_info(GENES_INFO):
    genes_info = pd.read_csv(GENES_INFO)
    assert all(column in genes_info.columns for column in ['chr', 'start','end','strand','Gene_name','Gene_ID','Transcript']) , 'Columns in genes_info are not correct'
    genes_info = genes_info[['chr', 'start','end','strand','Gene_name','Gene_ID','Transcript']].copy()
    genes_info['tss'] = genes_info.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)
    return genes_info

def select_promoter_regions(genes_info, promoter_range, active_promoters, active_promoter_data_type):
    promoter_regions = genes_info[['chr', 'tss','Transcript']].copy()
    promoter_regions['ENST'] = promoter_regions['Transcript'].apply(lambda x: x.split('/')[0])
    promoter_regions['start'] = (promoter_regions['tss'] - promoter_range).clip(lower=0)
    promoter_regions['end'] = promoter_regions['tss'] + promoter_range
    
    if active_promoters:
        active_promoters = pd.read_csv(active_promoters, sep='\t', header=None)
        if active_promoter_data_type == 'ENST':
            col_number = [i for i in range(active_promoters.shape[1]) if 'ENST' in str(active_promoters.iloc[0,i])]
            assert len(col_number)!=0, 'No column with ENST IDs in provided file.'
            promoter_regions = promoter_regions[promoter_regions['ENST'].isin(active_promoters[col_number[0]])]
        if active_promoter_data_type == 'TSS':
            promoter_regions = promoter_regions[promoter_regions['tss'].isin(active_promoters[1]) and promoter_regions['chr'].isin(active_promoters[0])]
    promoter_regions = promoter_regions.drop(columns=['tss','ENST'])
    return promoter_regions

def index_input_vcf(INPUT_VCF, GATK):
    vcf_path = '/'.join(INPUT_VCF.split('/')[0:-1])
    vcf_file_name = INPUT_VCF.split('/')[-1]

    if (vcf_file_name+'.tbi' not in os.listdir(vcf_path)):
        print('Indexing vcf file')
        cmd = GATK + ' IndexFeatureFile -I ' + INPUT_VCF
        try:
            result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,
                                                universal_newlines=True).splitlines()
            print('VCF file indexed')
        except subprocess.CalledProcessError as e:
            print(f"Command '{cmd}' failed with error code {e.returncode}:\n{e.output}")


#counting and selectin variants
def count_variants(GATK, INPUT_VCF):

    cmd = GATK + ' CountVariants -V ' + INPUT_VCF

    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,
                                         universal_newlines=True).splitlines()
        print("Number of variants:", result[-1])
    except subprocess.CalledProcessError as e:
        print(f"Command '{cmd}' failed with error code {e.returncode}:\n{e.output}")
        result = None


def select_biallelic_inside_regulatory_regions(GATK, INPUT_VCF, promoter_regions,PROMOTER_REGIONS_BED_PATH, ENHANCER_REGIONS, OUTPUT):
    select_logs = []
    count_logs = []
    result_files={}

    # create bed file with promoter_regions
    promoter_regions[['chr','start','end','Transcript']].to_csv(PROMOTER_REGIONS_BED_PATH, sep='\t', index=False, header=False)

    for r, regions in [("promoter", PROMOTER_REGIONS_BED_PATH), ("enhancer", ENHANCER_REGIONS)]:
        result_files[r] = f'{OUTPUT}/{r}_intermediate_result.vcf'
        command1 = f"{GATK} SelectVariants -V {INPUT_VCF} -L {regions} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {result_files[r]}"
        print('Selecting biallelic variants for', r)
        try:
            log1 = subprocess.check_output(command1, shell=True, stderr=subprocess.STDOUT,
                                           universal_newlines=True).splitlines()
            select_logs.append(str(log1))
        except subprocess.CalledProcessError as e:
            print(f"Command '{command1}' failed with error code {e.returncode}:\n{e.output}")

        command2 = f"{GATK} CountVariants -V {result_files[r]}"
        log2 = subprocess.check_output(command2, shell=True, stderr=subprocess.STDOUT,
                                       universal_newlines=True).splitlines()
        count_logs.append(log2)
        print(f"Number of biallelic SNPs in {r} regions:{log2[-1]} \n")
    logs = [select_logs, count_logs]
    return result_files


#
# #ANNOVAR - annotate frequencies
# #download database
def prepare_hg38_annovar_database(ANNOVAR):
    downloaded_flag=1
    if ('humandb' in os.listdir(ANNOVAR)):
        if ('hg38_gnomad_genome.txt' in os.listdir(ANNOVAR + '/humandb')):
            print('You already have human genome 38 database ready to use.')
            return 0
        if ('hg38_gnomad_genome.txt.gz' in os.listdir(ANNOVAR + '/humandb')):
            print('You already have human genome 38 database. It need to be unpacked. Do you want to unpack it now? (y/n)')
            answer = input()
            if answer == 'y':
                command = f'gunzip {ANNOVAR}/humandb/hg38_gnomad_genome.txt.gz'
                try:
                    logs = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                                    universal_newlines=True).splitlines()
                    print('Annovar database human genome 38 unpacked under path:', ANNOVAR + '/humandb/')
                    return 0
                except subprocess.CalledProcessError as e:
                    print('ANNOVAR database could not be unpacked. Please check if you have enough space on your disk.')
                    return 1
            else:
                return 1
            
    command = f'perl {ANNOVAR}annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome {ANNOVAR}/humandb/'
    try:
        print('Database need to be downloaded')
        logs = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
        print('Annovar database human genome 38 downloaded under path:', ANNOVAR + '/humandb/') #'hg38_gnomad_genome.txt' is 17G heavy
        downloaded_flag = 0
    except subprocess.CalledProcessError as e:
        print(f"Command '{command}' failed with error code {e.returncode}:\n{e.output}")
    return downloaded_flag # 0 - database downloaded, 1 - database not downloaded

def prepare_mm39_annovar_database(ANNOVAR):
    downloaded_flag=1
    if ('mousedb' in os.listdir(ANNOVAR)):
        if ('mm39_refGene.txt' in os.listdir(ANNOVAR + '/mousedb')):
            print('You already have mouse genome 39 database ready to use.')
            return 0

    command1 = f'perl {ANNOVAR}/annotate_variation.pl -downdb -buildver mm39 gene {ANNOVAR}/mousedb'
    command2 = f'perl {ANNOVAR}/annotate_variation.pl --buildver mm39 --downdb seq {ANNOVAR}/mousedb/mm39_seq'
    command3 = f'perl {ANNOVAR}/retrieve_seq_from_fasta.pl {ANNOVAR}/mousedb/mm39_refGene.txt -seqdir {ANNOVAR}/mousedb/mm39_seq -format refGene -outfile {ANNOVAR}/mousedb/mm39_refGeneMrna.fa'

    print('Database mm39 need to be downloaded and built.')
    try:
        logs = subprocess.check_output(command1, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Command '{command1}' failed with error code {e.returncode}:\n{e.output}")
    try:
        logs = subprocess.check_output(command2, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Command '{command2}' failed with error code {e.returncode}:\n{e.output}")
    try:
        logs = subprocess.check_output(command3, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
        print(f'Annovar database mouse genome 39 downloaded under path: {ANNOVAR}/mousedb/')
        downloaded_flag = 0
    except subprocess.CalledProcessError as e:
        print(f"Command '{command3}' failed with error code {e.returncode}:\n{e.output}")
    return downloaded_flag # 0 - database downloaded, 1 - database not downloaded

def prepare_dm3_annovar_database(ANNOVAR):
    downloaded_flag=1
    if ('flydb' in os.listdir(ANNOVAR)):
        if ('dm3_refGene.txt' in os.listdir(ANNOVAR + '/flydb')):
            print('You already have drosophila genome dm3 database ready to use.')
            return 0

    command1 = f'perl {ANNOVAR}/annotate_variation.pl -downdb -buildver dm3 gene {ANNOVAR}/flydb'
    command2 = f'perl {ANNOVAR}/annotate_variation.pl --buildver dm3 --downdb seq {ANNOVAR}/flydb/dm3_seq'
    command3 = f'perl {ANNOVAR}/retrieve_seq_from_fasta.pl {ANNOVAR}/flydb/dm3_refGene.txt -seqdir {ANNOVAR}/flydb/dm3_seq -format refGene -outfile {ANNOVAR}/flydb/dm3_refGeneMrna.fa'

    print('Database dm3 need to be downloaded and built.')
    try:
        logs = subprocess.check_output(command1, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Command '{command1}' failed with error code {e.returncode}:\n{e.output}")
    try:
        logs = subprocess.check_output(command2, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Command '{command2}' failed with error code {e.returncode}:\n{e.output}")
    try:
        logs = subprocess.check_output(command3, shell=True, stderr=subprocess.STDOUT,
                                            universal_newlines=True).splitlines()
        print(f'Annovar database drosophila genome 3 downloaded under path: {ANNOVAR}/flydb/')
        downloaded_flag = 0
    except subprocess.CalledProcessError as e:
        print(f"Command '{command3}' failed with error code {e.returncode}:\n{e.output}")
    return downloaded_flag # 0 - database downloaded, 1 - database not downloaded

def annotate_freq(variants_vcf_file_dict, ANNOVAR, genome_version):
    annovar_logs = []
    result_vcf_files_dict={}
    for r in ["promoter", "enhancer"]:
        out_file = variants_vcf_file_dict[r].split('.vcf')[0]
        result_vcf_files_dict[r] = f'{out_file}.{genome_version}_multianno.vcf'
        
        db_path = {'hg38': f'{ANNOVAR}/humandb/', 'mm39': f'{ANNOVAR}/mousedb/'}[genome_version]
        protocol = {'hg38': 'gnomad_genome', 'mm39': 'refGene'}[genome_version]

        command = f"perl {ANNOVAR}table_annovar.pl {variants_vcf_file_dict[r]} {db_path} -buildver {genome_version} -remove -protocol {protocol} -operation f -nastring . -vcfinput -out {out_file} -polish" 
        try:
            log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
            annovar_logs.append(log)
        except subprocess.CalledProcessError as e:
            print(f"Command '{command}' failed with error code {e.returncode}:\n{e.output}")
        print(f"Done: frequencies annotated for variants in {r} regions")
    return result_vcf_files_dict

# filter snps based on frequency in popultaion
# treating gaps - how to treat missing variants in annovar database
# target - rare or common snps
# cutoff - frequency cutoff
# for example: target='r', cutoff=0.01 will select snps with frequency below 0.01
# for example: target='c', cutoff=0.99 will select snps with frequency above 0.99
# population - list of populations to be used in frequency filter
def select_snps_by_freq(annotated_variants_file_dict, GATK, population, treating_gaps, target='r', cutoff=0.01):
    result_files_dict = annotated_variants_file_dict.copy()
    gaps_dict = {'r': f'={cutoff/2}', 'c': '=100.0'}
    ineq_sign = {'r': '<', 'c': '>'}
    target_full_word = {'r': 'rare', 'c': 'common'}
    population = list(set(population+['ALL']))

    annotations = ['gnomAD_genome_' + p for p in population]
    

    print('There will be selected variants which are', target_full_word[target], 'among these populations', population,
          'for frequency cutoff equal to', cutoff, 'and missing variants in annovar database treated as', target_full_word[treating_gaps])

    for r in ["promoter", "enhancer"]:
        with open(annotated_variants_file_dict[r], 'r') as i:
            lines = i.readlines()

        for i in range(len(lines)):
            for el in annotations:
                if el + '=.' in lines[i]:
                    lines[i] = lines[i].replace(el + '=.', el + gaps_dict[treating_gaps])

        with open(result_files_dict[r], 'w') as o:
            o.writelines(lines)

    select_rare_logs = []

    select_condition = '"'

    for a in annotations:
        select_condition = select_condition + a + ' ' + ineq_sign[target] + ' ' + str(cutoff) + ' && '
    select_condition = select_condition[0:-3]
    select_condition = select_condition + '"'
    result_filtered_files_dict = {}
    for r in ["promoter", "enhancer"]:
        result_filtered_files_dict[r] = result_files_dict[r].replace('.vcf','.filtered.vcf')

        # Select SNPs significantly enriched in analyzed cohort at specified frequency cutoff
        command = f'{GATK} SelectVariants -V {result_files_dict[r]} -select {select_condition} -O {result_filtered_files_dict[r]}'
        log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
        select_rare_logs.append(log)
        print("Counting selected", target_full_word[target], "variants for", r)
        count_variants(GATK, result_filtered_files_dict[r])

    return result_filtered_files_dict


# Selecting SNPs enriched in the analyzed cohort compared to population

def calc_binom_pval(row, p_col, target='r'):
    x = row['AC']
    n = row['AN']
    p = float(row[p_col])
    assert p <= 1.0 and p >= 0.0, "Population frequency must be between 0 and 1"
    assert target in ['r', 'c'], "Target must be 'r' or 'c'"
    if target == 'r':
        return binomtest(x, n, p, alternative='greater')
    else:
        return binomtest(x, n, p, alternative='less')


def select_enriched_snps(filtered_variants_files_dict, GATK, reference_population, bh_alpha=0.05, target='r'):
    # vcf file to csv table
    reference_population_cmd = f"gnomAD_genome_{reference_population}"
    totable_logs = []
    for r in ["promoter", "enhancer"]:
        command = f"{GATK} VariantsToTable  -V {filtered_variants_files_dict[r]} " \
                  "-F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN " \
                  f"-F gnomAD_genome_ALL -F {reference_population_cmd} " \
                  f"-O {filtered_variants_files_dict[r].replace('.vcf', '.csv')}"
        log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
        totable_logs.append(log)
    rare_promoter_snps = pd.read_csv(filtered_variants_files_dict['promoter'].replace('.vcf', '.csv'), sep='\t')
    rare_enhancer_snps = pd.read_csv(filtered_variants_files_dict['enhancer'].replace('.vcf', '.csv'), sep='\t')

    # execute one-sided binom test:
    # no. successes = no. ALT alleles in cohort
    # no. trials = no. all identified alleles
    # prob. of success = annotated popultaion frequency.
    # alt. hypetesis = observed freq is greated than expected

    for df in [rare_enhancer_snps, rare_promoter_snps]:
        df['binom_pval'] = df.apply(calc_binom_pval, axis=1, args=(reference_population_cmd,target,))

        # apply correction for multiple hypothesis testing with the Benjamini-Hochberg procedure, use FDR = bh_alpha
        multipletests_correction = multipletests(df['binom_pval'], alpha=bh_alpha,
                                                 method='fdr_bh', is_sorted=False, returnsorted=False)
        df['B-H_reject_H0'] = multipletests_correction[0]
        df['corrected_binom_pval'] = multipletests_correction[1]

        #round p-values to 3 decimal places
        df['binom_pval'] = df['binom_pval'].apply(lambda x: round(x, 4))
        df['corrected_binom_pval'] = df['corrected_binom_pval'].apply(lambda x: round(x, 4))

    rare_enriched_promoter_snps = rare_promoter_snps[rare_promoter_snps["B-H_reject_H0"]]
    rare_enriched_enhancer_snps = rare_enhancer_snps[rare_enhancer_snps["B-H_reject_H0"]]
    rare_enriched_promoter_snps = rare_enriched_promoter_snps.drop(columns='B-H_reject_H0')
    rare_enriched_enhancer_snps = rare_enriched_enhancer_snps.drop(columns='B-H_reject_H0')

    #left only one column with popultaion frequency
    cols_to_drop_enh = [col for col in rare_enriched_enhancer_snps.columns if ('gnomAD' in col) & (col !=reference_population_cmd)]
    cols_to_drop_prom = [col for col in rare_enriched_promoter_snps.columns if ('gnomAD' in col) & (col !=reference_population_cmd)]


    if cols_to_drop_enh and all(col in rare_enriched_enhancer_snps.columns for col in cols_to_drop_enh):
        rare_enriched_enhancer_snps = rare_enriched_enhancer_snps.drop(cols_to_drop_enh, axis=1)
    if cols_to_drop_prom and all(col in rare_enriched_promoter_snps.columns for col in cols_to_drop_prom):
        rare_enriched_promoter_snps = rare_enriched_promoter_snps.drop(cols_to_drop_prom, axis=1)

    print(len(rare_enriched_promoter_snps), "SNPs in promoters are enriched in analyzed cohort.")
    print(len(rare_enriched_enhancer_snps), "SNPs in enhancers are enriched in analyzed cohort.")
    return rare_enriched_promoter_snps, rare_enriched_enhancer_snps

# MOTIFS

def snps_to_bed_file(rare_enriched_promoter_snps, rare_enriched_enhancer_snps, OUTPUT):
    snps_bed_files = []
    for snps_df, r in [(rare_enriched_promoter_snps, "promoter"), (rare_enriched_enhancer_snps, "enhancer")]:
        snps_bed = pd.DataFrame()
        snps_bed["chromosome"] = snps_df["CHROM"]
        snps_bed["start"] = snps_df["POS"] - 1
        snps_bed["end"] = snps_df["POS"]
        snps_bed["name"] = snps_df["CHROM"] + ":" + snps_df["POS"].astype(str) + ":" + snps_df["REF"] + ":" + snps_df[
            "ALT"]
        snps_bed["score"] = 0
        snps_bed["strand"] = "+"
        output_bed_path = f"{OUTPUT}/{r}_rare_enriched_SNPs_intermediate_result.bed"
        snps_bed.to_csv(output_bed_path, sep="\t", index=False, header=False)
        snps_bed_files.append(output_bed_path)
    return snps_bed_files

def prepare_motifs_object():
    import rpy2.robjects as robjects
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)
    robjects.r('''
        library('motifbreakR')
        library('BSgenome.Hsapiens.UCSC.hg38')
        library('MotifDb')
        motifs <- query(MotifDb, andStrings=c("hocomocov11", "hsapiens"))    # czy to jest ok?
    ''')


def score_motifs(snps_bed_files, genomic_version):
    import rpy2.robjects as robjects
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)
    genome = {'hg38': 'BSgenome.Hsapiens.UCSC.hg38', 'mm38': 'BSgenome.Mmusculus.UCSC.mm39', 'dm3': 'BSgenome.Dmelanogaster.UCSC.dm3'}[genomic_version]
    robjects.r(f'''
        library('motifbreakR')
        library('{genome}')
        library('MotifDb')
        score_snps <- function(snps_file, out_file) {{
            #read SNPs from input bed file
            snps.mb.frombed <- snps.from.file(file = snps_file, search.genome = {genome}, format = "bed")
            #calculate scores
            results_log <- motifbreakR(snpList = snps.mb.frombed, filterp = TRUE,
                                   pwmList = motifs,
                                   threshold = 1e-5,
                                   method = "log",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = BiocParallel::bpparam())
        
            #reformat results to dataframe and save to file
            results_log_df <- data.frame(results_log)
            results_log_df <- apply(results_log_df,2,as.character)
            write.table(results_log_df, out_file, quote=F, sep="\t", row.names=F)
        }} 
    ''')

    score_snps_py = robjects.globalenv['score_snps']
    for snp_bed in snps_bed_files:
        snp_scores_csv = snp_bed.replace(".bed", "_motifbreakR-scores.csv")
        print(f"Calculating motif scores for {'promoters' if 'promoter' in snp_bed else 'enhancers'}.")
        snp_scores = score_snps_py(snp_bed, snp_scores_csv)
    return snp_scores


def find_best_matching_motif(group):
    # find motif with highest score (either for REF or ALT)
    scoreBest = "scoreRef" if max(group["scoreRef"])>max(group["scoreAlt"]) else "scoreAlt"
    best_score_motif = group[group[scoreBest] == max(group[scoreBest])]["providerId"].values[0]

    # find motif with highest abs(diff) between alleles
    best_alleleDiff = max(group["alleleDiff"].abs())
    best_alleleDiff_motif = group[group["alleleDiff"].abs() == best_alleleDiff]["providerId"].values[0]
    best_alleleDiff = group[group['providerId'] == best_alleleDiff_motif]["alleleDiff"].values[0]

    return f'{best_score_motif}:{round(max(group[scoreBest]),2)}', f'{best_alleleDiff_motif}:{round(best_alleleDiff,2)}'


def select_motif_results(rare_enriched_promoter_snps, rare_enriched_enhancer_snps, snps_bed_files):
    promoter_SNPs_motifbreakr = pd.read_csv(snps_bed_files[0].replace(".bed", "_motifbreakR-scores.csv"), sep="\t")
    enhancer_SNPs_motifbreakr = pd.read_csv(snps_bed_files[1].replace(".bed", "_motifbreakR-scores.csv"), sep="\t")

    # add columns with info about best matches: motif with highest score and alleleDiff
    # information about best motifs for each SNP will be stored in a dict in which SNP_ids will be keys
    best_motifs_dict = {}

    for df in [enhancer_SNPs_motifbreakr, promoter_SNPs_motifbreakr]:
        for snp_id, snp_records in df.groupby("SNP_id"):
            best_match, highest_diff = find_best_matching_motif(snp_records)
            best_motifs_dict[snp_id] = {"best_match": best_match, "highest_diff": highest_diff}
    

    # extract information from the dict to fill appropriate columns in enhancer and promoter SNPs dataframes
    enhancer_SNPs_motifbreakr.loc[:, "motif_best_match"] = enhancer_SNPs_motifbreakr.loc[:,"SNP_id"].apply(
                                                                    lambda x: best_motifs_dict[x]["best_match"])
    enhancer_SNPs_motifbreakr.loc[:, "motif_highest_diff"] = enhancer_SNPs_motifbreakr.loc[:,"SNP_id"].apply(
                                                                    lambda x: best_motifs_dict[x]["highest_diff"])

    promoter_SNPs_motifbreakr.loc[:, "motif_best_match"] = promoter_SNPs_motifbreakr.loc[:,"SNP_id"].apply(
                                                                    lambda x: best_motifs_dict[x]["best_match"])
    promoter_SNPs_motifbreakr.loc[:, "motif_highest_diff"] = promoter_SNPs_motifbreakr.loc[:,"SNP_id"].apply(
                                                                    lambda x: best_motifs_dict[x]["highest_diff"])

    # extract information about SNP location and best motifs, drop duplicates
    promoter_SNPs_motifbreakr_snps_only = promoter_SNPs_motifbreakr[
        ["seqnames", "start", "REF", "ALT", "motif_best_match", "motif_highest_diff"]].drop_duplicates()
    enhancer_SNPs_motifbreakr_snps_only = enhancer_SNPs_motifbreakr[
        ["seqnames", "start", "REF", "ALT", "motif_best_match", "motif_highest_diff"]].drop_duplicates()

    # change column names to keep the convention
    promoter_SNPs_motifbreakr_snps_only = promoter_SNPs_motifbreakr_snps_only.rename(
        columns={"seqnames": "CHROM", "start": "POS"})
    enhancer_SNPs_motifbreakr_snps_only = enhancer_SNPs_motifbreakr_snps_only.rename(
        columns={"seqnames": "CHROM", "start": "POS"})

    rare_enriched_promoter_snps_motif = pd.merge(rare_enriched_promoter_snps,
                                                 promoter_SNPs_motifbreakr_snps_only, how="left",
                                                 on=["CHROM", "POS", "REF", "ALT"]).replace(np.nan, '.')
    rare_enriched_enhancer_snps_motif = pd.merge(rare_enriched_enhancer_snps,
                                                 enhancer_SNPs_motifbreakr_snps_only, how="left",
                                                 on=["CHROM", "POS", "REF", "ALT"]).replace(np.nan, '.')
    print('Done: finding motifs influenced by SNPs.')
    return rare_enriched_promoter_snps_motif, rare_enriched_enhancer_snps_motif



def find_motifs(rare_enriched_promoter_snps, rare_enriched_enhancer_snps, OUTPUT, genomic_version):
    snps_bed_files = snps_to_bed_file(rare_enriched_promoter_snps, rare_enriched_enhancer_snps, OUTPUT)
    prepare_motifs_object()
    score_motifs(snps_bed_files, genomic_version)
    rare_enriched_promoter_snps_motif, rare_enriched_enhancer_snps_motif = select_motif_results(rare_enriched_promoter_snps, rare_enriched_enhancer_snps, snps_bed_files)

    return rare_enriched_promoter_snps_motif, rare_enriched_enhancer_snps_motif


def assign_genes_to_promoter_snps(rare_enriched_promoter_snps, PROMOTER_REGIONS_BED_PATH):
    rare_enriched_promoter_snps["genomic element"] = "promoter"

    # create BedTool object from promoter regions bed
    promoters_info = pbt.BedTool(PROMOTER_REGIONS_BED_PATH)

    # create BedTool object from dataframe with selected promoter SNPs
    rare_enriched_promoter_snps["POS-1"] = rare_enriched_promoter_snps["POS"] - 1
    rare_enriched_promoter_snps_bedtool = pbt.BedTool.from_dataframe(
        rare_enriched_promoter_snps[["CHROM", "POS-1", "POS"]])
    rare_enriched_promoter_snps = rare_enriched_promoter_snps.drop(labels=["POS-1"], axis=1)

    # intersect promoters and SNPs
    rare_enriched_promoter_snps_intersection = rare_enriched_promoter_snps_bedtool.intersect(promoters_info,
                                                                                                         wa=True,
                                                                                                         wb=True)

    # create a dataframe from the intersection results, keep only columns with SNP location and gene(s) name(s)
    rare_enriched_promoter_snps_intersection_df = rare_enriched_promoter_snps_intersection.to_dataframe(
        names=["CHROM", "POS","prom_start","prom_end", "Transcript"], usecols=[0, 2, 4,5,6]).drop_duplicates()
    rare_enriched_promoter_snps_gene = pd.merge(rare_enriched_promoter_snps,
                                                      rare_enriched_promoter_snps_intersection_df, how="left",
                                                      on=["CHROM", "POS"])
    rare_enriched_promoter_snps_gene['relation'] = 'containing'
    print('Done: assigning genes to promoters')
    return rare_enriched_promoter_snps_gene

# assign contacting transcripts using chromatin loops list
def assign_chromatin_contacting_gene_with_loops(snps, genes_info, CHROMATIN_LOOPS,OUTPUT, region):
    short_name = 'prom' if region=='promoter' else 'enh'

    looplist = pd.read_csv(CHROMATIN_LOOPS, sep='\t')
    assert all(column in looplist.columns for column in ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']) , 'Columns in looplist are not correct'
    # add prefix 'chr' to chromosome names if missing
    if looplist['chr1'][0][:3] != 'chr':
        looplist['chr1'] = looplist['chr1'].apply(lambda x: 'chr'+x)
        looplist['chr2'] = looplist['chr2'].apply(lambda x: 'chr'+x)
    looplist = looplist[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']].copy()

    distinct_snps = snps[['CHROM','POS','REF','ALT',f'{short_name}_start',f'{short_name}_end']].drop_duplicates()
    anchors = pd.concat([looplist[['chr1', 'x1', 'x2']].rename(columns={'x1':'anchor_start', 'x2':'anchor_end'}), looplist[['chr2', 'y1', 'y2']].rename(columns={'chr2':'chr1', 'y1':'anchor_start', 'y2':'anchor_end'})], ignore_index=True)
    anchors = anchors.drop_duplicates()

    # find regions which are in contact with transcript
    # it means that there is loop which one anchor intersects with regulatory region and second with transcript

    # Convert the DataFrame to BED format and save it as a file
    region_bed = distinct_snps[['CHROM', f'{short_name}_start', f'{short_name}_end']]
    region_bed.to_csv(OUTPUT+'/regulatory_region_intermediate_result.bed', sep='\t', header=False, index=False)

    anchor_bed = anchors[['chr1', 'anchor_start', 'anchor_end']]
    anchor_bed.to_csv(OUTPUT+'/anchor_intermediate_result.bed', sep='\t', header=False, index=False)

    genes_bed = genes_info[['chr', 'start', 'end']]
    genes_bed.to_csv(OUTPUT+'/genes_intermediate_result.bed', sep='\t', header=False, index=False)


    # Create BedTool objects
    region_bdt = pbt.BedTool(OUTPUT+'/regulatory_region_intermediate_result.bed')
    anchor_bdt = pbt.BedTool(OUTPUT+'/anchor_intermediate_result.bed')
    transcript_bdt = pbt.BedTool(OUTPUT+'/genes_intermediate_result.bed')

    # Intersect anchor with regulatory region
    anchor_region_intersection = anchor_bdt.intersect(region_bdt, wa=True, wb=True)

    # Intersect anchor with transcripts
    anchor_transcripts_intersection = anchor_bdt.intersect(transcript_bdt, wa=True, wb=True)

    # Convert the intersections to DataFrames
    anchor_region_df = anchor_region_intersection.to_dataframe(names=['chr1', 'anchor_start', 'anchor_end', 'CHROM', f'{short_name}_start', f'{short_name}_end'])
    anchor_transcripts_df = anchor_transcripts_intersection.to_dataframe(names=['chr1', 'anchor_start', 'anchor_end', 'chr', 'transcript_start', 'transcript_end'])

    # Merge regions that intersect with at least one anchor
    ll_region = looplist.merge(anchor_region_df, right_on=['chr1', 'anchor_start', 'anchor_end'],left_on=['chr1', 'x1', 'x2'], how='left', suffixes=('', '_left'))
    ll_region = ll_region.merge(anchor_region_df, right_on=['chr1', 'anchor_start', 'anchor_end'],left_on=['chr2', 'y1', 'y2'], how='left', suffixes=('', '_right'))
    ll_region = ll_region[ll_region['CHROM'].notna() | ll_region['CHROM_right'].notna()]

    # Merge transcripts that intersect with at least one anchor
    ll_region_transc = ll_region.merge(anchor_transcripts_df, right_on=['chr1', 'anchor_start', 'anchor_end'],left_on=['chr1', 'x1', 'x2'], how='left', suffixes=('', '_left'))
    ll_region_transc = ll_region_transc.merge(anchor_transcripts_df, right_on=['chr1', 'anchor_start', 'anchor_end'],left_on=['chr2', 'y1', 'y2'], how='left', suffixes=('', '_right'))
    ll_region_transc = ll_region_transc[ll_region_transc['chr'].notna() | ll_region_transc['chr_right'].notna()]

    col_to_drop = [col for col in ll_region_transc.columns if 'anchor' in col]
    ll_region_transc = ll_region_transc.drop(columns=col_to_drop)

    # Select pairs region-transcript where one anchor intersects with region and second with transcript
    ll_region_transc_left_right = ll_region_transc[(ll_region_transc[f'{short_name}_start'].notna() & ll_region_transc['transcript_start_right'].notna())]
    ll_region_transc_left_right = ll_region_transc_left_right[['chr1', f'{short_name}_start', f'{short_name}_end', 'chr2', 'transcript_start_right', 'transcript_end_right']]
    ll_region_transc_left_right = ll_region_transc_left_right.rename(columns={'transcript_start_right': 'transcript_start', 'transcript_end_right': 'transcript_end','chr1':f'{short_name}_chr','chr2':'transcript_chr'})
    ll_region_transc_right_left = ll_region_transc[(ll_region_transc[f'{short_name}_start_right'].notna() & ll_region_transc['transcript_start'].notna())]
    ll_region_transc_right_left = ll_region_transc_right_left[['chr2', f'{short_name}_start_right', f'{short_name}_end_right', 'chr1', 'transcript_start', 'transcript_end']]
    ll_region_transc_right_left = ll_region_transc_right_left.rename(columns={f'{short_name}_start_right': f'{short_name}_start', f'{short_name}_end_right': f'{short_name}_end','chr2':f'{short_name}_chr','chr1':'transcript_chr'})  
    found_pairs = pd.concat([ll_region_transc_left_right, ll_region_transc_right_left], ignore_index=True)
    found_pairs = found_pairs.drop_duplicates()

    found_pairs = found_pairs.merge(genes_info[['chr', 'start','end','Transcript']], left_on=['transcript_chr', 'transcript_start', 'transcript_end'], right_on=['chr', 'start', 'end'], how='left')
    found_pairs = found_pairs.rename(columns={'Transcript':'contacting_transcript'})

    # Merged found contacting genes with SNPs
    rare_enriched_region_snps_gene_chromatin = snps.merge(found_pairs[[f'{short_name}_chr',f'{short_name}_start',f'{short_name}_end','contacting_transcript']], left_on=['CHROM', f'{short_name}_start',f'{short_name}_end'], right_on=[f'{short_name}_chr', f'{short_name}_start',f'{short_name}_end'], how='inner')
    rare_enriched_region_snps_gene_chromatin = rare_enriched_region_snps_gene_chromatin.drop(columns=[f'{short_name}_chr','Transcript'])
    rare_enriched_region_snps_gene_chromatin = rare_enriched_region_snps_gene_chromatin.rename(columns={'contacting_transcript':'Transcript'})
    rare_enriched_region_snps_gene_chromatin['relation'] = 'contacting'
    rare_enriched_promoter_snps_gene_final = pd.concat([snps, rare_enriched_region_snps_gene_chromatin], ignore_index=True)
    return rare_enriched_promoter_snps_gene_final


# assign genes to snps which are located inside intronic enhancers 
def assign_genes_intronic_enhancer_snps(rare_enriched_enhancer_snps_df, ENHANCER_REGIONS, genes_info):
    enhancers = pd.read_csv(ENHANCER_REGIONS, sep='\t', header=None)
    enhancers.rename(columns={0: "enh_chr", 1: "enh_start", 2: "enh_end"}, inplace=True)

    # Intersect information about enhancers with SNPs to assign gene names to SNPs.
    # prepare bedtool objects
    rare_enriched_enhancer_snps_df["POS-1"] = rare_enriched_enhancer_snps_df["POS"] - 1
    gnomAD_col = [col for col in rare_enriched_enhancer_snps_df.columns if 'gnomAD' in col] # it has one element
    rare_enriched_enhancer_snps_bedtool = pbt.BedTool.from_dataframe(
        rare_enriched_enhancer_snps_df[["CHROM", "POS-1", "POS"]])
    enh_bedtool = pbt.BedTool.from_dataframe(enhancers)
    # intersect 
    rare_enriched_enhancer_snps_intersection = rare_enriched_enhancer_snps_bedtool.intersect(
        enh_bedtool, wa=True, wb=True, loj=True)
    # reformat intersection to dataframe, keep columns with enhancer coordinates - they will be usefull in the next step
    #dropping POS-1 column and enh_chr
    rare_enriched_enhancer_snps_gene = rare_enriched_enhancer_snps_intersection.to_dataframe(
        usecols=[0, 2, 4, 5],
        names=["CHROM", "POS", "enh_start", "enh_end"]).drop_duplicates()
    rare_enriched_enhancer_snps_gene = pd.merge(rare_enriched_enhancer_snps_df, rare_enriched_enhancer_snps_gene, how="left",
                                                on=["CHROM", "POS"]).replace(np.nan, '.')
    
    # assign transcript names by intersecting enhancer regions with transcripts
    transcripts_bedtool = pbt.BedTool.from_dataframe(genes_info[['chr', 'start', 'end', 'Transcript']])
    rare_enriched_enhancer_snps_gene_bedtool = pbt.BedTool.from_dataframe(rare_enriched_enhancer_snps_gene[['CHROM', 'enh_start', 'enh_end']])
    rare_enriched_enhancer_snps_gene_transcripts = rare_enriched_enhancer_snps_gene_bedtool.intersect(transcripts_bedtool, wa=True, wb=True)
    rare_enriched_enhancer_snps_gene_transcripts_df = rare_enriched_enhancer_snps_gene_transcripts.to_dataframe(
        names=['CHROM', 'enh_start', 'enh_end', 'Transcript'], usecols=[0, 1, 2, 6])
    rare_enriched_enhancer_snps_gene = pd.merge(rare_enriched_enhancer_snps_gene,
                                                rare_enriched_enhancer_snps_gene_transcripts_df, 
                                                how="left", on=["CHROM", "enh_start", "enh_end"]).replace(np.nan, '.')
    rare_enriched_enhancer_snps_gene.loc[rare_enriched_enhancer_snps_gene['Transcript'] != '.', 'genomic element'] = 'intronic enhancer'
    rare_enriched_enhancer_snps_gene.loc[rare_enriched_enhancer_snps_gene['Transcript'] == '.', 'genomic element'] = 'intergenic enhancer'
    return rare_enriched_enhancer_snps_gene


# assign genes to enhancers based on distance, save all in case of tie
def assign_closest_gene_to_enhancers(rare_enriched_enhancer_snps_gene, genes_info):
    genes_info_tss_bed = pbt.BedTool.from_dataframe(genes_info[['chr', 'tss', 'tss', 'Transcript']])
    genes_info_tss_bed_sorted = genes_info_tss_bed.sort()
    enhancers_bed = pbt.BedTool.from_dataframe(
        rare_enriched_enhancer_snps_gene[['CHROM', 'enh_start', 'enh_end']].drop_duplicates()).sort()

    tss_closest_to_enh = enhancers_bed.closest(genes_info_tss_bed_sorted, t='all', d=True)
    tss_closest_to_enh_df = tss_closest_to_enh.to_dataframe(
        names=['CHROM', 'enh_start', 'enh_end', "closest transcript"],
        usecols=[0, 1, 2, 6])
    
    rare_enriched_enhancer_snps_gene_closest = pd.merge(rare_enriched_enhancer_snps_gene,
                                                              tss_closest_to_enh_df, how="left",
                                                              on=["CHROM", "enh_start", "enh_end"])
    return rare_enriched_enhancer_snps_gene_closest

    
def reformat_target_genes_enh(rare_enriched_enhancer_snps_gene_closest_contacts, genes_info):
    genes_info['transcript_id'] = genes_info['transcript_id'].apply(lambda x: x.split('/')[0])
    snps_contacting = rare_enriched_enhancer_snps_gene_closest_contacts[rare_enriched_enhancer_snps_gene_closest_contacts['relation']=='contacting'].copy()
    snps_contacting.drop(columns=['closest transcript'], inplace=True)
    snps_contacting['Gene'] = snps_contacting['Transcript'].apply(lambda x: x.split('/')[1])
    snps_contacting['Transcript'] = snps_contacting['Transcript'].apply(lambda x: x.split('/')[0])
    snps_contacting['Gene_ID'] = snps_contacting.apply(lambda x: genes_info[genes_info['transcript_id'] == x['Transcript']]['Gene_ID'].values[0], axis=1)

    rare_enriched_enhancer_snps_gene_closest_contacts['Transcript'] = rare_enriched_enhancer_snps_gene_closest_contacts['Transcript'].apply(lambda x: x+'(containing)' if x!='.' else x)
    rare_enriched_enhancer_snps_gene_closest_contacts['closest transcript'] = rare_enriched_enhancer_snps_gene_closest_contacts['closest transcript'].apply(lambda x: x+'(closest)' if x!='.' else x)
    rare_enriched_enhancer_snps_gene_closest_contacts = rare_enriched_enhancer_snps_gene_closest_contacts[rare_enriched_enhancer_snps_gene_closest_contacts['relation']!='contacting']
    new_table=pd.DataFrame()
    new_rows = []
    for index, row in rare_enriched_enhancer_snps_gene_closest_contacts.iterrows():
        transcripts = [i for i in row[['Transcript','closest transcript']] if i!='.']
        for transcript in transcripts:
            new_row = row.copy()
            new_row = new_row.drop(['Transcript','closest transcript'])
            new_row['relation'] = transcript.split('(')[1].split(')')[0]
            new_row['Transcript'] = transcript.split('(')[0]
            new_row['Gene'] = new_row['Transcript'].split('/')[1]
            new_row['Transcript'] = new_row['Transcript'].split('/')[0]
            new_row['Gene_ID'] = genes_info[genes_info['transcript_id'] == new_row['Transcript']]['Gene_ID'].values[0]
            new_rows.append(new_row)
    new_table = pd.DataFrame(new_rows)
    new_table.drop_duplicates(inplace=True)
    new_table = new_table.groupby(new_table.columns.difference(['relation']).tolist(), as_index=False).agg({'relation': ';'.join})
    new_table['relation'] = new_table['relation'].apply(lambda x: ';'.join(list(set(x.split(';')))))
    new_table = pd.concat([new_table, snps_contacting], ignore_index=True)
    print('Done: assigning putative genes to enhancers')
    return new_table


def change_table_format_promoter(promoter_snps,genes_info):
    genes_info['transcript_id'] = genes_info['Transcript'].apply(lambda x: x.split('/')[0])

    promoter_snps['Gene'] = promoter_snps['Transcript'].apply(lambda x:x.split('/')[1])
    promoter_snps['Transcript'] = promoter_snps['Transcript'].apply(lambda x:x.split('/')[0])
    promoter_snps['Gene_ID'] = promoter_snps.apply(lambda x:genes_info[genes_info['transcript_id'] == x['Transcript']]['Gene_ID'].values[0], axis=1)
    promoter_snps = promoter_snps.groupby(promoter_snps.columns.difference(['relation']).tolist(), as_index=False).agg({'relation': ';'.join})
    promoter_snps['relation'] = promoter_snps['relation'].apply(lambda x: ';'.join(list(set(x.split(';')))))
    return promoter_snps

# Calculate correlation between enhancer activity and gene expression
def calculate_correlation_enh_gene(row, sample_names, counts, threshold):
    enh_act_vector = row[sample_names].values
    transcript = row['Transcript']
    gene_expr_rows = counts[counts["Transcript"] == (transcript+'_'+row['Gene'])]

    if len(gene_expr_rows) == 1:
        #calculate correlation for transcript
        expr_row = gene_expr_rows.iloc[0]
        expr_vector = expr_row[sample_names].values

        # Check if enh_act_vector and expr_vector has sufficient variability
        if np.std(enh_act_vector) != 0 and np.std(expr_vector) != 0:
            rho, pval = spearmanr(enh_act_vector, expr_vector)

            if str(rho) != 'nan' and rho > 0 and pval < threshold:
                return round(pval,3)
            
    elif len(gene_expr_rows) > 1:
        print('Transcript id in gene expression data is not unique: ',transcript)
    
    return "."


def check_signal_gene_expression_correlation_enhancer(rare_enriched_enhancer_snps_genes_collected,
                                                      ENHANCER_ACTIVITY, GENE_EXPRESSION, threshold):
    h3k27ac_cov = pd.read_csv(ENHANCER_ACTIVITY, sep="\t")
    h3k27ac_cov = h3k27ac_cov.rename(columns={"chr": "CHROM",
                                              "start": "enh_start",
                                              "end": "enh_end"})
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    # merge SNPs with H3K27ac coverage on enhancers
    rare_enriched_enhancer_snps_genes_collected_coverage = pd.merge(
        rare_enriched_enhancer_snps_genes_collected,
        h3k27ac_cov, on=["CHROM", "enh_start", "enh_end"], how="left")

    samples = h3k27ac_cov.columns.drop(["CHROM", "enh_start", "enh_end"])
    rare_enriched_enhancer_snps_genes_collected_coverage[
        "H3K27ac-expression correlation p-values"] = rare_enriched_enhancer_snps_genes_collected_coverage.apply(
        calculate_correlation_enh_gene, args=(samples, counts, threshold), axis=1)
    rare_enriched_enhancer_snps_genes_collected_corelations = rare_enriched_enhancer_snps_genes_collected_coverage.drop(
        labels=samples, axis=1)
    rare_enriched_enhancer_snps_genes_collected_corelations = rare_enriched_enhancer_snps_genes_collected_corelations[rare_enriched_enhancer_snps_genes_collected_corelations["H3K27ac-expression correlation p-values"] != "."].copy()
    print('Done: calculating correlation between enhancer activity and gene expression')
    return rare_enriched_enhancer_snps_genes_collected_corelations

#remove files with intermediate results, which constain string 'intermediate' in filename
def remove_unnecessary_files(OUTPUT):
    for file in os.listdir(OUTPUT):
        if 'intermediate' in file:
            os.remove(OUTPUT+'/'+file)


# add sample data to promoter and enhancer snps
def import_vcf_sample_level(INPUT_VCF, OUTPUT, GATK, promoter_snps, enhancer_snps):
    command = f'{GATK} VariantsToTable -V {INPUT_VCF}  \
        -F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN -GF AD -GF GT \
            -O {OUTPUT}/input_vcf_sample_level_intermediate_result.csv'
    log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                        universal_newlines=True).splitlines()
    all_snps = pd.read_csv(OUTPUT+'/input_vcf_sample_level_intermediate_result.csv', sep='\t')
    promoter_snps = promoter_snps.join(all_snps.set_index(['CHROM', 'POS', 'ALT', 'REF']), on=['CHROM', 'POS', 'ALT', 'REF'], how='left', rsuffix='_all_snps')
    enhancer_snps = enhancer_snps.join(all_snps.set_index(['CHROM', 'POS', 'ALT', 'REF']), on=['CHROM', 'POS', 'ALT', 'REF'], how='left', rsuffix='_all_snps')

    samples = [sample.split('.')[0] for sample in promoter_snps.columns if '.AD' in sample] 
    for sample in samples:
        promoter_snps[sample+'.var'] = promoter_snps.apply(lambda x: str(x[sample+'.GT']).split('/').count(x['ALT']) if x[sample+'.GT'] != './.' else './.', axis=1)
        enhancer_snps[sample+'.var'] = enhancer_snps.apply(lambda x: str(x[sample+'.GT']).split('/').count(x['ALT']) if x[sample+'.GT'] != './.' else './.', axis=1)

    var_columns = [sample+'.var' for sample in samples]

    dataframes = {'enhancer_snps': enhancer_snps, 'promoter_snps': promoter_snps}

    for name, df in dataframes.items():
        df['Num homref'] = (df[var_columns] == 0).sum(axis=1)
        df['Num het'] = (df[var_columns] == 1).sum(axis=1)
        df['Num homalt'] = (df[var_columns] == 2).sum(axis=1)

        #at least two genotype groups should be represented
        df = df[(df[['Num homref', 'Num het', 'Num homalt']] < 3).sum(axis=1) <= 1].copy()
        dataframes[name] = df

    enhancer_snps, promoter_snps = dataframes['enhancer_snps'], dataframes['promoter_snps']
    return promoter_snps, enhancer_snps
        
printed_warnings = set()
# calculate correlation between genotype and gene expression, returns string with results of correlation tests
def calculate_correlation_genotype_gene(row, counts):
    global printed_warnings

    samples_expr = counts.columns.drop(['Transcript', 'Gene'])
    samples_vcf = [col.split('.')[0] for col in row.index if '.var' in col]
    samples_vcf = [el for el in samples_vcf if row[el+'.var']!='./.']
    samples = list(set(samples_vcf) & set(samples_expr))    
    if len(samples) == 0:
        raise ValueError(f"No samples found for gene expression and genotype data. Check if samples names are correct.\n"
                 f"Gene expression samples: {samples_expr}\n"
                 f"Genotype samples: {samples_vcf}")
    transcript = row["Transcript"]

    if transcript == '.':
        return pd.Series(['.', '.'])
    transcript_counts = counts[counts['Transcript'] == row['Transcript']+'_'+row['Gene']]
    if len(transcript_counts) == 0:
        if transcript+'_'+row['Gene'] not in printed_warnings:
            print('Gene expression not found for transcript ', transcript+'_'+row['Gene'])
            printed_warnings.add(transcript+'_'+row['Gene'])
        return pd.Series(['.', '.'])

    expression_vector = transcript_counts[samples].values[0]

    genotype_vector = row[[sample+'.var' for sample in samples]]
    if np.std(expression_vector) != 0: 
        rho, pval = spearmanr(genotype_vector, expression_vector)
        if str(rho) != 'nan':
            return pd.Series([round(pval,3), rho>0])
    return pd.Series(['.', '.'])


def check_gene_genotype_correlation(GENE_EXPRESSION, promoter_snps, enhancer_snps):
    # Genes/transcripts normalized counts
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    counts['Gene'] = counts['Transcript'].apply(lambda x: x.split('_')[1])
    
    #calculate correlation between genotype and gene expression
    enhancer_snps[['gene_expr_correlations_pval','gene_expr_correlations_sign']] = enhancer_snps.apply(calculate_correlation_genotype_gene, args=(counts,), axis=1)
    promoter_snps[['gene_expr_correlations_pval','gene_expr_correlations_sign']] = promoter_snps.apply(calculate_correlation_genotype_gene, args=(counts,), axis=1)

    print('Done: calculating correlation between genotype and gene expression')
    return promoter_snps, enhancer_snps

def calculate_correlation_genotype_signal(row,samples_act):
    samples_var = [col.split('.')[0] for col in row.index if '.var' in col]
    samples_var = [el for el in samples_var if row[el+'.var']!='./.']
    samples = list(set(samples_var) & set(samples_act))

    if len(samples) == 0:
        raise ValueError(f"No samples found for enhancer activity and genotype data. Check if samples names are correct.\n"
                    f"Enhancer activity samples: {samples_act}\n"
                    f"Genotype samples: {samples_var}")
    
    enh_act_vector = row[samples].values
    genotype_vector = row[[sample+'.var' for sample in samples]]
    genotype_vector = [float(x) for x in genotype_vector]

    if len([i for i in [0,1,2] if genotype_vector.count(i)>2])<2:
        return pd.Series(['.', '.'])
        
    if np.std(enh_act_vector) != 0 and np.std(genotype_vector) != 0:
        rho, pval = spearmanr(genotype_vector, enh_act_vector)
        return pd.Series([round(pval,3), rho>0])
    else:
        return pd.Series(['.', '.'])
        

def check_genotype_signal_correlation(enhancer_snps, promoter_snps, ENHANCER_ACTIVITY, PROMOTER_ACTIVITY):
    # Calculations for enhancers
    h3k27ac_enh = pd.read_csv(ENHANCER_ACTIVITY, sep="\t")
    h3k27ac_enh = h3k27ac_enh.rename(columns={"chr": "CHROM",
                                              "start": "enh_start",
                                              "end": "enh_end"})
    
    # merge SNPs with H3K27ac coverage on enhancers
    samples_act = list(h3k27ac_enh.columns.drop(["CHROM", "enh_start", "enh_end"]))
    samples_var = [sample for sample in enhancer_snps.columns if '.var' in sample]
    merged = pd.merge(enhancer_snps, h3k27ac_enh, on=["CHROM", "enh_start", "enh_end"], how="left")
    merged = merged[['CHROM','enh_start','enh_end','POS']+samples_act+samples_var]
    merged.set_index(enhancer_snps.index, inplace=True)
    enhancer_snps[['genotype_act_corr_pval','genotype_act_corr_sign']] = merged.apply(calculate_correlation_genotype_signal, args=(samples_act,), axis=1)
    print('Done: calculating correlation between genotype and enhancer activity')

    # Calculations for promoters
    h3k27ac_prom = pd.read_csv(PROMOTER_ACTIVITY, sep="\t")
    h3k27ac_prom = h3k27ac_prom.rename(columns={"ID": "Transcript"})
    samples_act = list(h3k27ac_prom.columns.drop(["Transcript"]))
    samples_var = [sample for sample in promoter_snps.columns if '.var' in sample]
    merged = pd.merge(promoter_snps, h3k27ac_prom, on=["Transcript"], how="left")
    merged = merged[['CHROM','POS']+samples_act+samples_var]
    merged.set_index(promoter_snps.index, inplace=True)
    promoter_snps[['genotype_act_corr_pval','genotype_act_corr_sign']]= merged.apply(calculate_correlation_genotype_signal, args=(samples_act,), axis=1)

    print('Done: calculating correlation between genotype and promoter h3k27ac signal')

    return enhancer_snps, promoter_snps


def visualize_results(promoter_snps, enhancer_snps, GENE_EXPRESSION, OUTPUT,ENHANCER_ACTIVITY,PROMOTER_ACTIVITY,threshold,motifs):
    # select significant results
    enhancer_snps = enhancer_snps[(enhancer_snps['genotype_act_corr_pval'] != '.') | (enhancer_snps['gene_expr_correlations_pval'] != '.')].copy()
    promoter_snps = promoter_snps[(promoter_snps['genotype_act_corr_pval'] != '.') | (promoter_snps['gene_expr_correlations_pval'] != '.')].copy()

    promoter_snps = promoter_snps[(promoter_snps['gene_expr_correlations_pval'] == '.') | (promoter_snps['gene_expr_correlations_pval'].replace('.', 1).astype(float) < threshold)]
    enhancer_snps = enhancer_snps[(enhancer_snps['gene_expr_correlations_pval'] == '.') | (enhancer_snps['gene_expr_correlations_pval'].replace('.', 1).astype(float) < threshold)]
    promoter_snps = promoter_snps[(promoter_snps['genotype_act_corr_pval'] == '.') | (promoter_snps['genotype_act_corr_pval'].replace('.', 1).astype(float) < threshold)]
    enhancer_snps = enhancer_snps[(enhancer_snps['genotype_act_corr_pval'] == '.') | (enhancer_snps['genotype_act_corr_pval'].replace('.', 1).astype(float) < threshold)]

    # select results where sign of correlation is the same for gene expression and enhancer activity (only if signs for both tests are available)
    mask = ((enhancer_snps['gene_expr_correlations_pval'] == '.') |
        (enhancer_snps['genotype_act_corr_pval'] == '.') | 
        (enhancer_snps['genotype_act_corr_sign'] == enhancer_snps['gene_expr_correlations_sign']))

    enhancer_snps = enhancer_snps[mask]

    mask = ((promoter_snps['gene_expr_correlations_pval'] == '.') |
        (promoter_snps['genotype_act_corr_pval'] == '.') | 
        (promoter_snps['genotype_act_corr_sign'] == promoter_snps['gene_expr_correlations_sign']))
    promoter_snps = promoter_snps[mask]

    if len(enhancer_snps) == 0 and len(promoter_snps) == 0:
        print('No significant results found.')
        return

    promoter_snps = promoter_snps.sort_values(by = 'gene_expr_correlations_pval')
    enhancer_snps = enhancer_snps.sort_values(by = 'gene_expr_correlations_pval')
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    counts['Gene'] = counts['Transcript'].apply(lambda x: x.split('_')[1])
    h3k27ac_enh = pd.read_csv(ENHANCER_ACTIVITY, sep="\t")
    h3k27ac_enh = h3k27ac_enh.rename(columns={"chr": "CHROM",
                                              "start": "enh_start",
                                              "end": "enh_end"})
    h3k27ac_prom = pd.read_csv(PROMOTER_ACTIVITY, sep="\t")
    h3k27ac_prom = h3k27ac_prom.rename(columns={"ID": "Transcript"})


    enhancer_snps = enhancer_snps.drop_duplicates()
    promoter_snps = promoter_snps.drop_duplicates()

    #save plots to one pdf file
    with PdfPages(OUTPUT+'/regulatory_snps_plots.pdf') as pdf:  
        for index, row in enhancer_snps.iterrows():
            if (row['gene_expr_correlations_pval'] != '.') & (row['genotype_act_corr_pval'] != '.'):
                # plots for enhancers
                fig = plt.figure(figsize=(12, 5))

                #produce plot genotype vs. gene expression
                plt.subplot(1, 2, 1)
                samples_vcf = [col.split('.')[0] for col in enhancer_snps.columns if '.var' in col]
                samples_vcf = [el for el in samples_vcf if row[el+'.var']!='./.']
                samples_expr = counts.columns.drop(["Transcript",'Gene'])
                samples = list(set(samples_vcf) & set(samples_expr))
                vars = row[[sample + '.var' for sample in samples]]
                vars = [int(var)+np.random.uniform(-0.1, 0.1) for var in vars]
                transcript_gene  = row['Transcript']+'_'+row['Gene']
                exprs = counts[counts["Transcript"]==transcript_gene][samples]
                plt.scatter(vars, exprs, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('Expression for transcript '+row['Transcript']+'/'+row['Gene'])
                genotype_counts = str(row['Num homref'])+'/'+str(row['Num het'])+'/'+str(row['Num homalt'])
                pval = row['gene_expr_correlations_pval']
                plt.title('Genotype vs. gene_expression pvalue '+ str(pval)+' '+genotype_counts)

                #produce plot genotype vs. H3K27ac coverage
                plt.subplot(1, 2, 2)
                samples_act = h3k27ac_enh.columns.drop(["CHROM", "enh_start", "enh_end"])
                samples = list(set(samples_vcf) & set(samples_act))
                vars = row[[sample + '.var' for sample in samples]]
                vars = vars.values.tolist()
                vars = [int(var) for var in vars]
                genotype_counts = str(vars.count(0))+'/'+str(vars.count(1))+'/'+str(vars.count(2))
                vars = [int(var)+np.random.uniform(-0.1, 0.1) for var in vars]
                acts = h3k27ac_enh[(h3k27ac_enh["CHROM"]==row['CHROM']) & (h3k27ac_enh['enh_start']==row['enh_start']) & (h3k27ac_enh['enh_end']==row['enh_end'])][samples]
                plt.scatter(vars, acts, s=20, marker='o')                
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('H3K27ac coverage for enhancer '+row['CHROM']+':'+str(row['enh_start'])+'-'+str(row['enh_end']))
                pval = row['genotype_act_corr_pval']
                plt.title('Genotype vs. enhancer activity pvalue '+ str(pval)+' '+genotype_counts)

                fig.suptitle('SNP in enhancer '+row['CHROM']+':'+str(row['POS'])+row['REF']+'>'+row['ALT'])
                pdf.savefig() 
                plt.close()
            elif row['gene_expr_correlations_pval'] != '.':
                # plots for enhancers
                fig = plt.figure(figsize=(6, 5))

                #produce only one plot: genotype vs. gene expression
                samples_vcf = [col.split('.')[0] for col in enhancer_snps.columns if '.var' in col]
                samples_vcf = [el for el in samples_vcf if row[el+'.var']!='./.']
                samples_expr = counts
                samples = list(set(samples_vcf) & set(samples_expr))
                vars = row[[sample + '.var' for sample in samples]]
                vars = [int(var)+np.random.uniform(-0.1, 0.1) for var in vars]
                transcript_gene  = row['Transcript']+'_'+row['Gene']
                exprs = counts[counts['Transcript']==transcript_gene][samples]
                plt.scatter(vars, exprs, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('Expression for transcript '+row['Transcript']+'/'+row['Gene'])
                genotype_counts = str(row['Num homref'])+'/'+str(row['Num het'])+'/'+str(row['Num homalt'])
                pval = row['gene_expr_correlations_pval']
                plt.title('SNP in enhancer '+row['CHROM']+':'+str(row['POS'])+row['REF']+'>'+row['ALT']+'\nGenotype vs. gene_expression pvalue '+ str(pval)+' '+genotype_counts)
                pdf.savefig() 
                plt.close()
            else:
                fig = plt.figure(figsize=(6, 5))

                #produce only one plot: genotype vs. H3K27ac coverage
                samples_act = h3k27ac_enh.columns.drop(["CHROM", "enh_start", "enh_end"])
                samples = list(set(samples_vcf) & set(samples_act))
                vars = row[[sample + '.var' for sample in samples]]
                vars = vars.values.tolist()
                vars = [int(var) for var in vars]
                genotype_counts = str(vars.count(0))+'/'+str(vars.count(1))+'/'+str(vars.count(2))
                vars = [int(var)+np.random.uniform(-0.1, 0.1) for var in vars]
                acts = h3k27ac_enh[(h3k27ac_enh["CHROM"]==row['CHROM']) & (h3k27ac_enh['enh_start']==row['enh_start']) & (h3k27ac_enh['enh_end']==row['enh_end'])][samples]
                plt.scatter(vars, acts, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('H3K27ac coverage for enhancer '+row['CHROM']+':'+str(row['enh_start'])+'-'+str(row['enh_end']))
                pval = row['genotype_act_corr_pval']
                plt.title('SNP in enhancer '+row['CHROM']+':'+str(row['POS'])+row['REF']+'>'+row['ALT']+'\nGenotype vs. enhancer activity pvalue '+ str(pval)+' '+genotype_counts)
                pdf.savefig()
                plt.close()


        for index, row in promoter_snps.iterrows():
            if (row['gene_expr_correlations_pval'] != '.') & (row['genotype_act_corr_pval'] != '.'):
    
                # plots for promoters
                fig = plt.figure(figsize=(12, 5))

                #produce plot genotype vs. gene expression
                plt.subplot(1, 2, 1)
                samples_vcf = [col.split('.')[0] for col in promoter_snps.columns if '.var' in col]
                samples_vcf = [el for el in samples_vcf if row[el+'.var']!='./.']
                samples_expr = counts.columns.drop(["Transcript",'Gene'])
                samples = list(set(samples_vcf) & set(samples_expr))
                vars = row[[sample + '.var' for sample in samples]]
                vars = [int(var)+np.random.uniform(-0.1, 0.1) for var in vars]
                transcript_gene  = row['Transcript']+'_'+row['Gene']
                exprs = counts[counts["Transcript"]==transcript_gene][samples]
                plt.scatter(vars, exprs, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('Expression for transcript '+row['Transcript']+'/'+row['Gene'])
                genotype_counts = str(row['Num homref'])+'/'+str(row['Num het'])+'/'+str(row['Num homalt'])
                pval = row['gene_expr_correlations_pval']
                plt.title('Genotype vs. gene_expression pvalue '+ str(pval)+' '+genotype_counts)

                #produce plot genotype vs. H3K27ac coverage
                plt.subplot(1, 2, 2)
                samples_act = h3k27ac_prom.columns.drop(["Transcript"])
                samples = list(set(samples_vcf) & set(samples_act))
                vars = row[[sample + '.var' for sample in samples]]
                vars = vars.values.tolist()
                vars = [int(var) for var in vars]
                genotype_counts = str(vars.count(0))+'/'+str(vars.count(1))+'/'+str(vars.count(2))
                vars = [var+np.random.uniform(-0.1, 0.1) for var in vars]
                acts = h3k27ac_prom[h3k27ac_prom["Transcript"]==row['Transcript']][samples]
                plt.scatter(vars, acts, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('H3K27ac coverage for promoter '+row['Transcript']+'/'+row['Gene'])
                pval = row['genotype_act_corr_pval']
                plt.title('Genotype vs. h3k27ac pvalue '+ str(pval)+' '+genotype_counts)

                fig.suptitle('SNP in promoter '+row['CHROM']+':'+str(row['POS'])+row['REF']+'>'+row['ALT'])
                pdf.savefig() 
                plt.close()
            elif row['gene_expr_correlations_pval'] != '.':
                # plots for promoters
                fig = plt.figure(figsize=(6, 5))

                #produce only one plot: genotype vs. gene expression
                samples_vcf = [col.split('.')[0] for col in promoter_snps.columns if '.var' in col]
                samples_vcf = [el for el in samples_vcf if row[el+'.var']!='./.']
                samples_expr = counts
                samples = list(set(samples_vcf) & set(samples_expr))
                vars = row[[sample + '.var' for sample in samples]]
                vars = [int(var)+np.random.uniform(-0.1, 0.1) for var in vars]
                transcript_gene  = row['Transcript']+'_'+row['Gene']
                exprs = counts[counts['Transcript']==transcript_gene][samples]
                plt.scatter(vars, exprs, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('Expression for transcript '+row['Transcript']+'/'+row['Gene'])
                genotype_counts = str(row['Num homref'])+'/'+str(row['Num het'])+'/'+str(row['Num homalt'])
                pval = row['gene_expr_correlations_pval']
                plt.title('SNP in promoter '+row['CHROM']+':'+str(row['POS'])+row['REF']+'>'+row['ALT']+'\nGenotype vs. gene_expression pvalue '+ str(pval)+' '+genotype_counts)
                pdf.savefig()
                plt.close()
            else:
                fig = plt.figure(figsize=(6, 5))

                #produce only one plot: genotype vs. H3K27ac coverage
                samples_act = h3k27ac_prom.columns.drop(["Transcript"])
                samples = list(set(samples_vcf) & set(samples_act))
                vars = row[[sample + '.var' for sample in samples]]
                vars = vars.values.tolist()
                vars = [int(var) for var in vars]
                genotype_counts = str(vars.count(0))+'/'+str(vars.count(1))+'/'+str(vars.count(2))
                vars = [var+np.random.uniform(-0.1, 0.1) for var in vars]
                acts = h3k27ac_prom[h3k27ac_prom["Transcript"]==row['Transcript']][samples]
                plt.scatter(vars, acts, s=20, marker='o')
                genotype = [row['REF']*2, row['REF']+row['ALT'], row['ALT']*2]
                plt.xticks([0, 1, 2], genotype)
                plt.xlabel('Genotype')
                plt.ylabel('H3K27ac coverage for promoter '+row['Transcript']+'/'+row['Gene'])
                pval = row['genotype_act_corr_pval']
                plt.title('SNP in promoter '+row['CHROM']+':'+str(row['POS'])+row['REF']+'>'+row['ALT']+'\nGenotype vs. h3k27ac pvalue '+ str(pval)+' '+genotype_counts)
                pdf.savefig()
                plt.close()
        
    # dropping not important columns
    cols_to_drop_enh = [col for col in enhancer_snps.columns if (('.AD' in col) or ('.GT' in col) or ('.var' in col))]
    cols_to_drop_prom = [col for col in promoter_snps.columns if (('.AD' in col) or ('.GT' in col) or ('.var' in col))]
    cols_to_drop_enh = cols_to_drop_enh +['enh_start', 'enh_end']
    cols_to_drop_enh = cols_to_drop_enh + [col for col in enhancer_snps.columns if ('all_snps' in col) or ('Unnamed' in col)]
    cols_to_drop_prom = cols_to_drop_prom + [col for col in promoter_snps.columns if ('all_snps' in col) or ('Unnamed' in col)]
    enhancer_snps = enhancer_snps.drop(columns = cols_to_drop_enh)
    promoter_snps = promoter_snps.drop(columns = cols_to_drop_prom)

    #format columns
    enhancer_snps = enhancer_snps.rename(columns = {'genotype_act_corr_pval' : 'p-value_for_correlation_genotype_vs_h3k27ac',
                                                    'gene_expr_correlations_pval':'p-value_for_correlation_genotype_vs_expression',
                                                    'H3K27ac-expression correlation p-values':'p-value_for_correlation_expression_vs_h3k27ac',
                                                    'Gene':'Gene_name', 'binom_pval': 'p-value_for_binomial_test', 
                                                    'corrected_binom_pval': 'corrected_p-value_for_binomial_test'})
    
    promoter_snps = promoter_snps.rename(columns = {'genotype_act_corr_pval' : 'p-value_for_correlation_genotype_vs_h3k27ac',
                                                    'gene_expr_correlations_pval':'p-value_for_correlation_genotype_vs_expression',
                                                    'Gene':'Gene_name', 'binom_pval': 'p-value_for_binomial_test', 
                                                    'corrected_binom_pval': 'corrected_p-value_for_binomial_test'})


    # change columns order
    gnomad_cols = [col for col in enhancer_snps if 'gnomAD' in col]
    motif_cols = ['motif_best_match','motif_highest_diff'] if motifs else []
    cols_order_enh = ['CHROM','POS','REF','ALT','AC','AF','AN','Num homref','Num het','Num homalt']+gnomad_cols+['genomic element','Gene_name','Gene_ID','Transcript','relation']+ [col for col in enhancer_snps.columns if ('p-value' in col)]+motif_cols
    enhancer_snps = enhancer_snps.reindex(columns=cols_order_enh)   

    cols_order_prom = ['CHROM','POS','REF','ALT','AC','AF','AN','Num homref','Num het','Num homalt']+gnomad_cols+['genomic element','Gene_name','Gene_ID','Transcript','relation']+ [col for col in promoter_snps.columns if ('p-value' in col)]+motif_cols
    promoter_snps = promoter_snps.reindex(columns=cols_order_prom)
    promoter_snps['p-value_for_correlation_expression_vs_h3k27ac'] = '.'
    promoter_snps.reset_index(drop=True, inplace=True)
    enhancer_snps.reset_index(drop=True, inplace=True)
    enhancer_snps = enhancer_snps.drop_duplicates()
    promoter_snps = promoter_snps.drop_duplicates()
    found_snps = pd.concat([enhancer_snps, promoter_snps], ignore_index=True)

    #saving results to csv
    enhancer_snps.to_csv(OUTPUT+'/final_regulatory_snps_enhancer.csv')
    promoter_snps.to_csv(OUTPUT+'/final_regulatory_snps_promoter.csv')
    found_snps.to_csv(OUTPUT+'/final_regulatory_snps.csv')
    print('Found regulatory variants are saved in files:\n', OUTPUT+'/final_regulatory_snps_enhancer.csv\n', OUTPUT+'/final_regulatory_snps_promoter.csv\n', 'Plots are saved in file:', OUTPUT+'/regulatory_snps_plots.pdf')


def save_limited_results(promoter_snps, enhancer_snps, OUTPUT, motifs):
    if len(promoter_snps) == 0 and len(enhancer_snps) == 0:
        print('No significant results found.')
        return
    # format columns
    gnomad_cols = [col for col in enhancer_snps if 'gnomAD' in col]
    pval_cols = [col for col in enhancer_snps if 'pval' in col]
    motif_cols = ['motif_best_match','motif_highest_diff'] if motifs else []
    enhancer_snps = enhancer_snps[['CHROM','POS','REF','ALT','AC','AF','AN']+gnomad_cols+['genomic element','Transcript','Gene','Gene_ID','relation']+pval_cols+motif_cols]
    promoter_snps = promoter_snps[['CHROM','POS','REF','ALT','AC','AF','AN']+gnomad_cols+['genomic element','Transcript','Gene','Gene_ID','relation']+pval_cols+motif_cols]
    enhancer_snps = enhancer_snps.rename(columns = {'Gene':'Gene_name', 'binom_pval': 'p-value_for_binomial_test', 'corrected_binom_pval': 'corrected_p-value_for_binomial_test'})
    promoter_snps = promoter_snps.rename(columns = {'Gene':'Gene_name', 'binom_pval': 'p-value_for_binomial_test', 'corrected_binom_pval': 'corrected_p-value_for_binomial_test'})

    promoter_snps.reset_index(drop=True, inplace=True)
    enhancer_snps.reset_index(drop=True, inplace=True)
    enhancer_snps = enhancer_snps.drop_duplicates()
    promoter_snps = promoter_snps.drop_duplicates()
    found_snps = pd.concat([enhancer_snps, promoter_snps], ignore_index=True)

    #save results to csv
    enhancer_snps.to_csv(OUTPUT+'/limited_final_regulatory_snps_enhancer.csv')
    promoter_snps.to_csv(OUTPUT+'/limited_final_regulatory_snps_promoter.csv')
    found_snps.to_csv(OUTPUT+'/limited_final_regulatory_snps.csv')
    print('Found regulatory variants are saved in files:\n', OUTPUT+'/limited_final_regulatory_snps_enhancer.csv\n', OUTPUT+'/limited_final_regulatory_snps_promoter.csv\n', OUTPUT+'/limited_final_regulatory_snps.csv')
