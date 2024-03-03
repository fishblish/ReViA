import pandas as pd
import numpy as np
from scipy.stats import binom_test, spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools as pbt
import os
import re
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
            print('ANNOVAR ran successfully, but returned a non-zero exit code (1). This can be ignored.')
            print('ANNOVAR is installed correctly.')
            flags[1] = 0

        else:
            print('ANNOVAR is not installed correctly in indicated localization.')
            print(f"Command '{''.join(cmd)}' failed with error code {e.returncode}:\n{e.output}")

    return flags

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


def select_biallelic_inside_regulatory_regions(GATK, INPUT_VCF, PROMOTER_REGIONS, ENHANCER_REGIONS, OUTPUT):
    select_logs = []
    count_logs = []
    result_files={}
    for r, regions in [("promoter", PROMOTER_REGIONS), ("enhancer", ENHANCER_REGIONS)]:
        result_files[r] = f'{OUTPUT}/{r}_intermediate_result.vcf'
        command1 = f"{GATK} SelectVariants -V {INPUT_VCF} -L {regions} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {result_files[r]}"
        print('Selecting biallelic variants for', r)
        log1 = subprocess.check_output(command1, shell=True, stderr=subprocess.STDOUT,
                                       universal_newlines=True).splitlines()
        print(f'Biallelic variants inside {r} regions are saved in file: {result_files[r]}')
        select_logs.append(str(log1))

        command2 = f"{GATK} CountVariants -V {result_files[r]}"
        print('Counting biallelic variants for', r)
        log2 = subprocess.check_output(command2, shell=True, stderr=subprocess.STDOUT,
                                       universal_newlines=True).splitlines()
        count_logs.append(log2)
        print(f"Number of biallelic SNPs in {r} regions:{log2[-1]} \n")
    logs = [select_logs, count_logs]
    return result_files


#
# #ANNOVAR - annotate frequencies
# #download database
def prepare_annovar_database(ANNOVAR):   # może dodać możliwość załączenia własnej bazy danych
    downloaded_flag=0
    if ('humandb' in os.listdir(ANNOVAR)):
        if ('hg38_gnomad_genome.txt' in os.listdir(ANNOVAR + 'humandb')):
            print('You already have human genome 38 database ready to use.')
            return 1
        if ('hg38_gnomad_genome.txt.gz' in os.listdir(ANNOVAR + 'humandb')):
            print('You already have human genome 38 database. It need to be unpacked. Do you want to unpack it now? (y/n)')
            answer = input()
            if answer == 'y':
                command = f'gunzip {ANNOVAR}humandb/hg38_gnomad_genome.txt.gz'
                try:
                    logs = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                                    universal_newlines=True).splitlines()
                    print('Annovar database human genome 38 unpacked under path:', ANNOVAR + 'humandb/')
                    downloaded_flag = 1
                except subprocess.CalledProcessError as e:
                    print('ANNOVAR database could not be unpacked. Please check if you have enough space on your disk.')
                    return 0
            

    if downloaded_flag == 0:
        command = f'perl {ANNOVAR}annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_genome {ANNOVAR}humandb/'
        try:
            print('Database need to be downloaded')
            logs = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                             universal_newlines=True).splitlines()
            print('Annovar database human genome 38 downloaded under path:', ANNOVAR + 'humandb/') #'hg38_gnomad_genome.txt' is 17G heavy
            downloaded_flag = 1
        except subprocess.CalledProcessError as e:
            print(f"Command '{command}' failed with error code {e.returncode}:\n{e.output}")
    return downloaded_flag # 1 - database downloaded, 0 - database not downloaded


def annotate_freq(variants_vcf_file_dict, ANNOVAR):
    annovar_logs = []
    result_vcf_files_dict={}
    for r in ["promoter", "enhancer"]:
        out_file = variants_vcf_file_dict[r].split('.vcf')[0]
        result_vcf_files_dict[r] = out_file+'.hg38_multianno.vcf'
        
        command = f"perl {ANNOVAR}table_annovar.pl {variants_vcf_file_dict[r]} {ANNOVAR}humandb/ -buildver hg38 -remove -protocol gnomad_genome -operation f -nastring . -vcfinput -out {out_file} -polish" 
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
    gaps_dict = {'r': '=0.0', 'c': '=100.0'}
    ineq_sign = {'r': '<', 'c': '>'}
    target_full_word = {'r': 'rare', 'c': 'common'}
    population = list(set(population))
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
        print("Done: selecting", target_full_word[target], "snps for", r)
        print("Counting selected", target_full_word[target], "variants for", r)
        count_variants(GATK, result_filtered_files_dict[r])

    return result_filtered_files_dict


# Selecting SNPs enriched in the analyzed cohort compared to population

def calc_binom_pval(row, p_col='gnomAD_genome_NFE', target='r'):
    x = row['AC']
    n = row['AN']
    p = float(row[p_col])
    assert p <= 1.0 and p >= 0.0, "Population frequency must be between 0 and 1"
    assert target in ['r', 'c'], "Target must be 'r' or 'c'"
    if target == 'r':
        return binom_test(x, n, p, alternative='greater')
    else:
        return binom_test(x, n, p, alternative='less')


def select_enriched_snps(filtered_variants_files_dict, GATK, reference_population, bh_alpha=0.01, target='r'):
    # vcf file to csv table
    reference_population_cmd = f"gnomAD_genome_{reference_population}"
    totable_logs = []
    for r in ["promoter", "enhancer"]:
        command = f"{GATK} VariantsToTable  -V {filtered_variants_files_dict[r]} " \
                  "-F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN " \
                  f"-F {reference_population_cmd} " \
                  f"-O {filtered_variants_files_dict[r].replace('.vcf', '.csv')}"
        log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
        totable_logs.append(log)
        print("Done: vcf file to table")
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
    rare_enriched_promoter_snps = rare_promoter_snps[rare_promoter_snps["B-H_reject_H0"]]
    rare_enriched_enhancer_snps = rare_enhancer_snps[rare_enhancer_snps["B-H_reject_H0"]]

    #left only one column with popultaion frequency
    cols_to_drop_enh = [col for col in rare_enriched_enhancer_snps.columns if 'gnomAD' in col].remove(reference_population_cmd)
    cols_to_drop_prom = [col for col in rare_enriched_promoter_snps.columns if 'gnomAD' in col].remove(reference_population_cmd)
    if cols_to_drop_enh and all(col in rare_enriched_enhancer_snps.columns for col in cols_to_drop_enh):
        rare_enriched_enhancer_snps = rare_enriched_enhancer_snps.drop(cols_to_drop_enh, axis=1)
    if cols_to_drop_prom and all(col in rare_enriched_promoter_snps.columns for col in cols_to_drop_prom):
        rare_enriched_promoter_snps = rare_enriched_promoter_snps.drop(cols_to_drop_prom, axis=1)

    print(len(rare_enriched_promoter_snps), "SNPs in promoters are enriched in analyzed cohort.")
    print(len(rare_enriched_enhancer_snps), "SNPs in enhancers are enriched in analyzed cohort.")
    return rare_enriched_promoter_snps, rare_enriched_enhancer_snps



def assign_genes_to_promoter_snps(rare_enriched_promoter_snps, PROMOTER_REGIONS):
    rare_enriched_promoter_snps["genomic element"] = "promoter"
    # create BedTool object from promoter regions bed
    promoters_info = pbt.BedTool(PROMOTER_REGIONS)

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
        names=["CHROM", "POS", "Gene"], usecols=[0, 2, 6]).drop_duplicates()
    rare_enriched_promoter_snps_gene = pd.merge(rare_enriched_promoter_snps,
                                                      rare_enriched_promoter_snps_intersection_df, how="left",
                                                      on=["CHROM", "POS"])
    return rare_enriched_promoter_snps_gene

# assign genes to snps which are located inside intronic enhancers 
def assign_genes_intronic_enhancer_snps(rare_enriched_enhancer_snps_df, ENHANCER_REGIONS):
    enh_genes = pd.read_csv(ENHANCER_REGIONS, sep='\t', names=['chr', 'start', 'end', 'Gene'])

    # reformat to have one gene ID in cell
    enh_one_gene = pd.DataFrame()
    for i, row in enh_genes.iterrows():
        genes = row['Gene'].split(',')
        if len(genes) == 1:
            enh_one_gene = pd.concat([enh_one_gene, pd.DataFrame([row])], ignore_index=True)
        else:
            for gene in genes:
                new_row = row
                new_row['Gene'] = gene
                enh_one_gene = pd.concat([enh_one_gene, pd.DataFrame([new_row])], ignore_index=True)
                enh_one_gene = enh_one_gene.reindex(enh_genes.columns, axis=1)
                enh_one_gene["start"] = enh_one_gene["start"].astype(int)
                enh_one_gene["end"] = enh_one_gene["end"].astype(int)

    # Intersect information about enhancers with SNPs to assign gene names to SNPs.
    # prepare bedtool objects
    rare_enriched_enhancer_snps_df["POS-1"] = rare_enriched_enhancer_snps_df["POS"] - 1
    gnomAD_cols = [col for col in rare_enriched_enhancer_snps_df.columns if 'gnomAD' in col]
    rare_enriched_enhancer_snps_bedtool = pbt.BedTool.from_dataframe(
        rare_enriched_enhancer_snps_df[["CHROM", "POS-1", "POS", "REF", "ALT", "AC", "AF", "AN"]+gnomAD_cols+
                                       [ "binom_pval","B-H_reject_H0", "corrected_binom_pval"]])
    enh_genes_bedtool = pbt.BedTool.from_dataframe(enh_one_gene)

    # intersect 
    rare_enriched_enhancer_snps_intersection = rare_enriched_enhancer_snps_bedtool.intersect(
        enh_genes_bedtool, wa=True, wb=True, loj=True)
    # reformat intersection to dataframe, keep columns with enhancer coordinates - they will be usefull in the next step
    #dropping POS-1 column
    rare_enriched_enhancer_snps_gene = rare_enriched_enhancer_snps_intersection.to_dataframe(
        usecols=[0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15],
        names=["CHROM", "POS", "REF", "ALT", "AC", "AF", "AN"]+gnomAD_cols+
        ["binom_pval","B-H_reject_H0", "corrected_binom_pval", "enh_start", "enh_end", "Gene"])
    rare_enriched_enhancer_snps_gene["genomic element"] = rare_enriched_enhancer_snps_gene.Gene.apply(
        lambda x: "enhancer intergenic" if x == "." else "enhancer intronic")
    return rare_enriched_enhancer_snps_gene


def find_tss(row):
    if row['strand'] == '+':
        return row['start']
    else:
        return row['end']


def prepare_genes_info(GENES_INFO):
    rows_to_skip = 0
    with open(GENES_INFO) as f:
        for line in f:
            if "#" in line:
                rows_to_skip += 1
            else:
                break

    full_annot = pd.read_csv(GENES_INFO, sep='\t', skiprows=rows_to_skip, usecols=[0, 2, 3, 4, 6, 8], 
                            names = ["chr", "type", "start", "end", "strand", 'info'])
    genes_info = full_annot[full_annot["type"] == "gene"]
    #reformat genes_info['info'] to dictionary
    dict_series = genes_info['info'].str.split(';')
    #choose nonempty elements
    dict_series = pd.Series([[el for el in row if len(el) > 1] for row in dict_series])
    genes_info['info'] = dict_series.apply(lambda x: dict([el.split(' ')[-2:] for el in x]))
    genes_info['info'] = genes_info['info'].apply(lambda x: {k: v.strip('"') for k, v in x.items()})
    genes_info['ID'] = genes_info['info'].apply(lambda x: x['gene_id'])
    genes_info['Gene'] = genes_info['info'].apply(lambda x: x['gene_id']+'/'+x['gene_name'])
    genes_info["tss"] = genes_info.apply(find_tss, axis=1)
    return genes_info

# assign genes to enhancers based on distance, save all in case of tie
def assign_closest_gene_to_enhancers(rare_enriched_enhancer_snps_gene, genes_info):

    genes_info_tss_bed = pbt.BedTool.from_dataframe(genes_info[['chr', 'tss', 'tss', 'Gene']])
    genes_info_tss_bed_sorted = genes_info_tss_bed.sort()

    enhancers_bed = pbt.BedTool.from_dataframe(
        rare_enriched_enhancer_snps_gene[['CHROM', 'enh_start', 'enh_end']].drop_duplicates()).sort()

    tss_closest_to_enh = enhancers_bed.closest(genes_info_tss_bed_sorted, t='all', d=True)
    tss_closest_to_enh_df = tss_closest_to_enh.to_dataframe(
        names=['CHROM', 'enh_start', 'enh_end', "closest gene", "distance to closest gene"],
        usecols=[0, 1, 2, 6, 7])
    
    rare_enriched_enhancer_snps_gene_closest = pd.merge(rare_enriched_enhancer_snps_gene,
                                                              tss_closest_to_enh_df, how="left",
                                                              on=["CHROM", "enh_start", "enh_end"])
    return rare_enriched_enhancer_snps_gene_closest


def add_gene_name(gene_id, genes_info):
    if len(genes_info[genes_info["ID"] == gene_id]) != 0:
        return genes_info[genes_info["ID"] == gene_id]["Gene"].values[0]
    else:
        return '.'

# assign putative genes to enhancers based on chromatin contacts
def assign_chromatin_contacting_gene_to_enhancer(rare_enriched_enhancer_snps_gene_closest,
                                                              genes_info, CHROMATIN_CONTACTS):
    # Read bed file with chromatin contacts, merge with SNPs table
    contacts = pd.read_csv(CHROMATIN_CONTACTS, sep=' ')
    contacts_to_genes = contacts[contacts["ENSG"] != '-']
    contacts_to_genes = contacts_to_genes.drop_duplicates(subset=["chr", "start", "end", "ENSG"])
    contacts_to_genes = contacts_to_genes.rename(columns={"chr": "CHROM",
                                                          "start": "enh_start",
                                                          "end": "enh_end",
                                                          "ENSG": "contacting gene"})

    rare_enriched_enhancer_snps_gene_closest_contacts = pd.merge(rare_enriched_enhancer_snps_gene_closest,
                                                                       contacts_to_genes[
                                                                           ["CHROM", "enh_start", "enh_end",
                                                                            "contacting gene"]],
                                                                       on=["CHROM", "enh_start", "enh_end"],
                                                                       how="left").fillna('.')

    rare_enriched_enhancer_snps_gene_closest_contacts["contacting gene"] = \
    rare_enriched_enhancer_snps_gene_closest_contacts["contacting gene"].apply(
        lambda x: add_gene_name(x, genes_info))
    return rare_enriched_enhancer_snps_gene_closest_contacts

# collect putative genes in one column named 'Genes'
def reformat_target_genes_enh(rare_enriched_enhancer_snps_gene_closest_contacts):
    rare_enriched_enhancer_snps_genes_collected = pd.DataFrame()

    for name, group in rare_enriched_enhancer_snps_gene_closest_contacts.groupby(["CHROM", "POS", "REF", "ALT"]):
        containing_genes = [gene + "(containing)" for gene in group["Gene"].unique() if gene != "."]
        closest_genes = [gene + "(closest)" for gene in group["closest gene"].unique() if gene != "."]
        contacting_genes = [gene + "(contacting)" for gene in group["contacting gene"].unique() if gene != "."]

        all_genes = containing_genes + closest_genes + contacting_genes
        gnomAD_col = [col for col in group.columns if 'gnomAD' in col]
        group["Gene"] = ";".join(all_genes)
        row = group[['CHROM', 'POS', 'REF', 'ALT',
                     'AC', 'AF', 'AN']+gnomAD_col+['binom_pval', 'B-H_reject_H0',
                     'corrected_binom_pval', 'enh_start', 'enh_end', 'Gene', 
                     'genomic element']].drop_duplicates()

        rare_enriched_enhancer_snps_genes_collected = pd.concat(
            [rare_enriched_enhancer_snps_genes_collected, pd.DataFrame(row)], ignore_index=True)
    return rare_enriched_enhancer_snps_genes_collected

# Calculate correlation between enhancer activity and gene expression
def calculate_correlation_enh_gene(row, sample_names, GENE_EXPRESSION):
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    counts['Gene'] = counts.Transcript.apply(lambda x: '_'.join(x.split('_')[1:]))

    correlations = ""

    enh_act_vector = row[sample_names].values

    genes = set([el.split("(")[0] for el in row["Gene"].split(';')])

    # iterate over all target genes assigned to this variant
    for gene in genes:
        gene_name = gene.split('/')[1]
        gene_expr_rows = counts[counts["Gene"] == gene_name]

        if len(gene_expr_rows) != 0:
            #calculate correlations for each transcript of the analyzed gene
            gene_correlations = {}
            for j, expr_row in gene_expr_rows.iterrows():
                expr_vector = expr_row[sample_names].values

                # Check if enh_act_vector and expr_vector has sufficient variability
                if np.std(enh_act_vector) > 0.05 and np.std(expr_vector) > 0.05:
                    rho, pval = spearmanr(enh_act_vector, expr_vector)

                    if str(rho) != 'nan' and rho > 0:
                        gene_correlations[pval] = [rho, expr_row["Transcript"]]

            # find best correlating transcript for this gene 
            if len(gene_correlations.keys()) > 0:
                min_pval = min(gene_correlations.keys())
                if min_pval < 0.15:
                    correlations += gene_name + "/" + gene_correlations[min_pval][1].split("_")[0] + "/" + str(min_pval) + ";"

        else:
            pass
    if len(correlations) == 0:
        return "."
    else:
        return correlations.rstrip(";")


def check_signal_gene_expression_correlation_enhancer(rare_enriched_enhancer_snps_genes_collected,
                                                      ENHANCER_ACTIVITY, GENE_EXPRESSION):
    h3k27ac_cov = pd.read_csv(ENHANCER_ACTIVITY, sep="\t")
    h3k27ac_cov = h3k27ac_cov.rename(columns={"chr": "CHROM",
                                              "start": "enh_start",
                                              "end": "enh_end"})
    # merge SNPs with H3K27ac coverage on enhancers
    rare_enriched_enhancer_snps_genes_collected_coverage = pd.merge(
        rare_enriched_enhancer_snps_genes_collected,
        h3k27ac_cov, on=["CHROM", "enh_start", "enh_end"], how="left")

    samples = h3k27ac_cov.columns.drop(["CHROM", "enh_start", "enh_end"])
    rare_enriched_enhancer_snps_genes_collected_coverage[
        "H3K27ac-expression correlation p-values"] = rare_enriched_enhancer_snps_genes_collected_coverage.apply(
        calculate_correlation_enh_gene, args=(samples, GENE_EXPRESSION), axis=1)
    rare_enriched_enhancer_snps_genes_collected_corelations = rare_enriched_enhancer_snps_genes_collected_coverage.drop(
        labels=samples, axis=1)
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
            -O {OUTPUT}/input_vcf_sample_level.csv'
    log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                        universal_newlines=True).splitlines()
    all_snps = pd.read_csv(OUTPUT+'/input_vcf_sample_level.csv', sep='\t')
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

        #each genotype group should either be empty or have at least 3 samples
        df = df[((df['Num homref'] == 0) | (df['Num homref'] > 2)) &
                ((df['Num het'] == 0) | (df['Num het'] > 2)) &
                ((df['Num homalt'] == 0) | (df['Num homalt'] > 2))]

        #at least two genotype groups should be represented
        df = df[(df[['Num homref', 'Num het', 'Num homalt']] == 0).sum(axis=1) < 2].copy()
        dataframes[name] = df

    enhancer_snps, promoter_snps = dataframes['enhancer_snps'], dataframes['promoter_snps']

    return promoter_snps, enhancer_snps
        
# calculate correlation between genotype and gene expression, returns string with results of correlation tests
def calculate_correlation_snps_gene(row, counts):
    correlations = ""
    
    variant_columns = [col for col in row.index if '.var' in col]
    homref_indices = [i for i,val in enumerate(row[variant_columns]) if val != './.' and int(float(val)) == 0]
    het_indices = [i for i,val in enumerate(row[variant_columns]) if val != './.' and int(float(val)) == 1]
    homalt_indices = [i for i,val in enumerate(row[variant_columns]) if val != './.' and int(float(val)) == 2]
    
    
    homref_patients = set([variant_columns[i].split('.')[0] for i in homref_indices])
    het_patients = set([variant_columns[i].split('.')[0] for i in het_indices])
    homalt_patients = set([variant_columns[i].split('.')[0] for i in homalt_indices])
    genotype_vector = [0]*len(homref_patients) + [1]*len(het_patients) + [2]*len(homalt_patients)
    
    transcripts = set([el for el in row["transcript_to_check"].split(';')])
    transcripts = ['_'.join(transcript.split('/')) for transcript in transcripts]
    
    #iterate over all target transcripts assigned to this variant
    if transcripts == ['.']:
        return "."
    for transcript in transcripts:
        transcript_name = transcript.replace('/', '_')
        transcript_counts = counts[counts['Transcript'] == transcript_name]
        counts_homref = transcript_counts[list(homref_patients)].values[0]
        counts_het = transcript_counts[list(het_patients)].values[0]
        counts_homalt = transcript_counts[list(homalt_patients)].values[0]

        expression_vector = list(counts_homref) + list(counts_het) + list(counts_homalt) 
        if np.std(expression_vector) > 0.05: 
            rho, pval = spearmanr(genotype_vector, expression_vector)
            if str(rho) != 'nan':
                correlations += transcript + "/pval=" + "%.5f" % pval + "/R=" + "%.5f" % rho + ";"
                        
        else:
            pass
    if len(correlations) == 0:
        return "."
    else:
        return correlations.rstrip(";")

# choose the best candidate target gene transcript based on the lowest pvalue
def find_best_candidate_target(putative_targets):
    if putative_targets != ".":
        putative_targets_list = putative_targets.split(';')
        pvalues = [float(target.split('=')[1].split('/')[0]) for target in putative_targets_list]
        min_pval = min(pvalues)
        best_candidate = putative_targets_list[pvalues.index(min_pval)]
        return pd.Series([min_pval, best_candidate], index=["best_expr_corr_pval", "best_corr_gene"]) 
    else:
        return pd.Series([np.nan, np.nan], index=["best_expr_corr_pval", "best_corr_gene"])
        
    
# get transcripts for each gene, returns string of transcripts separated by ';'
def get_gene_transcript(row, counts):
    if row['Gene'] != '.':
        gene_names = [x.split('/')[1] for x in row['Gene'].split(',')]

        transcripts = ';'.join(counts[counts['Gene'].isin(gene_names)]['Transcript'].values)
        if transcripts == '':
            print('Transcripts not found for gene ', gene_names)
            return '.'
        return transcripts
    else:
        return '.'

def check_gene_snps_correlation(GENE_EXPRESSION, promoter_snps, enhancer_snps):
    # Genes/transcripts normalized counts
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    counts['Gene'] = counts['Transcript'].apply(lambda x: x.split('_')[1])
    
    #assign transcripts to snps related to enhancers and promoters
    enhancer_snps['transcript_to_check'] = enhancer_snps['H3K27ac-expression correlation p-values'].apply(lambda x: x.split('/')[1]+'_'+x.split('/')[0] if x != '.' else '.')
    promoter_snps['transcript_to_check'] = promoter_snps.apply(get_gene_transcript, args=(counts,), axis=1)

    #calculate correlation between genotype and gene expression
    enhancer_snps["gene_expr_correlations"] = enhancer_snps.apply(calculate_correlation_snps_gene, args=(counts,), axis=1)
    promoter_snps["gene_expr_correlations"] = promoter_snps.apply(calculate_correlation_snps_gene, args=(counts,), axis=1)
    enhancer_snps[["best_expr_corr_pval", "best_corr_gene"]] = enhancer_snps["gene_expr_correlations"].apply(find_best_candidate_target)
    promoter_snps[["best_expr_corr_pval", "best_corr_gene"]] = promoter_snps["gene_expr_correlations"].apply(find_best_candidate_target)
    return promoter_snps, enhancer_snps

def visualize_results(promoter_snps, enhancer_snps, GENE_EXPRESSION, OUTPUT, threshold = 0.1):
    promoter_snps = promoter_snps[promoter_snps['best_expr_corr_pval'] < threshold]
    enhancer_snps = enhancer_snps[enhancer_snps['best_expr_corr_pval'] < threshold]
    promoter_snps = promoter_snps.sort_values(by = 'best_expr_corr_pval')
    enhancer_snps = enhancer_snps.sort_values(by = 'best_expr_corr_pval')
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    counts['Gene'] = counts['Transcript'].apply(lambda x: x.split('_')[1])
    
    #save plots to one pdf file
    with PdfPages(OUTPUT+'/plots.pdf') as pdf:  
        for index, row in enhancer_snps.iterrows():
            samples = [col.split('.')[0] for col in enhancer_snps.columns if '.var' in col]
            samples = [el for el in samples if row[el+'.var']!='./.']
            vars = row[[sample + '.var' for sample in samples]]
            transcript  = row['best_corr_gene'].split('/')[0]
            acts = counts[counts["Transcript"]==transcript][samples]
            plt.scatter(vars, acts, s=5, marker='*')
            plt.xlabel('Genotype')
            plt.ylabel('Expression for transcript '+transcript)
            genotype = str(row['Num homref'])+'/'+str(row['Num het'])+'/'+str(row['Num homalt'])
            plt.title('SNP in enhancer: POS '+str(row['POS'])+' '+row['CHROM']+' '+row['REF']+'/'+row['ALT']+' pvalue '+ str(row['best_expr_corr_pval'])+' '+genotype)
            pdf.savefig() 
            plt.close()

        for index, row in promoter_snps.iterrows():
            samples = [col.split('.')[0] for col in promoter_snps.columns if '.var' in col]
            samples = [el for el in samples if row[el+'.var']!='./.']
            vars = row[[sample + '.var' for sample in samples]]
            transcript  = row['best_corr_gene'].split('/')[0]
            acts = counts[counts["Transcript"]==transcript][samples]
            plt.scatter(vars, acts, s=5, marker='*')
            plt.xlabel('Genotype')
            plt.ylabel('Expression for transcript '+transcript)
            genotype = str(row['Num homref'])+'/'+str(row['Num het'])+'/'+str(row['Num homalt'])
            plt.title('SNP in promoter: POS '+str(row['POS'])+' '+row['CHROM']+' '+row['REF']+'/'+row['ALT']+' pvalue '+ str(row['best_expr_corr_pval'])+' '+genotype)
            pdf.savefig()  
            plt.close()
        
        #saving results to csv
        enhancer_snps.to_csv(OUTPUT+'/final_regulatory_snps_enhancer.csv', sep='\t')
        promoter_snps.to_csv(OUTPUT+'/final_regulatory_snps_promoter.csv', sep='\t')


