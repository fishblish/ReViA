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
        result_files[r] = f'{OUTPUT}/{r}_SNPs.vcf' #path to files where results will be saved
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
        variants_file=variants_vcf_file_dict[r]
        out_file = variants_file.split('.')[-2]
        
        command = f"perl {ANNOVAR}table_annovar.pl {variants_file} {ANNOVAR}humandb/ -buildver hg38 -remove -protocol gnomad_genome -operation f -nastring . -vcfinput -out {out_file} -polish" 
        try:
            log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
            annovar_logs.append(log)
            result_vcf_files_dict[r]=f'{out_file}.hg38_multianno.vcf'
        except subprocess.CalledProcessError as e:
            print(f"Command '{command}' failed with error code {e.returncode}:\n{e.output}")
        print(f"Done: frequencies annotated for variants in {r} regions")
    return result_vcf_files_dict


def select_snps_by_freq(annotated_variants_file_dict, GATK, population, treating_gaps, target='r', cutoff=0.01):
    result_files_dict = {}
    gaps_dict = {'r': '=0.0', 'c': '=100.0'}
    ineq_sign = {'r': '<', 'c': '>'}
    target_full_word = {'r': 'rare', 'c': 'common'}

    annotations = ['gnomAD_genome_' + p for p in population]

    print('There will be selected variants which are', target_full_word[target], 'among these populations', population,
          'for frequency cutoff equal to', cutoff, 'and missing variants in annovar database treated as', target_full_word[treating_gaps])

    for r in ["promoter", "enhancer"]:
        result_file_name = annotated_variants_file_dict[r].split('.')
        result_file_name.insert(-1, 'nomissing')
        result_files_dict[r]='.'.join(result_file_name)
        with open(result_files_dict[r], 'w') as o:
            for line in open(annotated_variants_file_dict[r]).readlines():
                for el in annotations:
                    if el + '=.' in line:
                        line = line.replace(el + '=.', el + gaps_dict[treating_gaps])
                o.write(line)

    select_rare_logs = []

    select_condition = '"'

    for a in annotations:
        select_condition = select_condition + a + ' ' + ineq_sign[target] + ' ' + str(cutoff) + ' && '
    select_condition = select_condition[0:-3]
    select_condition = select_condition + '"'
    result_filtered_files_dict = {}
    below_above={'r': 'below', 'c': 'above'}
    for r in ["promoter", "enhancer"]:

        command = f"{GATK} VariantsToTable  -V {result_files_dict[r]} " \
            "-F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN " \
            f"-F gnomAD_genome_ALL " \
            f"-O {result_files_dict[r].replace('.vcf', '.csv')}"
        log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
        
        table= pd.read_csv(result_files_dict[r].replace('.vcf', '.csv'), sep='\t')
        print('UWAGA1!!!!!!!!!', table[table['POS'] == 172400207])
        
        result_filtered_files_names = result_files_dict[r].split('.')
        result_filtered_files_names.insert(-1, 'gnomad_'+below_above[target]+'_'+str(cutoff))
        result_filtered_files_dict[r]='.'.join(result_filtered_files_names)
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

def calc_binom_pval(row, p_col='gnomAD_genome_NFE'):
    x = row['AC']
    n = row['AN']
    p = float(row[p_col])
    if p > 1.0:
        print(row)

    return binom_test(x, n, p, alternative='greater')


def select_enriched_snps(filtered_variants_files_dict, GATK, OUTPUT, reference_population="NFE", bh_alpha=0.01):
    # vcf file to csv table
    reference_population_cmd = ""
    if reference_population != "ALL":
        reference_population_cmd = f"-F gnomAD_genome_{reference_population} "
    totable_logs = []
    for r in ["promoter", "enhancer"]:
        command = f"{GATK} VariantsToTable  -V {filtered_variants_files_dict[r]} " \
                  "-F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN " \
                  f"-F gnomAD_genome_ALL {reference_population_cmd}" \
                  f"-O {filtered_variants_files_dict[r].replace('.vcf', '.csv')}"
        log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                      universal_newlines=True).splitlines()
        totable_logs.append(log)
        print("Done: vcf file to table")
    rare_promoter_snps = pd.read_csv(filtered_variants_files_dict['promoter'].replace('.vcf', '.csv'), sep='\t')
    rare_enhancer_snps = pd.read_csv(filtered_variants_files_dict['enhancer'].replace('.vcf', '.csv'), sep='\t')
    print('UWAGA!!!!!!!!!', rare_enhancer_snps[rare_enhancer_snps['POS'] == 172400207])

    # execute one-sided binom test:
    # no. successes = no. ALT alleles in cohort
    # no. trials = no. all identified alleles
    # prob. of success = annotated popultaion frequency.
    # alt. hypetesis = observed freq is greated than expected

    for df in [rare_enhancer_snps, rare_promoter_snps]:
        df['binom_pval'] = df.apply(calc_binom_pval, axis=1)

        # apply correction for multiple hypothesis testing with the Benjamini-Hochberg procedure, use FDR = bh_alpha
        multipletests_correction = multipletests(df['binom_pval'], alpha=bh_alpha,
                                                 method='fdr_bh', is_sorted=False, returnsorted=False)
        df['B-H_reject_H0'] = multipletests_correction[0]
        df['corrected_binom_pval'] = multipletests_correction[1]
    rare_enriched_promoter_snps = rare_promoter_snps[rare_promoter_snps["B-H_reject_H0"]]
    rare_enriched_enhancer_snps = rare_enhancer_snps[rare_enhancer_snps["B-H_reject_H0"]]
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


def assign_genes_intronic_enhancer_snps(rare_enriched_enhancer_snps_df, ENHANCER_REGIONS):
    enh_genes = pd.read_csv(ENHANCER_REGIONS, sep='\t', names=['chr', 'start', 'end', 'Gene'])

    # reformat to have one gene ID in cell
    enh_gene = pd.DataFrame()
    for i, row in enh_genes.iterrows():
        genes = row['Gene'].split(',')
        if len(genes) == 1:
            enh_gene = pd.concat([enh_gene, pd.DataFrame([row])], ignore_index=True)
        else:
            for gene in genes:
                new_row = row
                new_row['Gene'] = gene
                enh_gene = pd.concat([enh_gene, pd.DataFrame([new_row])], ignore_index=True)
                enh_gene = enh_gene.reindex(enh_genes.columns, axis=1)
                enh_gene["start"] = enh_gene["start"].astype(int)
                enh_gene["end"] = enh_gene["end"].astype(int)

    # Intersect information about enhancers with SNPs to assign gene names to SNPs.
    # prepare bedtool objects
    rare_enriched_enhancer_snps_df["POS-1"] = rare_enriched_enhancer_snps_df["POS"] - 1
    rare_enriched_enhancer_snps_bedtool = pbt.BedTool.from_dataframe(
        rare_enriched_enhancer_snps_df[["CHROM", "POS-1", "POS", "REF", "ALT", "AC", "AF", "AN",
                                           "gnomAD_genome_ALL", "gnomAD_genome_NFE", "binom_pval",
                                           "B-H_reject_H0", "corrected_binom_pval"]])
    enh_genes_bedtool = pbt.BedTool.from_dataframe(enh_gene)

    # intersect 
    rare_enriched_enhancer_snps_intersection = rare_enriched_enhancer_snps_bedtool.intersect(
        enh_genes_bedtool, wa=True, wb=True, loj=True)
    # reformat intersection to dataframe, keep columns with enhancer coordinates - they will be usefull in the next step
    #dropping POS-1 column
    rare_enriched_enhancer_snps_gene = rare_enriched_enhancer_snps_intersection.to_dataframe(
        usecols=[0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16],
        names=["CHROM", "POS", "REF", "ALT", "AC", "AF", "AN",
               "gnomAD_genome_ALL", "gnomAD_genome_NFE", "binom_pval",
               "B-H_reject_H0", "corrected_binom_pval", "enh_start", "enh_end", "Gene"])
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


def reformat_target_genes_enh(rare_enriched_enhancer_snps_gene_closest_contacts):
    rare_enriched_enhancer_snps_genes_collected = pd.DataFrame()

    for name, group in rare_enriched_enhancer_snps_gene_closest_contacts.groupby(["CHROM", "POS", "REF", "ALT"]):
        containing_genes = [gene + "(containing)" for gene in group["Gene"].unique() if gene != "."]
        closest_genes = [gene + "(closest)" for gene in group["closest gene"].unique() if gene != "."]
        contacting_genes = [gene + "(contacting)" for gene in group["contacting gene"].unique() if gene != "."]

        all_genes = containing_genes + closest_genes + contacting_genes

        group["Gene"] = ";".join(all_genes)
        row = group[['CHROM', 'POS', 'REF', 'ALT',
                     'AC', 'AF', 'AN', 'gnomAD_genome_ALL',
                     'gnomAD_genome_NFE', 'binom_pval', 'B-H_reject_H0',
                     'corrected_binom_pval', 'enh_start', 'enh_end', 'Gene', 
                     'genomic element']].drop_duplicates()

        rare_enriched_enhancer_snps_genes_collected = pd.concat(
            [rare_enriched_enhancer_snps_genes_collected, pd.DataFrame(row)], ignore_index=True)
    return rare_enriched_enhancer_snps_genes_collected


def calculate_correlation_enh_gene(row, sample_names, GENE_EXPRESSION):
    # Calculate correlation between enhancer activity and gene expression
    # Genes/transcripts normalized counts
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
                if np.std(enh_act_vector) > 0.05 and np.std(expr_vector) > 0.1:
                    rho, pval = spearmanr(enh_act_vector, expr_vector)

                    if str(rho) != 'nan' and rho > 0:
                        gene_correlations[pval] = [rho, expr_row["Transcript"]]

                else:
                    print('Too low variance in gene exprression (', np.std(expr_vector), ')  for gene ', gene_name, 'transcript',
                          expr_row["Transcript"])

            # find best correlating transcript
            if len(gene_correlations.keys()) > 0:
                min_pval = min(gene_correlations.keys())
                if min_pval < 0.15:
                    correlations += gene_name + "/" + gene_correlations[min_pval][1].split("_")[
                        0] + "/" + "%.5f" % min_pval + ";"

        else:
            pass
    if len(correlations) == 0:
        return "."
    else:
        return correlations.rstrip(";")


def find_best_candidate_target1(putative_targets):
    if putative_targets != ".":
        putative_targets_list = putative_targets.split(';')
        pvalues = [float(target.split('/')[2]) for target in putative_targets_list]
        min_pval = min(pvalues)
        best_candidate = putative_targets_list[pvalues.index(min_pval)]
        print(best_candidate)
        return best_candidate
    else:
        return '.'


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
    rare_enriched_enhancer_snps_genes_collected_corelations["Putative target with highest correlation"] = \
    rare_enriched_enhancer_snps_genes_collected_corelations["H3K27ac-expression correlation p-values"].apply(
        find_best_candidate_target1)
    return rare_enriched_enhancer_snps_genes_collected_corelations


def get_gene_names(genes_string):
    # promoters will have ENSG00000136026/CKAP4, comma separated
    # enhancers will have ";"-separated lists with the following format: ENSG00000171735/CAMTA1(containing)
    if genes_string:
        if "(" not in genes_string:
            return [el.split('/')[1] for el in genes_string.split(',')]
        else:
            return [el.split('/')[1].split('(')[0] for el in genes_string.split(';')]
    else:
        return ""


def check_expression_in_brain(genes, GTEX):
    gtex = pd.read_csv(GTEX, sep='\t', skiprows=[0, 1])
    brain_columns = [col for col in list(gtex.columns) if "Brain" in col]
    gene_names_list = get_gene_names(genes)
    expression_list = []
    for gene in gene_names_list:
        if gene != "" and gene != ".":
            try:
                mean_median_tpm = sum(gtex[gtex['Description'] == gene.strip()][brain_columns].values[0]) / float(
                    len(brain_columns))
                expression_list.append(gene + ':' + "%.2f" % mean_median_tpm)
            except:
                print("no gtex brain data for:", gene)

    return ",".join(expression_list)


def assign_median_tpm_gtex(rare_enriched_promoter_snps_gene,
                           rare_enriched_enhancer_snps_genes_collected_corelations, GTEX):
    rare_enriched_promoter_snps_gene[
        "Median TPM in brain tissues in GTEx"] = rare_enriched_promoter_snps_gene.Gene.apply(
        check_expression_in_brain, args=(GTEX,))
    rare_enriched_enhancer_snps_genes_collected_corelations[
        "Median TPM in brain tissues in GTEx"] = rare_enriched_enhancer_snps_genes_collected_corelations.Gene.apply(
        check_expression_in_brain, args=(GTEX,))
    return rare_enriched_promoter_snps_gene, rare_enriched_enhancer_snps_genes_collected_corelations

def import_vcf_sample_level(INPUT_VCF, OUTPUT, GATK, promoter_snps, enhancer_snps):
    command = f'{GATK} VariantsToTable -V {INPUT_VCF}  \
        -F CHROM -F POS -F REF -F ALT -F AC -F AF -F AN -GF AD -GF GT \
            -O {OUTPUT}\input_vcf_sample_level.csv'
    log = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT,
                                        universal_newlines=True).splitlines()
    all_snps = pd.read_csv(OUTPUT+'input_vcf_sample_level.csv', sep='\t')
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

        # Assign the modified dataframe back to the original variable
        dataframes[name] = df

    # Update the original dataframes
    enhancer_snps, promoter_snps = dataframes['enhancer_snps'], dataframes['promoter_snps']

    return promoter_snps, enhancer_snps
        
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
    #print('transcripts: ',transcripts)
    #iterate over all target transcripts assigned to this variant
    if transcripts == ['.']:
        return "."
    for transcript in transcripts:
        transcript_name = transcript.replace('/', '_')
        transcript_counts = counts[counts['Transcript'] == transcript_name]
        #print('transcript: ',transcript)
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
    
def find_best_candidate_target(putative_targets):
    if putative_targets != ".":
        putative_targets_list = putative_targets.split(';')
        pvalues = [float(target.split('=')[1].split('/')[0]) for target in putative_targets_list]
        min_pval = min(pvalues)
        best_candidate = putative_targets_list[pvalues.index(min_pval)]
        return pd.Series([min_pval, best_candidate], index=["best_expr_corr_pval", "best_corr_gene"]) 
    else:
        return pd.Series([np.nan, np.nan], index=["best_expr_corr_pval", "best_corr_gene"])
        
    

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
    enhancer_snps['transcript_to_check'] = enhancer_snps['H3K27ac-expression correlation p-values'].apply(lambda x: x.split('/')[1]+'_'+x.split('/')[0] if x != '.' else '.')

    #add transcripts to promoter
    print('Assigning transcripts to promoters')
    promoter_snps['transcript_to_check'] = promoter_snps.apply(get_gene_transcript, args=(counts,), axis=1)

    #promoter_snps['transcript_to_check'] = promoter_snps['Gene']
    enhancer_snps["gene_expr_correlations"] = enhancer_snps.apply(calculate_correlation_snps_gene, args=(counts,), axis=1)
    promoter_snps["gene_expr_correlations"] = promoter_snps.apply(calculate_correlation_snps_gene, args=(counts,), axis=1)
    enhancer_snps[["best_expr_corr_pval", "best_corr_gene"]] = enhancer_snps["gene_expr_correlations"].apply(find_best_candidate_target)
    promoter_snps[["best_expr_corr_pval", "best_corr_gene"]] = promoter_snps["gene_expr_correlations"].apply(find_best_candidate_target)
    print(promoter_snps.sort_values(by = "best_expr_corr_pval")[:10])
    print(enhancer_snps.sort_values(by = "best_expr_corr_pval")[:10])
    return promoter_snps, enhancer_snps

def visualize_results(promoter_snps, enhancer_snps, GENE_EXPRESSION, OUTPUT, threshold = 0.1):
    promoter_snps = promoter_snps[promoter_snps['best_expr_corr_pval'] < threshold]
    enhancer_snps = enhancer_snps[enhancer_snps['best_expr_corr_pval'] < threshold]
    promoter_snps = promoter_snps.sort_values(by = 'best_expr_corr_pval')
    enhancer_snps = enhancer_snps.sort_values(by = 'best_expr_corr_pval')
    counts = pd.read_csv(GENE_EXPRESSION, sep='\t')
    counts['Gene'] = counts['Transcript'].apply(lambda x: x.split('_')[1])
    

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
            pdf.savefig()  # saves the current figure into a pdf page
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
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
        enhancer_snps.to_csv(OUTPUT+'/regulatory_snps_enhancer.csv', sep='/t')
        promoter_snps.to_csv(OUTPUT+'/regulatory_snps_promoter.csv', sep='/t')

#motywy

# Annotate snps with predicted transcript factor binding site
# Identifing SNPs which can destroy or create a TF binding site

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
        output_bed_path = "%s/%s_rare_enriched_SNPs.bed" % (OUTPUT, r)
        snps_bed.to_csv(output_bed_path, sep="\t", index=False, header=False)
        snps_bed_files.append(output_bed_path)
    return snps_bed_files


def search_motifs(path_to_r_packages='.libPaths("/home/julia/miniconda3/lib/R/library/")'):
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)

    robjects.r(path_to_r_packages)

    robjects.r('''
        library('motifbreakR')
        library('BSgenome.Hsapiens.UCSC.hg38')
        library('MotifDb')
        motifs <- query(MotifDb, andStrings=c("hocomocov11", "hsapiens"))
        print(paste('Number of motifs is', length(motifs)))

    ''')


def score_motifs(snps_bed_files, path_to_r_packages='.libPaths("/home/julia/miniconda3/lib/R/library/")'):
    import rpy2.robjects as robjects
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)

    robjects.r(path_to_r_packages)

    robjects.r('''
        library('motifbreakR')
        library('BSgenome.Hsapiens.UCSC.hg38')
        library('MotifDb')
        
        score_snps <- function(snps_file, out_file) {
            #read SNPs from input bed file
            snps.mb.frombed <- snps.from.file(file = snps_file, search.genome = BSgenome.Hsapiens.UCSC.hg38, format = "bed")
        
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
        } 
    ''')

    score_snps_py = robjects.globalenv['score_snps']
    for snp_bed in snps_bed_files:
        snp_scores_csv = snp_bed.replace(".bed", "_motifbreakR-scores.csv")
        print("Calculate scores for input: %s, save output to: %s" % (snp_bed, snp_scores_csv))
        snp_scores = score_snps_py(snp_bed, snp_scores_csv)
        print("Done")
    return snp_scores


def find_best_matching_motif(group):
    # find motif with highest pct score (either for REF or ALT)
    best_pctRef_score = max(group["pctRef"])
    best_pctAlt_score = max(group["pctAlt"])
    if best_pctRef_score > best_pctAlt_score:
        best_pct_score_motif = group[group["pctRef"] == best_pctRef_score]["providerId"].values[0]
    else:
        best_pct_score_motif = group[group["pctAlt"] == best_pctAlt_score]["providerId"].values[0]

    # find motif with highest abs(diff) between alleles
    best_alleleDiff = max(group["alleleDiff"].abs())
    best_alleleDiff_motif = group[group["alleleDiff"].abs() == best_alleleDiff]["providerId"].values[0]

    return best_pct_score_motif + ":" + "%.2f" % max(best_pctRef_score,
                                                     best_pctAlt_score), best_alleleDiff_motif + ":" + "%.2f" % best_alleleDiff


def select_motif_results(OUTPUT, rare_enriched_promoter_snps, rare_enriched_enhancer_snps):
    promoter_SNPs_motifbreakr = pd.read_csv("%s/promoter_rare_enriched_SNPs_motifbreakR-scores.csv" % OUTPUT, sep="\t")
    enhancer_SNPs_motifbreakr = pd.read_csv("%s/enhancer_rare_enriched_SNPs_motifbreakR-scores.csv" % OUTPUT, sep="\t")
    # select records with "strong" effect
    promoter_SNPs_motifbreakr_strong = promoter_SNPs_motifbreakr[promoter_SNPs_motifbreakr["effect"] == "strong"]
    enhancer_SNPs_motifbreakr_strong = enhancer_SNPs_motifbreakr[enhancer_SNPs_motifbreakr["effect"] == "strong"]
    # add columns with info about best matches: motif with highest pct score and alleleDiff

    # information about best motifs for each SNP will be stored in a dict in which SNP_ids will be keys
    best_motifs_dict = {}

    for df in [enhancer_SNPs_motifbreakr_strong, promoter_SNPs_motifbreakr_strong]:
        for snp_id, snp_records in df.groupby("SNP_id"):
            best_match, highest_diff = find_best_matching_motif(snp_records)
            best_motifs_dict[snp_id] = {"best_match": best_match, "highest_diff": highest_diff}

    enhancer_SNPs_motifbreakr_strong = enhancer_SNPs_motifbreakr_strong.copy()
    promoter_SNPs_motifbreakr_strong = promoter_SNPs_motifbreakr_strong.copy()
    # extract information from the dict to fill appropriate columns in enhancer and promoter SNPs dataframes
    enhancer_SNPs_motifbreakr_strong.loc[:, "motif_best_match"] = enhancer_SNPs_motifbreakr_strong.loc[:,
                                                                  "SNP_id"].apply(
        lambda x: best_motifs_dict[x]["best_match"])
    enhancer_SNPs_motifbreakr_strong.loc[:, "motif_highest_diff"] = enhancer_SNPs_motifbreakr_strong.loc[:,
                                                                    "SNP_id"].apply(
        lambda x: best_motifs_dict[x]["highest_diff"])

    promoter_SNPs_motifbreakr_strong.loc[:, "motif_best_match"] = promoter_SNPs_motifbreakr_strong.loc[:,
                                                                  "SNP_id"].apply(
        lambda x: best_motifs_dict[x]["best_match"])
    promoter_SNPs_motifbreakr_strong.loc[:, "motif_highest_diff"] = promoter_SNPs_motifbreakr_strong.loc[:,
                                                                    "SNP_id"].apply(
        lambda x: best_motifs_dict[x]["highest_diff"])

    # extract information about SNP location and best motifs, drop duplicates
    promoter_SNPs_motifbreakr_strong_snps_only = promoter_SNPs_motifbreakr_strong[
        ["seqnames", "start", "REF", "ALT", "motif_best_match", "motif_highest_diff"]].drop_duplicates()
    enhancer_SNPs_motifbreakr_strong_snps_only = enhancer_SNPs_motifbreakr_strong[
        ["seqnames", "start", "REF", "ALT", "motif_best_match", "motif_highest_diff"]].drop_duplicates()

    # change column names to keep the convention used in the whole notebook
    promoter_SNPs_motifbreakr_strong_snps_only = promoter_SNPs_motifbreakr_strong_snps_only.rename(
        columns={"seqnames": "CHROM", "start": "POS"})
    enhancer_SNPs_motifbreakr_strong_snps_only = enhancer_SNPs_motifbreakr_strong_snps_only.rename(
        columns={"seqnames": "CHROM", "start": "POS"})

    print(len(promoter_SNPs_motifbreakr_strong_snps_only), "and", len(enhancer_SNPs_motifbreakr_strong_snps_only),
          "SNPs in promoters and enhancers have predicted strong effect of motif binding (out of %s and %s, respectively)." % (
          str(len(rare_enriched_promoter_snps)), str(len(rare_enriched_enhancer_snps))))

    rare_enriched_promoter_snps_motif = pd.merge(rare_enriched_promoter_snps,
                                                 promoter_SNPs_motifbreakr_strong_snps_only, how="right",
                                                 on=["CHROM", "POS", "REF", "ALT"])
    rare_enriched_enhancer_snps_motif = pd.merge(rare_enriched_enhancer_snps,
                                                 enhancer_SNPs_motifbreakr_strong_snps_only, how="right",
                                                 on=["CHROM", "POS", "REF", "ALT"])
    return rare_enriched_promoter_snps_motif, rare_enriched_enhancer_snps_motif
