import pandas as pd
import statistics
import gzip
import tempfile
import psutil
import os

def calc_roh(roh, fai, output):
    dffai = pd.read_table(fai, sep='\t', header=None)
    dffai = dffai[dffai[1] > 10000000]  # only include chromosomes gr 1Mbp
    chroms = dffai[0].values

    glength = dffai[1].sum()

    dfroh = pd.read_table(roh, sep='\t', header=0)
    dfroh = dfroh[dfroh['[6]Length (bp)'] > 500000]  # only look at 500kbp ROH and longer
    dfroh = dfroh[dfroh['[3]Chromosome'].isin(chroms)]

    dffroh = dfroh.groupby(
        ['[2]Sample']
    ).agg(
        {
            '[6]Length (bp)': "sum"
        }
    ).div(glength)
    dffroh = dffroh.reset_index()

    if dffroh.empty:
        print("No individuals found with ROH greater than 500kbp.")
        with open(output, 'w'):
            pass
        with open(output.replace(".froh", "_top.froh"), "w"):
            pass
        return False

    dffroh.to_csv(output, sep='\t', index=False, header=None)
    print(dffroh)

    if len(dffroh) < 8:
        combined_df = dffroh
    else:
        # get the largest 2 ROH
        df_lg = dffroh.nlargest(2, '[6]Length (bp)')
        # get 8 random rows so we have 10 individuals total
        random_rows = dffroh.sample(n=8)
        combined_df = pd.concat([df_lg, random_rows], axis=0)

    combined_df.to_csv(output.replace(".froh", "_top.froh"), sep='\t', index=False, header=False, columns=[combined_df.columns[0]])
    return True

def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss

def vcf2dosage(vcf_file, output_file):
    """
    Converts a VCF file to genotype dosage format.

    Args:
        vcf_file (str): Path to the input VCF file.
        output_file (str): Path to the output genotype dosage file.
    """
    snp_list = []
    temp_files = []

    with gzip.open(vcf_file, 'rb') as vcf_f:
        for line in vcf_f:
            line = line.decode('utf-8').strip()

            if line.startswith('#CHROM'):
                samples = line.split('\t')[9:]
                for _ in samples:
                    temp_files.append(tempfile.TemporaryFile(mode='w+'))
            elif not line.startswith('#'):
                fields = line.split('\t')

                # Skip non-biallelic SNPs
                if ',' in fields[4]:
                    continue

                chrom = fields[0]
                pos = fields[1]
                ref_allele = fields[3]
                alt_allele = fields[4]
                genotypes = fields[9:]
                snp_id = f'{chrom}_{pos}'
                snp_list.append(snp_id)

                # Convert genotype fields to genotype dosage values
                for i, gt in enumerate(genotypes):
                    gt = gt.split(':')[0]
                    if gt == '0/0':
                        dosage = '0'
                    elif gt == '0/1' or gt == '1/0':
                        dosage = '1'
                    elif gt == '1/1':
                        dosage = '2'
                    else:
                        dosage = 'NA'
                    temp_files[i].write(dosage + '\t')

    # Write genotype dosage file
    with gzip.open(output_file, 'wb') as dosage_f:
        # Write header line
        dosage_header = 'Sample_ID\t' + '\t'.join(snp_list) + '\n'
        dosage_f.write(dosage_header.encode('utf-8'))

        # Write genotype dosage values
        for i, temp_file in enumerate(temp_files):
            temp_file.seek(0)
            dosage_row = samples[i] + '\t' + temp_file.read().strip() + '\n'
            dosage_f.write(dosage_row.encode('utf-8'))
            temp_file.close()


# def vcf2dosage(vcf_file, output_file):
#     """
#     Converts a VCF file to genotype dosage format.

#     Args:
#         vcf_file (str): Path to the input VCF file.
#         output_file (str): Path to the output genotype dosage file.
#     """
#     snp_list = []
#     sample_dosages = {}
    
#     with gzip.open(vcf_file, 'rb') as vcf_f:
#         for line in vcf_f:
#             line = line.decode('utf-8').strip()

#             if line.startswith('#CHROM'):
#                 samples = line.split('\t')[9:]
#                 for sample in samples:
#                     sample_dosages[sample] = []
#             elif not line.startswith('#'):
#                 fields = line.split('\t')

#                 # Skip non-biallelic SNPs
#                 if ',' in fields[4]:
#                     continue

#                 chrom = fields[0]
#                 pos = fields[1]
#                 ref_allele = fields[3]
#                 alt_allele = fields[4]
#                 genotypes = fields[9:]
#                 snp_id = f'{chrom}_{pos}'
#                 snp_list.append(snp_id)

#                 # Convert genotype fields to genotype dosage values
#                 for i, gt in enumerate(genotypes):
#                     gt = gt.split(':')[0]
#                     if gt == '0/0':
#                         dosage = '0'
#                     elif gt == '0/1' or gt == '1/0':
#                         dosage = '1'
#                     elif gt == '1/1':
#                         dosage = '2'
#                     else:
#                         dosage = 'NA'
#                     sample_dosages[samples[i]].append(dosage)

#     # Write genotype dosage file
#     with gzip.open(output_file, 'wb') as dosage_f:
#         # Write header line
#         dosage_header = 'Sample_ID\t' + '\t'.join(snp_list) + '\n'
#         dosage_f.write(dosage_header.encode('utf-8'))

#         # Write genotype dosage values
#         for sample, dosages in sample_dosages.items():
#             dosage_row = sample + '\t' + '\t'.join(dosages) + '\n'
#             dosage_f.write(dosage_row.encode('utf-8'))


# def vcf2dosage(vcf_file, output_file):
#     """
#     Converts a VCF file to genotype dosage format.

#     Args:
#         vcf_file (str): Path to the input VCF file.
#         output_file (str): Path to the output genotype dosage file.
#     """
#     with gzip.open(vcf_file, 'rb') as f:
#         # Read and decode header lines
#         header = ''
#         for line in f:
#             line = line.decode('utf-8').strip()
#             if line.startswith('#CHROM'):
#                 header = line
#                 break
#         # Parse sample IDs
#         samples = header.split('\t')[9:]
#         num_samples = len(samples)
#         # Initialize genotype dosage matrix
#         geno_dosage = [[] for i in range(num_samples)]
#         # Parse genotypes
#         for line in f:
#             fields = line.decode('utf-8').strip().split('\t')
#             # Skip non-biallelic SNPs
#             if ',' in fields[4]:
#                 continue
#             # Parse CHROM and POS fields
#             chrom = fields[0]
#             pos = fields[1]
#             # Parse reference and alternate alleles
#             ref_allele = fields[3]
#             alt_allele = fields[4]
#             # Parse genotype fields
#             genotypes = fields[9:]
#             # Convert genotype fields to genotype dosage values
#             for i in range(num_samples):
#                 gt = genotypes[i].split(':')[0]
#                 if gt == '0/0':
#                     dosage = 0
#                 elif gt == '0/1' or gt == '1/0':
#                     dosage = 1
#                 elif gt == '1/1':
#                     dosage = 2
#                 else:
#                     dosage = 'NA'
#                 geno_dosage[i].append((chrom, pos, dosage))
#         # Write genotype dosage file
#         with gzip.open(output_file, 'wb') as f:
#             # Write header line
#             f.write(('Sample_ID\t' + '\t'.join([chrom + '_' + pos for chrom, pos, dosage in geno_dosage[0]]) + '\n').encode('utf-8'))
#             # Write genotype dosage values
#             for i in range(num_samples):
#                 f.write((samples[i] + '\t' + '\t'.join([str(dosage) for chrom, pos, dosage in geno_dosage[i]]) + '\n').encode('utf-8'))