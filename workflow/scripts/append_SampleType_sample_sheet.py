import pandas as pd
import shutil
import argparse

# create a single column text file listing sample names to exclude
# e.g. after running QC, identify problem samples, add them to a file 
# then this script will create a new sample sheet that will prompt post processing to run

def append_col(sample_list, excl_list):
    # Load the CSV file
    csv_file = pd.read_csv(sample_list)
    #csv_file=csv_file.drop(['SampleType'], axis = 1)
    print(csv_file)
    # copy the original as backup
    shutil.copy2(sample_list, sample_list+'.bak')

    # Load the tab-delimited file
    tab_file = pd.read_csv(excl_list, sep='\t', header=None, names=['BioSample'])
    # add a column for exclude
    tab_file['SampleType'] = 'exclude'
    # Merge the two dataframes on the 'BioSample' column
    merged_df = pd.merge(csv_file, tab_file, on='BioSample', how='left')
    # Save the updated CSV file
    merged_df.to_csv(sample_list, index=False)

def main() -> None:

    parser = argparse.ArgumentParser(description='Write sample files.')
    parser.add_argument('-s', '--sample_list', dest='samp', required=True, help="Specify path to sample list")
    parser.add_argument('-e', '--excl', dest='excl', required=True, help="Specify path to sample list")
    args = parser.parse_args()

    sample_list = args.samp
    excl_list = args.excl

    append_col(sample_list, excl_list)

main()