import pandas as pd
import statistics

def calc_roh(roh, fai, output):
    dffai = pd.read_table(fai, sep='\t', header = None)
    dffai = dffai[dffai[1] > 10000000] # only include chromosomes gr 1Mbp
    chroms = dffai[0].values

    glength = dffai[1].sum()

    dfroh = pd.read_table(roh, sep='\t', header = 0)
    dfroh = dfroh[dfroh['[6]Length (bp)'] > 500000] # only look at 500kbp ROH and longer
    dfroh = dfroh[dfroh['[3]Chromosome'].isin(chroms)]

    dffroh = dfroh.groupby(
        ['[2]Sample']
        ).agg(
            {
                '[6]Length (bp)': "sum" 
            }
        ).div(glength)
    dffroh = dffroh.reset_index()
    dffroh.to_csv(output, sep='\t', index=False, header = None)
    print(dffroh)

    if len(dffroh) < 8:
        combined_df = dffroh
    else:
        #get the largest 2 ROH
        df_lg = dffroh.nlargest(2, '[6]Length (bp)')
        #get 8 random rows so we have 10 individuals total
        random_rows = dffroh.sample(n=8)
        combined_df = pd.concat([df_lg, random_rows], axis=0)

    combined_df.to_csv(output.replace(".froh","_top.froh"), sep='\t', index=False, header = False, columns = [combined_df.columns[0]])


