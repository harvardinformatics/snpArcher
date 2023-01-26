import pandas as pd
import statistics

def calc_ruh(roh, fai, output):
    dffai = pd.read_table(fai, sep='\t', header = None)
    glength = dffai[1].sum()

    dfroh = pd.read_table(roh, sep='\t', header = 0)

    dffroh = dfroh.groupby(
        ['[2]Sample']
        ).agg(
            {
                '[6]Length (bp)': "sum" 
            }
        ).div(glength)
    dffroh = dffroh.reset_index()
    dffroh.to_csv(output, sep='\t', index=False, header = None)

    dffroh = dffroh.nlargest(2, '[6]Length (bp)')
    dffroh.to_csv(output.replace(".froh","_top.froh"), sep='\t', index=False, header = False, columns = [dffroh.columns[0]])

