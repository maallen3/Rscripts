import pandas
import os
import sys

def main(fn):
    rootname = os.path.basename(fn)
    rootname= rootname.strip(".bed")
    outdir = os.path.dirname(fn)
    outfile = outdir+"/"+rootname+".txt"
    df = pandas.read_csv(fn, sep="\t",usecols=[3], names=["transcript_genename"])
    df['genename'] = df['transcript_genename'].str.split("_").str[-1]
    df = df[["genename"]].drop_duplicates()
    df = df.sort_values(by="genename")
    df.to_csv(outfile, header=False, index=False)


main(sys.argv[1])
