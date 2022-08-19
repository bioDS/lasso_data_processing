#!/usr/bin/env python

# Takes the output from predicttargets.py and a list of siRNAs and genes (1 gene and siRNA per line,
# respectively) and condenses the information to a sparse Matrix (mtx format) of siRNAs by genes with raw 
# target relation scores in the cells. Later this matrix can be read into R using the sparse matrix package "Matrix".
#
# Author: Fabian Schmich (fabian.schmich@bsse.ethz.ch)
#

# Libraries
import sys
import os.path
import optparse
import subprocess
from os import remove
import tempfile

DEBUG = True
CPS_DELIM = "\t"
CPS_HAS_HEADER = True
CPS_POS_GENE = 0
CPS_POS_SCORE = 4

def main():
    # Parse Arguments
    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option(
        "--cps",
        action="store",
        type="string",
        dest="cps",
        help="context+ scores folder"
    )
    parser.add_option(
        "--sirnas",
        action="store",
        type="string",
        dest="siRNAs",
        help="siRNAs file"
    )
    parser.add_option(
        "--genes",
        action="store",
        type="string",
        dest="genes",
        help="genes file"
    )
    parser.add_option(
        "--out",
        action="store",
        type="string",
        dest="output",
        help="output file"
    )
    (options, args) = parser.parse_args(sys.argv)

    # Open input
    if options.cps is None or options.genes is None or options.siRNAs is None:
        exit("ERROR: Missing input")
    if not (os.path.isfile(options.genes) or os.path.isfile(options.siRNAs)):
        exit("ERROR: Not a file")
    sif = open(options.siRNAs, "r")
    genf = open(options.genes, "r")

    # Read genes file
    genpos = dict()
    i = 1
    for line in genf:
        gen = line.strip()
        if gen in genpos: 
            print "WARN: gene %s already hashed!" % gen
        else: 
            genpos[gen] = i
            i += 1

    # Read siRNA file
    sipos = dict()
    i = 1
    for line in sif:
        si = line.strip()
        if si in sipos:
            print "WARN: siRNA  %s already hashed" % si
        else:
            print "adding siRNA %s" % si
            sipos[si] = i
            i += 1

    # Write Matrix
    if DEBUG: print("Writing matrix entries")
    out_mattemp = tempfile.NamedTemporaryFile("w", prefix="mat_", delete=False)
    elemcount = 0
    for si in sipos.keys():
        #print(si)
        line = {}
        count = 0
        cpsfile = options.cps + "/" + si + ".txt"
        while os.path.isfile(cpsfile):
            cpsf = open(cpsfile, "r")
            print("opened {}".format(cpsfile))
            if CPS_HAS_HEADER: cpsf.readline()
            for sline in cpsf:
                # print("sline: " + sline)
                lsp = sline.split(CPS_DELIM)
                gene = lsp[CPS_POS_GENE].strip()
                #print(lsp)
                #print(CPS_POS_SCORE)
                score = float(lsp[CPS_POS_SCORE])
                #out_mattemp.write(" ".join([str(sipos[si]), str(genpos[gene]), str(score)]) + "\n")
                line[gene] = line.get(gene, 0) + score;
                elemcount += 1
            cpsf.close()
            cpsfile = options.cps + "/" + si + ".txt" + str(count)
            count = count + 1
        else:
            print("WARN: No scores file found for %s.%d" % (si, count))
        for (gene, score) in line.items():
            out_mattemp.write(" ".join([str(sipos[si]), str(genpos[gene]), str(score)]) + "\n")

    out_mattemp.close()

    # Write matrix header
    if DEBUG: print("Writing out matrix header")
    out_mathead = tempfile.NamedTemporaryFile("w", prefix="header_", delete=False)
    out_mathead.write("%%MatrixMarket matrix coordinate real general\n")
    dims = [len(sipos.values()), len(genpos.keys()), elemcount]
    out_mathead.write(" ".join([str(e) for e in dims])+"\n")
    out_mathead.close()

    # Write sirnas and genes
    if DEBUG: print("Writing out rownames and colnames")
    out_sir = open(options.output + ".rownames", "w")
    out_gen = open(options.output + ".colnames", "w")
    for key, value in sipos.items():
        out_sir.write("\t".join([str(value), key]) + "\n")
    for key, value in genpos.items():
        out_gen.write("\t".join([str(value), key]) + "\n")
    out_sir.close()
    out_gen.close()

    # Join header and matrix
    if DEBUG: print("Writing out final matirx")
    p = subprocess.Popen("cat %s %s > %s" % (out_mathead.name, out_mattemp.name, options.output+".mtx"), shell=True)
    p.wait()
    remove(out_mattemp.name)
    remove(out_mathead.name)


if __name__ == "__main__":
    main()
