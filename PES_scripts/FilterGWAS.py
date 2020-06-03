import csv
import argparse
import sys
from copy import copy

def check_args(args=None):
    """Parse arguments"""
    parser = argparse.ArgumentParser(description="Filter a GWAS by a list of genomic regions.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gwas', help="Tab-separated GWAS summary statistics file containing at least the chromosome and position of the variants. Accepted chromosome codes: 1-22 for autosomal chromosomes, X or 23 for the X chromosome, Y or 24 for the Y chromosome and M, MT or 25 for mitochondrial DNA.", required=True)
    parser.add_argument('-f', '--filter', help="Tab-separated filter list, containing three columns: chromosome, start position and end position. Expects 1-based inclusive coordinates.", required=True)
    parser.add_argument('-c', '--chr', help="Chromosome column number in GWAS file.", required=True)
    parser.add_argument('-b', '--bp', help="Position column number in GWAS file.", required=True)
    parser.add_argument('-o', '--output', help="Output file name for filtered GWAS.", required=True)
    return(parser.parse_args(args))

def filterGwas(gwasFile, filterFile, chrCol, bpCol):
    c = int(chrCol) - 1
    b = int(bpCol) - 1
    filterDict = {}
    for i in range(1, 26):
        filterDict[i] = []
    gwas = []
    with open(filterFile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            chr = line[0]
            try:
                chr = int(chr)
            except ValueError:
                if chr == 'X':
                    chr = 23
                elif chr == 'Y':
                    chr = 24
                elif chr == 'M' or chr == 'MT':
                    chr = 25
                else:
                    continue
            start = int(line[1])
            end = int(line[2]) + 1
            filterDict[chr].append(range(start, end))
    with open(gwasFile, 'r') as f:
        filteredGwasList = []
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for line in reader:
            chr = copy(line[c])  # chr column
            try:
                chr = int(chr)
            except ValueError:
                if chr == 'X':
                    chr = 23
                elif chr == 'Y':
                    chr = 24
                elif chr == 'M' or chr == 'MT':
                    chr = 25
                else:
                    continue
            bp = int(line[b])  # bp column
            if any(bp in i for i in filterDict[chr]):
                filteredGwasList.append(line)
    return(filteredGwasList)

arguments = check_args(sys.argv[1:])
newGwas = filterGwas(arguments.gwas, arguments.filter, arguments.chr, arguments.bp)
newGwas = ['\t'.join([str(i) for i in j]) for j in newGwas]
with open(arguments.output, 'w') as o:
    o.write('\n'.join(newGwas) + '\n')
