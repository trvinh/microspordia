import sys
import os
import getopt
import glob
import re
import time

def main(argv):
    inFile = ''
    taxaFile = ''
    try:
        opts, args = getopt.getopt(argv,'i:t:h',['in','taxa','help'])
    except getopt.GetoptError:
        print('python subsetMatrix.py -i phyloprofile.matrix -t taxa.list')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h','--help'):
            print('python subsetMatrix.py -i phyloprofile.matrix -t taxa.list')
            sys.exit()
        elif opt in ('-i','--in'):
            inFile = arg
        elif opt in ('-t','--taxa'):
            taxaFile = arg

    ## get taxa list
    taxaList = ''
    with open(taxaFile, 'r') as taxaF:
        taxaList = taxaF.read().replace('\n', ';')
    taxaList = taxaList + 'geneID;'
    # print(taxaList)
    # time.sleep(3)

    ## read main matrix and get only info from taxaList
    ## print to a new output matrix file
    matrixOut = open(inFile+'.subset', 'w')
    header = []
    with open(inFile) as f:
        for line in f:
            if re.search('geneID',line):
                #print(line)
                header = re.split("\t",line.rstrip())

            item = re.split("\t",line.rstrip())
            newLine = []
            for i in range(0,len(item)):
                if re.search(header[i],taxaList):
                    # print('%s - %s' % (i,item[i]))
                    newLine.append(item[i])
                    # time.sleep(1)

            out = '\t'.join(newLine)
            matrixOut.write(out+"\n")


    matrixOut.close()

if __name__ == "__main__":
    if len(sys.argv[1:]) < 4:
        print('python subsetMatrix.py -i phyloprofile.matrix -t taxa.list')
        sys.exit(2)
    else:
        main(sys.argv[1:])
