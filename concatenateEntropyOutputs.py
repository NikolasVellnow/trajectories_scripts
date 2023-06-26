from argparse import ArgumentParser
import time
import pprint
import os


parser = ArgumentParser()
parser.add_argument("-p", "--path", required=False, default=None)
args = parser.parse_args()

path = args.path

timeStr = time.strftime("%Y%m%d_%H%M%S")

run=0
generations = list(range(10000,10055,5))

concatData = 'run\tgeneration\tpseudoH\tdelta_pseudoH\n'

# read in all output files
for filename in os.listdir(path):
    if filename.endswith('delta_entropy.txt'):
        strList = filename.split("_")
        run = strList[2]
        currentFile = open(filename, 'r')
        genCounter = 0
        for line in currentFile:
            if not line.startswith('H'):
                concatData = concatData + str(run) + '\t' + str(generations[genCounter]) + '\t' + line
                genCounter += 1
        currentFile.close()

fileString = timeStr + "_concatenated_delta_entropy.txt"
outputFile = open(fileString, 'w')
outputFile.write(concatData)
outputFile.close()
