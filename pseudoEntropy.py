from argparse import ArgumentParser
from decimal import * # needed to avoid weird rounding errros/bias
from pprint import pprint
import subprocess
import time

getcontext().prec = 5 # setting decimal precision to 5 decimals

def isGenLine(line):
    if(line.startswith("\"generation: ")):
        return True
    else:
        return False

def isIdLine(line):
    if(line.startswith("MutationIds: ")):
        return True
    else:
        return False

def isFreqLine(line):
    if(line.startswith("MutationFreqs: ")):
        return True
    else:
        return False

# calcultes pseudo-entropy from dictionary with mutation frequencies
def calculateH(mutFreqs):
    l = len(mutFreqs)
    result = Decimal(0.0)
    sum = Decimal(0.0)
    for v in mutFreqs.values():
        freq = Decimal(v)
        if (freq > Decimal(0.0) and freq < Decimal(1.0)):
            sum += freq * Decimal(freq).ln() + (1-freq) * Decimal(1-freq).ln()
    result = -sum/l
    return result

parser = ArgumentParser()
parser.add_argument("-i", "--input_file", required=True)
parser.add_argument("-n", "--num_runs", required=False, type=int, default=1)
args = parser.parse_args()

numRuns = args.num_runs

print("Total number or runs to be executed: " + str(numRuns))
for i in range(numRuns):
    runNum = str(i+1)
    print("Run number " + runNum)
    timeStr = time.strftime("%Y%m%d_%H%M%S") + '_' + str(i+1)
    # build command to run SLiM simulation
    command = "slim " + args.input_file

    # save output stream from slim run to work with
    output_stream = subprocess.run(command, capture_output=True, shell=True, text=True)
    output = output_stream.stdout

    # save slim out put to file
    fileString = timeStr + "_entropy_rec_0.5_original_slim" + ".txt"
    slimFile = open(fileString, 'w')
    for line in output:
        slimFile.write(line)
    slimFile.close()

    # process data from slim output
    generationDict = {}
    sampledGens = []
    for line in output.split('\n'):
        if(isGenLine(line)):
            # save current generation as integer
            currentGen = int(line[13:-1])
            sampledGens.append(currentGen)

        if(isIdLine(line)):
            # save current mutation IDs as list of strings
            currentIds = line.split()[1:]

        if(isFreqLine(line)):
            # save current mutation frequencies as list of 'decimals'
            currentFreqs = list(map(Decimal,line.split()[1:]))
            
            # save current mutation IDs and frequencies in dictionary
            currentFreqDict = {}
            for indx in range(len(currentIds)):
                currentFreqDict.setdefault(currentIds[indx], currentFreqs[indx])
            # save the current generation with the current frequency dictionary in the main dictionary
            generationDict.setdefault(currentGen, currentFreqDict)

    #### clean up data ######
    t_0_Dict = generationDict.get(10000)

    for id in t_0_Dict.keys():  # loop through mutation ids (str)
        for gen in generationDict.keys():   # loop through generation numbers (int)
            currentGenDict = generationDict.get(gen)
            if(id not in currentGenDict.keys()):
                generationDict.get(gen).setdefault(id, Decimal(0.0))    # insert frequencies of 0.0 for lost mutations
            for k in list(currentGenDict.keys()):
                if(k not in t_0_Dict.keys()):
                    currentGenDict.pop(k)   # remove mutation that was not in population at the start of sampling period
    
    # save mutation frequencies of run to output file
    fileString = timeStr + "_entropy_rec_0.5_freqs_run" + ".txt"
    freqsFile = open(fileString, 'w')

    for id in t_0_Dict.keys():
        freqsFile.write(id + '\t')
        for gen in generationDict.keys():
            currentGenDict = generationDict.get(gen)
            currentFreq = currentGenDict.get(id, "NA")
            currentStr = '{x:.4f}'.format(x = currentFreq) + '\t'
            freqsFile.write(currentStr) 
        freqsFile.write('\n')
    freqsFile.close()

    # calculate delta and save it to output file
    fileString2 = timeStr + "_entropy_rec_0.5_delta_entropy" + ".txt"
    hFile = open(fileString2, 'w')

    h_0 = calculateH(generationDict.get(10000))
    currentStr = 'H(0): {x:.4f}'.format(x = h_0)
    hFile.write(currentStr)
    hFile.write('\n')
    for gen in sampledGens:
        h_t = calculateH(generationDict.get(gen)) # !!! columns not rows in freqs file
        hFile.write('{x:.4f}'.format(x = h_t))
        hFile.write('\t')
        delta_H = '{x:.4f}'.format(x = h_t - h_0)
        hFile.write(str(delta_H))
        hFile.write('\n')
    hFile.close()
    
    

