from argparse import ArgumentParser
from decimal import * # needed to avoid weird rounding errros/bias
from pprint import pprint
import os
import subprocess
import time

getcontext().prec = 5 # setting decimal precision to 5 decimals

# functions to read through Slim output
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

def isMutNotFound(line):
    if(line.startswith("\"No mutation ")):
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
rList = ['1e-10','1e-8','1e-6','1e-4']
cList = ['0.25', '0.5', '0.75', '1.0']

domCoeffList = ['0.5']



timeStr = time.strftime("%Y%m%d_%H%M%S")
modelStr = 'pseudoEntropySimsNegFreqDepSweep_neutral_sel_n400'
# prepare file to store Slim output
slimFileString = modelStr+ "_entropy_original_slim" + ".txt"
slimFile = open(slimFileString, 'w')

# prepare file to store mutation frequency output
freqsFileString = modelStr + "_entropy_freqs_run" + ".txt"
freqsFile = open(freqsFileString, 'w')
freqsFile.write('dominance_coeff\tscaling_constant\trecombination_rate\trun\n')

# prepare file to store pseudo-entropy output
hFileString = modelStr + "_entropy_delta_entropy" + ".txt"
hFile = open(hFileString, 'w')
hFile.write('dominance_coeff\tscaling_constant\trecombination_rate\trun\tnum_trajectories\tgeneration\tpseudo_entropy\tdelta_pseudo_entropy\n')


# read in template script
with open(args.input_file, 'r') as file:
    templateFile = file.read()
#templateFile = open(args.input_file, 'r').read()

print("Total number of runs to be executed: " + str(numRuns*len(rList)*len(cList)*len(domCoeffList)))

# loop through list of dominance coefficients
for currentDomCoeff in domCoeffList:
    print("Dominance coefficient " + currentDomCoeff)

    # loop through list of scaling constants
    for currentC in cList:
        print("Scaling constant " + currentC)

    # loop through list of recombination rates
        for currentR in rList:
            print("Recombination rate " + currentR)
            # replace place holder in template file
            replacedFile = templateFile.replace('$1', currentR)
            replacedFile = replacedFile.replace('$2', currentC)
            replacedFile = replacedFile.replace('$3', currentDomCoeff)

            # write new simulation file with replaced placeholder
            simulationFileName = args.input_file.split('.')[0] + '_' + currentDomCoeff + '_'+ currentR + '_' + currentC + '.' + args.input_file.split('.')[1]
            simulationFile = open(simulationFileName, 'w')
            simulationFile.write(replacedFile)
            simulationFile.close()

            # loop that starts replicate runs with the same parameters
            for i in range(numRuns):
                runNum = str(i+1)
                print("Run number " + runNum)
                mutNotFound = False

                # write new job shell script
                jobScriptName = simulationFileName + '_job.sh'
                jobScript = open(jobScriptName, 'w')
                jobScript.write("/vol/biotools/bin/slim " + simulationFileName)
                jobScript.close()

                # build command to run SLiM simulation
                command = 'sh ' + jobScriptName

                # save output stream from slim run to work with
                output_stream = subprocess.run(command, capture_output=True, shell=True, text=True)
                output = output_stream.stdout

                # save original slim output to file
                slimFile.write('Dominance coefficient: ' + currentDomCoeff + '\n')
                slimFile.write('Scaling constant: ' + currentC + '\n')
                slimFile.write('Recombination rate: ' + currentR + '\n')
                slimFile.write('Run number: ' + str(runNum) + '\n')
                for line in output:
                    slimFile.write(line)
                slimFile.write('\n')

                # process data from slim output
                generationDict = {}
                sampledGens = []
                for line in output.split('\n'):
                    if(isMutNotFound(line)):
                        print("mut not found")
                        mutNotFound=True
                        print(str(mutNotFound))

                    if(isMutNotFound(line)):        # check if mutation with desired frequency has been found
                        break

                    elif(isGenLine(line)):
                        # save current generation as integer
                        currentGen = int(line[13:-1])
                        sampledGens.append(currentGen)

                    elif(isIdLine(line)):
                        # save current mutation IDs as list of strings
                        currentIds = line.split()[1:]

                    elif(isFreqLine(line)):
                        # save current mutation frequencies as list of 'decimals'
                        currentFreqs = list(map(Decimal,line.split()[1:]))
                        
                        # save current mutation IDs and frequencies in dictionary
                        currentFreqDict = {}
                        for indx in range(len(currentIds)):
                            currentFreqDict.setdefault(currentIds[indx], currentFreqs[indx])
                        # save the current generation with the current frequency dictionary in the main dictionary
                        generationDict.setdefault(currentGen, currentFreqDict)

                # if mutation was not found we restart the loop iteration
                if(mutNotFound==True):
                    i-=1
                else:

                    #### clean up data ######
                    t_0_Dict = generationDict.get(10000)
                    if(t_0_Dict is None):
                        continue

                    for id in t_0_Dict.keys():  # loop through mutation ids (str)
                        for gen in generationDict.keys():   # loop through generation numbers (int)
                            currentGenDict = generationDict.get(gen)
                            if(id not in currentGenDict.keys()):
                                generationDict.get(gen).setdefault(id, Decimal(0.0))    # insert frequencies of 0.0 for lost mutations
                            for k in list(currentGenDict.keys()):
                                if(k not in t_0_Dict.keys()):
                                    currentGenDict.pop(k)   # remove mutation that was not in population at the start of sampling period
                    
                    # save mutation frequencies of run to output file
                    for id in t_0_Dict.keys():      # loops through mutation ids found in population
                        freqsFile.write(currentDomCoeff + '\t')
                        freqsFile.write(currentC + '\t')
                        freqsFile.write(currentR + '\t')
                        freqsFile.write(str(runNum) + '\t')
                        freqsFile.write(id + '\t')
                        for gen in generationDict.keys():       # loops through sampled generations
                            currentGenDict = generationDict.get(gen)
                            currentFreq = currentGenDict.get(id, "NA")
                            currentStr = '{x:.4f}'.format(x = currentFreq) + '\t'
                            freqsFile.write(currentStr) 
                        freqsFile.write('\n')
                    
                    # calculate delta and save it to output file
                    h_0 = calculateH(generationDict.get(10000))
                    for gen in sampledGens:
                        hFile.write(currentDomCoeff + '\t')
                        hFile.write(currentC + '\t')
                        hFile.write(currentR + '\t')
                        hFile.write(str(runNum) + '\t')
                        hFile.write(str(len(generationDict.get(gen))) + '\t')
                        hFile.write(str(gen) + '\t')
                        h_t = calculateH(generationDict.get(gen))
                        hFile.write('{x:.4f}'.format(x = h_t))
                        hFile.write('\t')
                        delta_H = '{x:.4f}'.format(x = h_t - h_0)
                        hFile.write(str(delta_H)+ '\n')
                
            # remove simulation file after use
            # remove job script after use
            os.remove(jobScriptName)
            os.remove(simulationFileName)
    
hFile.close()
slimFile.close()
freqsFile.close()
print('Skript finished!')
