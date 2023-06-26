

def isMutNotFound(line):
    if(line.startswith("\"No mutation ")):
        return True
    else:
        return False




output = 'line 1' + '\n' + 'line 2' + '\n' + 'line 3' + '\n' + '\"No mutation that fullfills condition found!' + '\n' + 'line 5' + '\n'

for line in output.split('\n'):

    if(isMutNotFound(line)):
        break
    
    print(line)