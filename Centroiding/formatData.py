path = "data/"

with open(path + 'fiber_values.txt') as f:
    content = f.readlines()[1:]
    
          
with open(path + 'fiberOut.txt', 'w') as f:
    f.writelines('eventid layer signal axis fiberid x/y z energy\n')
    for line in content:
        tempStr = ""
        
        for indexer in line:
            
            if(indexer != " "):
                
                tempStr = tempStr + indexer
                
            else :
                f.writelines(tempStr + '\n')
                tempStr = ""
                
        f.writelines(tempStr)
        
        
with open(path + 'raw_geant.txt') as f:
    rawData = f.readlines()[1:]
    
with open(path + 'geantOut.txt', 'w') as f:
    f.writelines('id x y z e deposit layer material\n')
    for line in rawData:
        tempStr = ""
        
        for indexer in line:
            
            if(indexer != " "):
                
                tempStr = tempStr + indexer
                
            else :
                f.writelines(tempStr + '\n')
                tempStr = ""
                
        f.writelines(tempStr)
