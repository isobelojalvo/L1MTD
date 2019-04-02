import das_client
from sys import argv
import os

inputFile = argv[1] #contains a list of files that you're running over
j = 0
with open(inputFile,'r') as f:
    for x in f:
        file = open('submit-'+str(j)+'.py','w')
        j = j+1
        x = x.rstrip()
        output='\nprocess.source.secondaryFileNames = cms.untracked.vstring( \n "'+x+'") \nprocess.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:'
        x = x.split('/store',1)[-1] 
        x = '/store'+x #gets rid of any xrootd info in the filename
        das_query_string = 'dasgoclient --query="lumi file=%s" --limit=0' % x
        das_output = os.popen(das_query_string).read() #query DAS via bash
        das_output = das_output.replace("root","root',")#commas between each secondary file
        das_output = das_output.replace("/store","'root://cmsxrootd.fnal.gov//store") #need xrootd info in secondary filename
        das_output = das_output.replace("\n","") #don't want newline in string
        das_output = das_output.replace("[","") #don't want in string
        das_output = das_output.replace("]","") #don't want in string
        #print das_output
        output = output+das_output+'")\n\n'
        file.write(output)
        file.close()
        print output

