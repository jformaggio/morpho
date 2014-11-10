import ROOT as root

import pystanLoad as pyL

import json
import sys
from array import array

import numpy as np

#import matplotlib
#matplotlib.use('TKAgg')
#import matplotlib.pyplot as plt

def theHack(theString,theVariable,theSecondVariable="",theThirdVariable=""):
    theResult = str(theString.format(theVariable,theSecondVariable,theThirdVariable))
    return theResult

theConfigFile = sys.argv[1]
print "Using configuration file :",theConfigFile

# Load in configuration file

json_file = open(theConfigFile).read()
json_data = json.loads(json_file)
config_file = pyL.readLabel(json_data,'stan','./run.json')

# Set up the stan model and cache directory

theModel = pyL.readLabel(config_file, 'model')
casheDirectory = pyL.readLabel(theModel, 'cache', './cache')
theModelFile =  pyL.readLabel(theModel,'file')

# Set up the input, output and plotting configurations
theDataFiles =  pyL.readLabel(config_file, 'data')
theSample = pyL.readLabel(config_file,'sample', None)
thePlots =  pyL.readLabel(config_file,'plot', None)
theOutput = pyL.readLabel(config_file,'output', None)

# Set up running conditions
theRunConditions =  pyL.readLabel(config_file,'run')
theAlgorithm =  pyL.readLabel(theRunConditions,'algorithim','NUTS')
nIter =   pyL.readLabel(theRunConditions,'iter',2000)
nChain =  pyL.readLabel(theRunConditions,'chain',4)
nThin  =  pyL.readLabel(theRunConditions,'thin',1)
nWarmup = pyL.readLabel(theRunConditions,'warmup',nIter/2)
nInitPar = pyL.readLabel(theRunConditions,'init','')

parList = pyL.readLabel(theRunConditions,'pars',None)
if parList is None:
    nPars = None
else :
    nPars  =  []
    for keys in parList:
	nPars.append(str(keys))

# Load in the data
theData = pyL.stan_data_files(theDataFiles)

# Execute the fit
theFit = pyL.stan_cache(model_code= theModelFile,
			cashe_dir= casheDirectory,
			data=theData,
			algorithm = theAlgorithm,
			iter=nIter, chains=nChain,
			thin=nThin, warmup=nWarmup,
			pars=nPars, init = nInitPar,
			sample_file=theSample)

print 'Analysis complete.'

# Output the data into a root file

if theOutput is not None:
    theOutputFile = pyL.readLabel(theOutput,'name',None)
    theOutputType = pyL.readLabel(theOutput,'format',None)
    theTreeName =  pyL.readLabel(theOutput,'tree',None)
    if theOutputType=='root':

	afile = root.TFile.Open(theOutputFile, "recreate")
	atree = root.TTree(theTreeName,theTreeName)
	theOutputVar = pyL.readLabel(theOutput,'branches', None)
	theOutputData = {}

	nBranches = len(theOutputVar)

	iBranch = 0
	for key in theOutputVar:

	    nSize = pyL.readLabel(key,'ndim',1)
	    exec(theHack("theVariable_{} = np.zeros({}, dtype=float)",str(key['root_alias']),nSize))

	    if (nSize == 1) :
		exec(theHack("atree.Branch(str(key['root_alias']), theVariable_{}, key['root_alias']+'/D')",str(key['root_alias'])))
	    else :
		exec(theHack("atree.Branch(str(key['root_alias']), theVariable_{}, key['root_alias']+'[{}]/D')",str(key['root_alias']),nSize))

	    theOutputData[iBranch] = theFit.extract(key['variable'])
	    nEvents = len(theOutputData[iBranch][key['variable']])
	    iBranch += 1

	for iEvent in range(0,nEvents):
	    iBranch = 0
	    for key in theOutputVar:
		theValue = theOutputData[iBranch][key['variable']][iEvent]
		nSize = pyL.readLabel(key,'ndim',1)
		if (nSize == 1) :
		    exec(theHack("theVariable_{}[0] = theValue",str(key['root_alias'])))
		else :
		    for iNum in range(0,nSize):
			exec(theHack("theVariable_{}[{}] = theValue[{}]",str(key['root_alias']),iNum,iNum))
		iBranch +=1
	    atree.Fill()

	atree.Write()
	afile.Close()

	print 'The root file has been written to', theOutputFile

    elif theOutputType=='R':

	theOutputVar = pyL.readLabel(theOutput,'branches', None)
	theOutputData = {}
	for key in theOutputVar:
	    theOutputData.append(theFit.extract(key['variable']))

	pystan.misc.stan_rdump(theOutputFile)

    else :
	print 'The output format is not supported:', theOutputType

# Print results and make plots

if thePlots is not None:
    for key in thePlots:
	thePlotData = theFit.extract()
	theFit.plot(key['variable'])
#    plt.show()
    print(theFit)
