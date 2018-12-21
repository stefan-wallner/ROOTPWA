#!/usr/bin/env python

import argparse
import sys
import os
import numpy as np

from pyRootPwa.utils import progressBar
import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Calculate weights for accepted (generated) events"
	                                )

	parser.add_argument("fitResult", type=str, metavar="fitResult", help="fitResult to get the production amplitudes")
	parser.add_argument("outputFolder", type=str, metavar="outputFolder", help="output folder")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="integralBin", default=0, help="integral-bin ID of fit (default: 0)")
	parser.add_argument("--generatedEvents", action="store_true", dest="generatedEvents", help="Calculate weights for generated events instead of accepted events")
	parser.add_argument("--outTreeName", dest="outTreeName", default=None, help="Output-tree name (default: as defined in the config file)")
	parser.add_argument("--removeZeroWeightsfiles", action="store_true", dest="removeZeroWeightsfiles",
	                    help="Remove the created output files for which  all weights are zero")
	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	if not os.path.isdir(args.outputFolder):
		os.makedirs(args.outputFolder)

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	if args.integralBin < 0:
		pyRootPwa.utils.printErr("bin < 0 (" + str(args.integralBin) + "). Aborting...")
		sys.exit(1)
	elif args.integralBin >= len(fileManager.binList):
		pyRootPwa.utils.printErr("bin out of range (" + str(args.integralBin) + ">=" + str(len(fileManager.binList)) + "). Aborting...")
		sys.exit(1)

	multiBin = fileManager.binList[args.integralBin]

	# get the matching fit resulst
	resultsinMultibins = pyRootPwa.core.getFitResultsFromFilesInMultibins([args.fitResult], config.fitResultTreeName, config.fitResultBranchName,
	                                                                      True, True, True)
	if multiBin in resultsinMultibins:
		fitResult = resultsinMultibins[multiBin][0]
	else:
		printErr("Cannot find multibin '{0}' in fit-result file '{1}'".format(multiBin, args.fitResult))
		sys.exit(1)

	waveNames = fitResult.waveNames() # determines the indexing
	waveIndexFlat = waveNames.index('flat')
	nWaves = len(waveNames)
	spinDensityMatrix = np.empty((nWaves,nWaves), dtype=np.complex)
	for i in xrange(nWaves):
		for j in xrange(nWaves):
			spinDensityMatrix[i,j] = fitResult.spinDensityMatrixElem(waveNames[i],waveNames[j])

	# get normalization integral
	psIntegralPath  = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.GENERATED)
	intFile = pyRootPwa.ROOT.TFile.Open(psIntegralPath, "READ")
	intMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(intFile)
	intMatrix = intMeta.getAmpIntegralMatrix()

	# determine decay amplitude normalization factors
	decayAmplitudeNormalizationFactors = np.empty(nWaves, dtype=np.float)
	for iWave in xrange(nWaves):
		if iWave != waveIndexFlat:
			decayAmplitudeNormalizationFactors[iWave] = 1.0/np.sqrt( intMatrix.element(waveNames[iWave], waveNames[iWave]).real )
		else:
			decayAmplitudeNormalizationFactors[iWave] = 1.0

	eventType = pyRootPwa.core.eventMetadata.GENERATED if args.generatedEvents else pyRootPwa.core.eventMetadata.ACCEPTED
	eventAndAmpFileDict = fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, eventType)
	if not eventAndAmpFileDict:
		printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)

	for dataFilename in sorted(eventAndAmpFileDict.keys()):
		printInfo("Calculating weights for '{0}'".format(dataFilename))
		# open event file
		eventFile, eventMeta = pyRootPwa.utils.openEventFile(dataFilename)
		binningVariables = {} # TODO: use the additionalVariables class when merged
		for variable in multiBin.boundaries.keys():
			if variable not in eventMeta.additionalTreeVariableNames():
				printErr("Binning variable '{0}' is not in events additional saved variabled".format(variable))
				sys.exit(1)
			binningVariables[variable] = np.zeros(1, dtype=np.float)
			eventMeta.eventTree().SetBranchAddress(variable, binningVariables[variable])

		# open output file
		outFilename = os.path.join(args.outputFolder, "{0}_{1}_weights.root".format(args.integralBin, os.path.splitext(os.path.basename(dataFilename))[0]))
		if os.path.exists(outFilename):
			printErr("Output file exists: '{0}'".format(outFilename))
			sys.exit(1)
		outFile = ROOT.TFile(outFilename, "CREATE")
		weight = np.zeros(1, dtype=np.float)
		outTreeName = args.outTreeName if args.outTreeName is not None else config.weightTreeName
		outTree = ROOT.TTree(outTreeName, outTreeName)
		outTree.Branch("weight", weight, "weight/D")

		# open amplitude files
		ampMetas = {}
		ampFiles = {}
		for waveName, ampFileName in eventAndAmpFileDict[dataFilename].iteritems():
			if not waveName in waveNames: # amplitudes for this wave are not needed
				continue
			if not ampFileName in ampFiles:
				ampFiles[ampFileName] = ROOT.TFile.Open(ampFileName, "READ")
			ampFile = ampFiles[ampFileName]
			if not ampFile:
				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'. Aborting...")
				sys.exit(1)
			ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
			if not ampMeta:
				pyRootPwa.utils.printErr("could not read metadata from amplitude file '" + ampFileName + "'. Aborting...")
				sys.exit(1)
			ampMetas[ampMeta.objectBaseName()] = ampMeta

		# get and check number of events
		nEvents = eventMeta.eventTree().GetEntries()
		for ampMeta in ampMetas.itervalues():
			if ampMeta.amplitudeTree().GetEntries() != nEvents:
				printErr("Amplitude files have different number of events")
				sys.exit(1)

		if len(ampMetas) != nWaves - 1: # -1 because of flat wave
			printErr("Different number of wave names in fit result and amplitude files!")
			sys.exit(1)

		amplitudeLeaves = [None] * nWaves
		amplitudeTrees= [None] * nWaves
		for waveIndex, waveName in enumerate(waveNames):
			if waveName == 'flat':
				continue
			if not waveName in ampMetas:
				printErr("Cannot find wave '{0}' in amplitude files".format(waveName))
				sys.exit(1)

			amplitudeTrees[waveIndex] = ampMetas[waveName].amplitudeTree()
			amplitudeLeaves[waveIndex] = pyRootPwa.core.amplitudeTreeLeaf()
			amplitudeLeaves[waveIndex].setBranchAddress(amplitudeTrees[waveIndex], ampMeta.amplitudeLeafName)

		# loop over all events
		eventsInIntegralMultibin = False
		waveIndices = [i for i in xrange(nWaves)]
		pbar = progressBar(maximum=nEvents-1)
		pbar.start()
		for iEvent in xrange(nEvents):
			eventMeta.eventTree().GetEntry(iEvent)
			if multiBin.inBin(binningVariables): # @todoTODO: use c++ function when merged
				# build decayAmplitudes
				decayAmplitudes = np.empty(nWaves, dtype=np.complex)
				for iWave in xrange(nWaves):
					ampl = 0.0j
					if amplitudeLeaves[iWave] is not None: # if it is not the flat wave
						amplitudeTrees[iWave].GetEntry(iEvent)
						ampl = amplitudeLeaves[iWave].incohSubAmp(0)
						# normalization of decay amplitude
						ampl *= decayAmplitudeNormalizationFactors[iWave]
					decayAmplitudes[iWave] = ampl
				decayAmplitudesConj = np.conj(decayAmplitudes)

				# build the amplitude matrix
				# waves with different reflectivity are handled by the spin-density matrix
				decayAmplitudeMatrix = np.empty((nWaves, nWaves), dtype=np.complex)
				for i in xrange(nWaves):
					decayAmplitudeMatrix[i,:] = decayAmplitudes[i] * decayAmplitudesConj[:]
				decayAmplitudeMatrix[waveIndexFlat, waveIndexFlat] = 1.0
				# calculate the weight
				weight[0] = pyRootPwa.core.spinDensityMatrixTimesAmplitudeMatrix(waveIndices, spinDensityMatrix, decayAmplitudeMatrix)
				eventsInIntegralMultibin = True
				# normalization of weights -> sum_{i=0}^{NaccEvents} w_i = NaccEvents
				weight[0] = weight[0] / intMatrix.nmbEvents()
			else: # event outside the integral bin
				weight[0] = 0.0
			outTree.Fill()
			pbar.update(iEvent)

		for ampFile in ampFiles.values():
			ampFile.Close()
		eventFile.Close()
		outFile.Write()
		outFile.Close()

		if not eventsInIntegralMultibin and args.removeZeroWeightsfiles:
			printInfo("No events from this data file in the integral multibin")
			printInfo("Deleting the weight file")
			os.remove(outFilename)
