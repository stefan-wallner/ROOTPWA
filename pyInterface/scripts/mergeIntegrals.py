#!/usr/bin/env python

import argparse
import sys
import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="merge integral matrices"
	                                )
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("--same-multibin", dest="sameMultibin", action="store_true",
	                    help="The integrals, which should be merged, are from the same multibin.")

	parser.add_argument("inputFileNames", type=str, metavar="inFileNames", nargs='*', help="path to input files")
	parser.add_argument("outputFileName", type=str, metavar="outFileName", help="path to output file")

	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	if os.path.exists(args.outputFileName):
		printErr("Output file '{0}' exists!")
		sys.exit(1)
	outFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "CREATE")
	if not outFile: # Do this up here. Without the output file, nothing else makes sense
		pyRootPwa.utils.printErr("could not open output file. Aborting...")
		sys.exit(1)

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)

	progressBar = pyRootPwa.utils.progressBar(0, len(args.inputFileNames))
	progressBar.start()
	metadataObject = None
	for i_inputFileName, inputFileName in enumerate(args.inputFileNames):
		inputFile = ROOT.TFile.Open(inputFileName, "READ")
		inputMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(inputFile)
		if metadataObject is None:
			metadataObject = inputMeta
		else:
			if not metadataObject.mergeIntegralMatrix(inputMeta, args.sameMultibin):
				printErr("Could not add matrix from file '{0}'".format(inputFileName))
				sys.exit(1)
		inputFile.Close()
		progressBar.update(i_inputFileName)

	if not metadataObject.writeToFile(outFile):
		pyRootPwa.utils.printErr("could not write integral objects to file. Aborting...")
		sys.exit(1)
	outFile.Close()

	printSucc("Successfully merged {0} integral files".format(len(args.inputFileNames)))
