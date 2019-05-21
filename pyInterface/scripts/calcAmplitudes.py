#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core
import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT


if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file"
	                                )

	parser.add_argument("-c", type=str, metavar="configFileName", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-b", type=int, metavar="eventFileId", default=-1, dest="eventFileId", help="event file id to be calculated (default: all)")
	parser.add_argument("-e", type=str, metavar="eventsType", default="all", dest="eventsType", help="events type to be calculated ('real', 'generated' or 'accepted', default: all)")
	parser.add_argument("-f", "--no-progress-bar", action="store_true", dest="noProgressBar", help="disable progress bars (decreases computing time)")
	parser.add_argument("-k", "--keyfileIndex", action='append', metavar="#", type=int,
	                    help="keyfile index to calculate amplitude for (overrides settings from the config file, index from 0 to number of keyfiles - 1)")
	parser.add_argument("-w", type=str, metavar="wavelistFileName", default="", dest="wavelistFileName", help="path to wavelist file (default: none)")
	parser.add_argument("--output", type=str, metavar="outputFileName", default="", dest="outputFileName",
	                    help="path to output file for all calculated amplitudes (default: from file manager)")
	args = parser.parse_args()

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	pyRootPwa.core.integralTableContainer.setDirectory(config.phaseSpaceIntegralDirectory)
	pyRootPwa.core.integralTableContainer.setUpperMassBound(config.phaseSpaceUpperMassBound)

	if not args.wavelistFileName == "" and not args.keyfileIndex == -1:
		pyRootPwa.utils.printErr("Setting both options -k and -w is conflicting. Aborting...")
		sys.exit(1)

	waveList = []
	if not args.wavelistFileName == "":
		waveList = [ i[0] for i in pyRootPwa.utils.getWaveDescThresFromWaveList(args.wavelistFileName, fileManager.getWaveDescriptions()) ]
	if args.keyfileIndex:
		allWaveNames = fileManager.getWaveNameList()
		for keyfileIndex in args.keyfileIndex:
			if not keyfileIndex < len(allWaveNames):
				pyRootPwa.utils.printErr("keyfileIndex '" + str(keyfileIndex) + "' from command line argument out of range. Maximum value is " + str(len(allWaveNames)-1) + ". Aborting...")
				sys.exit(1)
			pyRootPwa.utils.printInfo("using keyfile index " + str(keyfileIndex) + " resulting in the wave name '" + allWaveNames[keyfileIndex] + "'.")
		waveList = [allWaveNames[keyfileIndex] for keyfileIndex in args.keyfileIndex]
	if not waveList:
		waveList = fileManager.getWaveNameList()

	eventsTypes = []
	if args.eventsType == "real":
		eventsTypes = [ pyRootPwa.core.eventMetadata.REAL ]
	elif args.eventsType == "generated":
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED ]
	elif args.eventsType == "accepted":
		eventsTypes = [ pyRootPwa.core.eventMetadata.ACCEPTED ]
	elif args.eventsType == "all":
		eventsTypes = [ pyRootPwa.core.eventMetadata.REAL,
		                pyRootPwa.core.eventMetadata.GENERATED,
		                pyRootPwa.core.eventMetadata.ACCEPTED ]
	else:
		pyRootPwa.utils.printErr("Invalid events type given ('" + args.eventsType + "'). Aborting...")
		sys.exit(1)

	outputFile = None
	if args.outputFileName:
		outputFile = ROOT.TFile.Open(args.outputFileName, "NEW")

	for waveName in waveList:
		for eventsType in eventsTypes:
			eventAmpFilePairs = fileManager.getEventAndAmplitudePairPathsForWave(eventsType, waveName)
			if args.eventFileId >= 0:
				if args.eventFileId >= len(eventAmpFilePairs):
					pyRootPwa.utils.printErr("event file id out of range (" + str(args.eventFileId) + ">=" + str(len(eventAmpFilePairs)) + "). Aborting...")
					sys.exit(1)
				eventAmpFilePairs = eventAmpFilePairs[args.eventFileId:args.eventFileId+1]
			for eventFilePath, amplitudeFilePath in eventAmpFilePairs:
				if outputFile is not None:
					amplitudeFilePath = outputFile
				if not pyRootPwa.calcAmplitude(eventFilePath, waveName, fileManager.getWaveDescription(waveName),
				                               amplitudeFilePath, not args.noProgressBar):
					pyRootPwa.utils.printWarn("could not calculate amplitude.")
					raise Exception()
