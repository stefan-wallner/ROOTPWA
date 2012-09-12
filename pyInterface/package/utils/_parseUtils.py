
import glob
import os

import pyRootPwa
import pyRootPwa.utils

def parseMassBinArgs(allMassBins, massBinArg):

	massBins = []
	if massBinArg == "all":
		massBins = allMassBins
	elif massBinArg.find("-") > 0 or massBinArg.find(",") > 0:
		rawMassBinIndices = massBinArg.split(",")
		massBinIndices = []
		for massBinIndex in rawMassBinIndices:
			if massBinIndex.find("-") > 0:
				(lb, tmp, ub) = massBinIndex.partition("-")
				try:
					lb = int(lb)
					ub = int(ub)
				except ValueError:
					return []
				for i in range(lb, ub+1):
					massBinIndices.append(i)
			else:
				try:
					mbi = int(massBinIndex)
				except ValueError:
					return []
				massBinIndices.append(mbi)
		for index in massBinIndices:
			try:
				massBins.append(allMassBins[index-1])
			except IndexError:
				pyRootPwa.utils.printErr("Mass bin command line option out of range. Aborting...")
				sys.exit(1)
	else:
		try:
			mbi = int(massBinArg)
			massBins.append(allMassBins[mbi-1])
		except ValueError:
			return []
		except IndexError:
			pyRootPwa.utils.printErr("Mass bin command line option out of range. Aborting...")
			sys.exit(1)

	return massBins

def _evtOrRoot(filename):
	if os.path.isfile(filename + ".root"):
		filename += ".root"
	elif os.path.isfile(filename + ".evt"):
		filename += ".evt"
	else:
		filename = ""
	return filename

def getListOfInputFiles(massBins):

	inputDataFiles = []
	inputPSFiles = []
	inputAccPSFiles = []

	dataFileExtensionQualifier = pyRootPwa.config.dataFileExtensionQualifier
	phaseSpaceEventFileExtenisonQualifier = pyRootPwa.config.phaseSpaceEventFileExtenisonQualifier
	accCorrPSEventFileExtensionQualifier = pyRootPwa.config.accCorrPSEventFileExtensionQualifier

	for massBin in massBins:
		inputFile = massBin + "/" + massBin.rsplit('/', 1)[-1]
		if dataFileExtensionQualifier != "":
			inputFile += "." + dataFileExtensionQualifier
		inputFile = _evtOrRoot(inputFile)
		if inputFile:
			inputDataFiles.append(inputFile)
		else:
			pyRootPwa.utils.printWarn('Mass bin "' + massBin + '" does not contain data input file "' + inputFile + '{.root/.evt}".')
		if phaseSpaceEventFileExtenisonQualifier != "":
			inputFile = massBin + "/" + massBin.rsplit('/', 1)[-1] + "." + phaseSpaceEventFileExtenisonQualifier
			inputFile = _evtOrRoot(inputFile)
			if inputFile:
				inputPSFiles.append(inputFile)
		if accCorrPSEventFileExtensionQualifier != "":
			inputFile = massBin + "/" + massBin.rsplit('/', 1)[-1] + "." + accCorrPSEventFileExtensionQualifier
			inputFile = _evtOrRoot(inputFile)
			if inputFile:
				inputAccPSFiles.append(inputFile)

	return (inputDataFiles, inputPSFiles, inputAccPSFiles)

def getListOfKeyfiles(keyfilePattern):
	keyfiles = []
	if os.path.isdir(keyfilePattern):
		keyfiles = glob.glob(keyfilePattern + "/*.key")
	elif os.path.isfile(keyfilePattern) and keyfilePattern.find(".key") > 0:
		keyfiles.append(keyfilePattern)
	else:
		globbedKeyfiles = glob.glob(keyfilePattern)
		for keyfile in globbedKeyfiles:
			if os.path.isfile(keyfile) and keyfile.find(".key") > 0:
				keyfiles.append(keyfile)
			else:
				pyRootPwa.utils.printWarn("Keyfile " + keyfile + " is not valid. Skipping...")
	return keyfiles

