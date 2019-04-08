
import os as _os

import pyRootPwa.core
import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def calcAmplitude(inputFileName,
                  waveName,
                  waveDescription,
                  outputFile,
                  printProgress = True):

	printInfo = pyRootPwa.utils.printInfo
	printSucc = pyRootPwa.utils.printSucc
	printWarn = pyRootPwa.utils.printWarn




	if 'ROOTPWA' not in _os.environ:
		printWarn("$ROOTPWA not set.")
		return False

	inputFile = ROOT.TFile.Open(inputFileName, "READ")
	if not inputFile:
		printWarn("could not open input file '" + inputFileName + "'.")
		return False
	eventMeta = pyRootPwa.core.eventMetadata.readEventFile(inputFile, True)
	if not eventMeta:
		printWarn("could not read metadata from input file '" + inputFileName + "'.")
		return False
	if isinstance(outputFile, ROOT.TFile):
		closeOutputFile = False
		outputFileName = outputFile.GetName()
	else:
		closeOutputFile = True
		outputFileName = outputFile
		outputFile = ROOT.TFile.Open(outputFileName, "NEW")
	if not outputFile:
		printWarn("could not open output file '" + outputFileName + "'.")
		return False

	printInfo("Calculating amplitude for wave '" + waveName + "'" +
	          " with input file '" + inputFileName + "', output file '" + outputFileName + "'.")

	(result, amplitude) = waveDescription.constructAmplitude()
	if not result:
		printWarn("could not construct amplitude for wave '" + waveName + "'.")
		outputFile.Close()
		return False

	ampFileWriter = pyRootPwa.core.amplitudeFileWriter()
	objectBaseName = waveDescription.waveNameFromTopology(amplitude.decayTopology())
	if not ampFileWriter.initialize(outputFile, [eventMeta], waveDescription.keyFileContent(), objectBaseName):
		printWarn("could not initialize amplitudeFileWriter.")
		outputFile.Close()
		return False
	amplitudes = pyRootPwa.core.calcAmplitude(eventMeta, amplitude, -1, printProgress)
	nEvents = eventMeta.eventTree().GetEntries()
	if not amplitudes and nEvents > 0:
		printWarn("could not calculate amplitudes.")
		outputFile.Close()
		return False
	if nEvents != len(amplitudes):
		printWarn("number of events (" + str(nEvents) +
		          ") does not match with number of amplitudes (" + str(len(amplitudes)) + ").")
		return False
	ampFileWriter.addAmplitudes(amplitudes)
	if not ampFileWriter.finalize():
		printWarn("could not finalize amplitudeFileWriter.")
		outputFile.Close()
		return False

	if closeOutputFile:
		outputFile.Close()
	printSucc("successfully calculated amplitude for " + str(nEvents) + " events.")
	return True
