#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="""
	                                                Generate phase-space Monte Carlo events and calculate the
	                                                weight of each event from the provided fit result. The events
	                                                in the output file are not deweighted.
	                                             """
	                                )

	parser.add_argument("reactionFile", type=str, metavar="reactionFile", help="reaction config file")
	parser.add_argument("fitResult", type=str, metavar="fitResult", help="fitResult to get the production amplitudes")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="output root file")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName",
	                    help="path to config file (default: '%(default)s')")
	parser.add_argument("-i", "--integralFile", type=str, metavar="integralFile", help="integral file")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=100, help="(max) number of events to generate (default: %(default)s)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random number generator seed (default: %(default)s)")
	parser.add_argument("-M", type=float, metavar="#", dest="massLowerBinBoundary",
	                    help="lower boundary of mass range in MeV (!) (overwrites values from reaction file)")
	parser.add_argument("-B", type=float, metavar="#", dest="massBinWidth", help="width of mass bin in MeV (!)")
	parser.add_argument("-u", "--auxString", type=str, metavar="#", dest="auxString", help="auxiliary string stored in metadata", default="phaseSpaceEvents")
	parser.add_argument("--massTPrimeVariableNames", type=str, dest="massTPrimeVariableNames", help="Name of the mass and t' variable (default: '%(default)s')",
	                    default="mass,tPrime")
	parser.add_argument("--noStoreMassTPrime", action="store_true", dest="noStoreMassTPrime", help="Do not store mass and t' variable of each event.")
	parser.add_argument("--beamfile", type=str, metavar="<beamFile>", dest="beamFileName", help="path to beam file (overrides values from config file)")
	parser.add_argument("--noRandomBeam", action="store_true", dest="noRandomBeam", help="read the events from the beamfile sequentially")
	parser.add_argument("--randomBlockBeam", action="store_true", dest="randomBlockBeam", help="like --noRandomBeam but with random starting position")

	args = parser.parse_args()

	# print some info
	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.instance.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	pyRootPwa.core.integralTableContainer.setDirectory(config.phaseSpaceIntegralDirectory)
	pyRootPwa.core.integralTableContainer.setUpperMassBound(config.phaseSpaceUpperMassBound)

	overrideMass = (args.massLowerBinBoundary is not None) or (args.massBinWidth is not None)

	if overrideMass and not ((args.massLowerBinBoundary is not None) and (args.massBinWidth is not None)):
		printErr("'-M' and '-B' can only be set together. Aborting...")
		sys.exit(2)

	printInfo("Setting random seed to " + str(args.seed))
	pyRootPwa.core.randomNumberGenerator.instance.setSeed(args.seed)

	generatorManager = pyRootPwa.core.generatorManager()
	if args.beamFileName is not None:
		generatorManager.overrideBeamFile(args.beamFileName)

	if not generatorManager.readReactionFile(args.reactionFile):
		printErr("could not read reaction file. Aborting...")
		sys.exit(1)

	if overrideMass:
		generatorManager.overrideMassRange(args.massLowerBinBoundary / 1000., (args.massLowerBinBoundary + args.massBinWidth) / 1000.)
	if args.noRandomBeam:
		generatorManager.readBeamfileSequentially()
	if args.randomBlockBeam:
		generatorManager.readBeamfileSequentially()
		generatorManager.randomizeBeamfileStartingPosition()

	if not generatorManager.initializeGenerator():
		printErr("could not initialize generator. Aborting...")
		sys.exit(1)

	massRange = generatorManager.getGenerator().getTPrimeAndMassPicker().massRange()
	massBinCenter = (massRange[0] + massRange[1]) / 2.
	fitResult = pyRootPwa.utils.getBestFitResultFromFile(fitResultFileName = args.fitResult,
	                                                     massBinCenter = massBinCenter,
	                                                     fitResultTreeName = config.fitResultTreeName,
	                                                     fitResultBranchName = config.fitResultBranchName)
	if not fitResult:
		printErr("could not find fit result in file '" + args.fitResult +
		         "' for mass bin " + str(massBinCenter) + ". Aborting...")
		sys.exit(1)

	waveNames = fitResult.waveNames()
	if 'flat' in waveNames:
		waveNames.remove('flat')  # ignore flat wave

	model = pyRootPwa.core.modelIntensity(fitResult)
	for waveName in waveNames:
		waveDescription = fileManager.getWaveDescription(waveName)
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			printErr('could not construct amplitude for wave "' + waveName + '".')
			sys.exit(1)
		if not model.addDecayAmplitude(amplitude):
			printErr('could not add amplitude for wave "' + waveName + '".')

	# overwrite integral matrix from fit result with one read from a file
	if args.integralFile:
		integralFile = pyRootPwa.ROOT.TFile.Open(args.integralFile, "READ")
		integralMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(integralFile)
		integral = integralMeta.getAmpIntegralMatrix()
		model.loadPhaseSpaceIntegral(integral)

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFile, "NEW")
	if not outputFile:
		printErr("could not open output file. Aborting...")
		sys.exit(1)

	fileWriter = pyRootPwa.core.eventFileWriter()

	printInfo("opened output root file: " + args.outputFile)
	try:
		printInfo(generatorManager)
		progressBar = pyRootPwa.utils.progressBar(0, args.nEvents, sys.stdout)
		progressBar.start()
		massTPrimeVariables = []
		attempts = 0
		decayKin = None
		prodKin = None
		first = True
		weight = 0.

		eventsGenerated = 0
		for eventsGenerated in range(args.nEvents):

			attempts += generatorManager.event()
			generator = generatorManager.getGenerator()
			beam = generator.getGeneratedBeam()
			finalState = generator.getGeneratedFinalState()

			if first:
				if not args.noStoreMassTPrime:
					if len(args.massTPrimeVariableNames.split(',')) == 2:
						massTPrimeVariables = args.massTPrimeVariableNames.split(',')
					else:
						printErr("Option --massTPrimeVariableNames has wrong format '" + args.massTPrimeVariableNames + "'. Aborting...")
						sys.exit(1)

				prodKinNames = [ beam.name ]
				decayKinNames = [ particle.name for particle in finalState]
				success = fileWriter.initialize(
				                                outputFile,
				                                args.auxString,
				                                pyRootPwa.core.eventMetadata.REAL,
				                                prodKinNames,
				                                decayKinNames,
# TODO: FILL THESE
				                                { "mass": (massRange[0], massRange[1]) },
				                                massTPrimeVariables + ["weight"]
				                                )
				if not success:
					printErr('could not initialize file writer. Aborting...')
					sys.exit(1)

				if not model.initDecayAmplitudes(prodKinNames, decayKinNames):
					printErr('could not initialize kinematics Data. Aborting...')
					sys.exit(1)
				first = False

			prodKin = [ beam.lzVec.Vect() ]
			decayKin = [ particle.lzVec.Vect() for particle in finalState ]

			weight = model.getIntensity(prodKin, decayKin)

			additionalVariables = []
			if not args.noStoreMassTPrime:
				additionalVariables = [ generator.getGeneratedXMass(), generator.getGeneratedTPrime() ]
			fileWriter.addEvent(prodKin, decayKin, additionalVariables + [weight])

			progressBar.update(eventsGenerated)
	except:
		raise
	finally:
		fileWriter.finalize()

	eventsGenerated += 1
	printSucc("generated " + str(eventsGenerated) + " events.")
	printInfo("attempts: " + str(attempts))
	printInfo("efficiency: " + str(100. * (float(eventsGenerated) / float(attempts))) + "%")

	pyRootPwa.utils.printPrintingSummary(pyRootPwa.utils.printingCounter)
