import pyRootPwa.utils
import pyRootPwa.core


def _getAmplitudes(waveDescriptions):
	amplitudes      = []
	for waveDescription in waveDescriptions:
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			pyRootPwa.utils.printErr('could not construct amplitude. Aborting...')
			return False, False
		amplitudes.append(amplitude)
	return amplitudes


def calcIntegralsOnTheFly(integralFileName, eventFileNames, keyFileNameList, multibinBoundaries = None,
	                      maxNmbEvents = -1, startEvent = 0, splitEventsInBunchesN = None, splitEventsInBunchesBunch = None):
	"""
	Calculate the integral matrix of the given multibin by calculating the amplitudes on the fly.

	The integral can be calculated from the whole input file or from a subsample of events in the input file
	by the following two combinations of two arguments.
	@param maxNmbEvents: Maximal number of events per input file
	@param startEvent: Index of first event processed per input file
	@param splitEventsInBunchesN: Number of bunches, into which each input file will be splitted
	@param splitEventsInBunchesBunch: Current bunch to include in the input file
	"""

	if splitEventsInBunchesN is not None or splitEventsInBunchesBunch is not None:
		if splitEventsInBunchesN is None or splitEventsInBunchesBunch is None:
			pyRootPwa.core.printErr("Both, 'splitEventsInBunchesN' and 'splitEventsInBunchesBunch', need to be defined!")
			return False
		if maxNmbEvents != -1 or startEvent != 0:
			pyRootPwa.core.printErr("'splitEventsInBunchesN' and 'splitEventsInBunchesBunch', cannot be defined together with 'maxNmbEvents', 'startEvent'!")
			return False
		splitEventsInBunchesN = int(splitEventsInBunchesN)
		splitEventsInBunchesBunch = int(splitEventsInBunchesBunch)


	if not isinstance(eventFileNames, list):
		eventFileNames = [eventFileNames]

	outFile = pyRootPwa.ROOT.TFile.Open(integralFileName, "CREATE")
	if not outFile: # Do this up here. Without the output file, nothing else makes sense
		pyRootPwa.utils.printErr("could not open output file. Aborting...")
		return False

	metadataObject =  pyRootPwa.core.ampIntegralMatrixMetadata()
	metadataObject.setGitHash(pyRootPwa.core.gitHash())
	if "mass" in multibinBoundaries:
		if multibinBoundaries["mass"][0] > 200.:
			multibinBoundaries["mass"] = (multibinBoundaries["mass"][0]/1000.,multibinBoundaries["mass"][1]/1000.)
	metadataObject.setMultibinBoundaries(multibinBoundaries)

	waveDescriptions = []
	for keyFileName in keyFileNameList:
		waveDescriptions += pyRootPwa.core.waveDescription.parseKeyFile(keyFileName)
	for waveDescription in waveDescriptions:
		if not metadataObject.addKeyFileContent(waveDescription.keyFileContent()):
			pyRootPwa.utils.printWarn("could not add keyfile content. Aborting...")
			return False

	integrals = []
	for  eventFileName in eventFileNames:
		eventFile = pyRootPwa.ROOT.TFile.Open(eventFileName, "READ")
		if not eventFile:
			pyRootPwa.utils.printErr("could not open event file. Aborting...")
			return False
		eventMeta  = pyRootPwa.core.eventMetadata.readEventFile(eventFile)
		amplitudes = _getAmplitudes(waveDescriptions)
		if not amplitudes:
			pyRootPwa.utils.printErr("could initialize amplitudes. Aborting...")
			return False
		eventTree = eventMeta.eventTree()
		nEvents   = eventTree.GetEntries()
		minEvent = startEvent
		maxEvent = nEvents
		if splitEventsInBunchesN is not None:
			quotient = nEvents / splitEventsInBunchesN
			modulus  = nEvents % splitEventsInBunchesN
			minEvent = quotient*splitEventsInBunchesBunch     + min(splitEventsInBunchesBunch, modulus)   # inclusive
			maxEvent = quotient*(splitEventsInBunchesBunch+1) + min(splitEventsInBunchesBunch+1, modulus) # exclusive
		if maxNmbEvents	> -1:
			maxEvent = min(maxEvent, startEvent + maxNmbEvents)
		if not metadataObject.addEventMetadata(eventMeta):
			pyRootPwa.utils.printErr("could not add event metadata to integral metadata. Aborting...")
			return False
		if not multibinBoundaries:
			multibinBoundaries = eventMeta.multibinBoundaries()
			if not multibinBoundaries:
				pyRootPwa.utils.printWarn("no binning map found.")
		integralMatrix, hashers = pyRootPwa.core.calcIntegralOnTheFly(eventMeta, amplitudes, multibinBoundaries, minEvent, maxEvent)
		if integralMatrix.nmbEvents() == 0:
			continue # no events from the multibin found in this event file -> skipping it
		if not integralMatrix or not hashers:
			pyRootPwa.utils.printErr("could not integrate. Aborting...")
			return False
		for hasher in hashers:
			if not metadataObject.addAmplitudeHash(hasher.hash()):
				pyRootPwa.utils.printWarn("could not add the amplitude hash.")
				# This error is not fatal, since in special cases the same hash can appear twice:
				# e.g. in freed-isobar analyses with spin zero, the angular dependences are constant
				# and the shape is either 0 or 1. If two such waves accidentally have the same number
				# of events, both will also have the same hash.
		integrals.append(integralMatrix)
		eventFile.Close()

	integralMatrix = integrals[0]
	if len(integrals) > 1:
		for integral in integrals[1:]:
			integralMatrix += integral
	if not metadataObject.setAmpIntegralMatrix(integralMatrix):
		pyRootPwa.utils.printErr("could not add the integral matrix to the metadata object. Aborting...")
		return False
	if not metadataObject.writeToFile(outFile):
		pyRootPwa.utils.printErr("could not write integral objects to file. Aborting...")
		return False
	outFile.Close()
	return True
