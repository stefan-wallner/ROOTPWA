import numpy
import pyRootPwa.utils
import pyRootPwa.core


def _getAmplitudes(keyFileNameList, prodNames, decayNames, integralMetaData, addKeyfilecontent):
	amplitudes      = []
	waveNames       = []
	for keyFile in keyFileNameList:
		waveDescription = pyRootPwa.core.waveDescription.parseKeyFile(keyFile[0])[keyFile[1]]
		if addKeyfilecontent:
			if not integralMetaData.addKeyFileContent(waveDescription.keyFileContent()):
				pyRootPwa.utils.printWarn("could not add keyfile content. Aborting...")
				return False, False
		else:
			if not integralMetaData.hasKeyFileContent(waveDescription.keyFileContent()):
				pyRootPwa.utils.printErr("keyfile content of additional eventFiledID missing in first eventFieldID.")
				return False, False

		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			pyRootPwa.utils.printErr('could not construct amplitude for keyfile "' + keyFile[0] + '" (ID '+str(keyFile[1])+'). Aborting...')
			return False, False
		amplitude.init()
		topo = amplitude.decayTopology()
		if not topo.initKinematicsData(prodNames, decayNames):
			pyRootPwa.utils.printErr("could not initialize the decay topology with the kinematics data. Aborting...")
			return False, False
		amplitudes.append(amplitude)
		waveNames.append(pyRootPwa.core.waveDescription.waveNameFromTopology(amplitude.decayTopology()))
	return amplitudes, waveNames


def _integrate(amplitudes, eventTree, waveNames, minEvent, maxEvent, multibinBoundaries):
	prodKinMomenta  = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	eventTree.SetBranchAddress(pyRootPwa.core.eventMetadata.productionKinematicsMomentaBranchName, prodKinMomenta)
	eventTree.SetBranchAddress(pyRootPwa.core.eventMetadata.decayKinematicsMomentaBranchName, decayKinMomenta)
	integralMatrix = pyRootPwa.core.ampIntegralMatrix()
	hashers = [pyRootPwa.core.hashCalculator() for _ in range(len(amplitudes))]
	integralMatrix.setWaveNames(waveNames)
	ampWaveNameMap   = {}
	binningVariables = {}
	for key in multibinBoundaries:
		binningVariables[key] = numpy.array(1, dtype = float)
		eventTree.SetBranchAddress(key, binningVariables[key])
	for waveName in waveNames:
		ampWaveNameMap[waveName] = 0.+0.j
	pyRootPwa.utils.printInfo("starting event loop.")
	skippedEvents = 0
	progressBar = pyRootPwa.utils.progressBar(minEvent, maxEvent)
	progressBar.start()
	for evt_i in range(minEvent, maxEvent):
		progressBar.update(evt_i)
		eventTree.GetEvent(evt_i)
		skipEvent = False
		for key in multibinBoundaries:
			if binningVariables[key] < multibinBoundaries[key][0] or  binningVariables[key] >= multibinBoundaries[key][1]:
				skippedEvents += 1
				skipEvent = True
				break
		if skipEvent:
			continue
		for amp_i, amplitude in enumerate(amplitudes):
			topo = amplitude.decayTopology()
			if not topo.readKinematicsData(prodKinMomenta, decayKinMomenta):
				pyRootPwa.utils.printErr("could not load kinematics data. Aborting...")
				return None, None
			ampl = amplitude()
			hashers[amp_i].Update(ampl)
			ampWaveNameMap[waveNames[amp_i]] = ampl
		if not integralMatrix.addEvent(ampWaveNameMap):
			pyRootPwa.utils.printErr("could not add event to integral matrix. Aborting...")
			return None, None
	pyRootPwa.utils.printInfo(str(skippedEvents) + " events rejected because they are outside the binning.")
	return integralMatrix, hashers


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
	integrals = []
	for iEventFileName, eventFileName in enumerate(eventFileNames):
		eventFile = pyRootPwa.ROOT.TFile.Open(eventFileName, "READ")
		if not eventFile:
			pyRootPwa.utils.printErr("could not open event file. Aborting...")
			return False
		eventMeta  = pyRootPwa.core.eventMetadata.readEventFile(eventFile)
		prodNames  = eventMeta.productionKinematicsParticleNames()
		decayNames = eventMeta.decayKinematicsParticleNames()
		amplitudes, waveNames = _getAmplitudes(keyFileNameList, prodNames, decayNames, metadataObject, addKeyfilecontent=iEventFileName==0)
		if not amplitudes or not waveNames:
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
		integralMatrix, hashers = _integrate(amplitudes, eventTree, waveNames, minEvent, maxEvent, multibinBoundaries)
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
