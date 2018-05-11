#include "calcIntegral.h"

#include <complex>
#include <map>
#include <boost/progress.hpp>

#include <TTree.h>
#include <TClonesArray.h>

#include <reportingUtils.hpp>

using namespace std;
using namespace rpwa;

std::pair<rpwa::ampIntegralMatrix, std::vector<rpwa::hashCalculator> >
rpwa::hli::calcIntegralOnTheFly(const rpwa::eventMetadata& eventMeta,
                                const std::vector<rpwa::isobarAmplitudePtr>& amplitudes,
                                rpwa::multibinBoundariesType& multibinBoundaries,
                                const long int minEvent,
                                long int maxEvent,
                                const long int treeCacheSize)
                                {

	TTree* eventTree = eventMeta.eventTree();
	TBranch* prodKinMomentaBr = 0;
	TBranch* decayKinMomentaBr = 0;
	TClonesArray* prodKinMomenta = 0;
	TClonesArray* decayKinMomenta = 0;

	// connect leaf variables to tree branches
	eventTree->SetBranchAddress(eventMetadata::productionKinematicsMomentaBranchName.c_str(), &prodKinMomenta, &prodKinMomentaBr);
	eventTree->SetBranchAddress(eventMetadata::decayKinematicsMomentaBranchName.c_str(), &decayKinMomenta, &decayKinMomentaBr);
	eventTree->SetCacheSize(treeCacheSize);
	eventTree->AddBranchToCache(eventMetadata::productionKinematicsMomentaBranchName.c_str(), true);
	eventTree->AddBranchToCache(eventMetadata::decayKinematicsMomentaBranchName.c_str(), true);
	eventTree->StopCacheLearningPhase();
	if (not prodKinMomenta or not decayKinMomenta) {
		printErr<< "at least one of the input data arrays is a null pointer: "
		<< "        production kinematics: " << "momenta = " << prodKinMomenta << endl
		<< "        decay kinematics:      " << "momenta = " << decayKinMomenta << endl
		<< "skipping event." << endl;
		throw;
	}

	// initialize amplitude and decay topology
	std::vector<std::string> waveNames;
	waveNames.reserve(amplitudes.size());
	for (const isobarAmplitudePtr& amp : amplitudes) {
		amp->init();
		const isobarDecayTopologyPtr& topo = amp->decayTopology();
		if (not topo->initKinematicsData(eventMeta.productionKinematicsParticleNames(), eventMeta.decayKinematicsParticleNames())) {
			printErr<< "problems initializing input data. cannot read input data." << endl;
			throw;
		}
		waveNames.push_back(waveDescription::waveNameFromTopology(*topo));
	}

	additionalTreeVariables binningVariables;
	binningVariables.setBranchAddresses(eventMeta);

	ampIntegralMatrix integralMatrix;
	vector<hashCalculator> hashers(amplitudes.size());
	integralMatrix.setWaveNames(waveNames);

	map<string, complex<double>> ampWaveNameMap;
	for (const auto& waveName : waveNames)
		ampWaveNameMap[waveName] = 0.0;

	printInfo<< "starting event loop." << endl;
	long int nmbEvents = eventTree->GetEntries();
	if (maxEvent < 0)
		maxEvent = nmbEvents;
	maxEvent = min(maxEvent, nmbEvents);
	long int skippedEvents = 0;
	boost::progress_display* progressIndicator = new boost::progress_display(maxEvent - minEvent, cout, "");
	for (long int event_i = minEvent; event_i < maxEvent; ++event_i) {
		++(*progressIndicator);
		eventTree->GetEntry(event_i);
		if (not binningVariables.inBoundaries(multibinBoundaries)) {
			++skippedEvents;
			continue;
		}

		for (size_t amp_i = 0; amp_i < amplitudes.size(); ++amp_i) {
			const isobarAmplitudePtr& amp = amplitudes[amp_i];
			const isobarDecayTopologyPtr& topo = amp->decayTopology();
			if (not topo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {
				printErr<< "could not load kinematics data. Aborting..." << endl;
				throw;
			}
			complex<double> ampValue = amp->amplitude();
			hashers[amp_i].Update(ampValue);
			ampWaveNameMap[waveNames[amp_i]] = ampValue;
		}

		if (not integralMatrix.addEvent(ampWaveNameMap)) {
			printErr<< "could not add event to integral matrix. Aborting..." << endl;
			throw;
		}
	}

	printInfo<< skippedEvents << " events skipped, because they are outside the binning." << endl;

	return make_pair(integralMatrix, hashers);
}
