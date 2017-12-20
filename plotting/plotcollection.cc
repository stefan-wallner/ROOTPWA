#include "reportingUtils.hpp"
#include "hashCalculator.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>

#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include "plotcollection.h"

using namespace std;
using namespace rpwa;

namespace {
	double phaseSpectrumDiscontinuity(componentPlot::plotType* plot);
	std::set<std::string> buildAllTotals(const std::vector<std::string>& waveNames);
	std::set<std::string> extractTotalPatternsFromWavename(std::string waveName, const bool oldNaming = false);
}


const std::vector<std::string> rpwa::plotComponents::names = {
        "massIndependent"
};


rpwa::multibinPlots::multibinPlots():
		_initialized(false),
		_negLogLikeSpectrum(nullptr){
}


rpwa::multibinPlots::multibinPlots(const std::vector<rpwa::fitResult>& fitresults, const std::string& label, const std::string& description) :
		_initialized(true),
		_negLogLikeSpectrum(nullptr){

	if (fitresults.size() == 0) {
		printErr<< "fitresults is empty" << std::endl;
		return;
	}

	_metadata.descriptions.push_back(description);
	_metadata.labels.push_back(label);

	std::map<double, std::vector<const fitResult*>> fitResultPtrsInMassbins;
	// sort results by mass bin
	for(const auto& fitresult: fitresults) {
		const double massBinCenter = fitresult.massBinCenter();
		fitResultPtrsInMassbins[massBinCenter].push_back(&fitresult);
	}

	// sort fitResults by convergence, neg log-like, and has valid cov. matrix
	for(auto& massFitResults: fitResultPtrsInMassbins) {
		const double massBinCenter = massFitResults.first;
		std::vector<const fitResult*>& fitResultPtrs = massFitResults.second;
		std::sort(fitResultPtrs.begin(), fitResultPtrs.end(),
				[](const fitResult* a, const fitResult* b) -> bool {
					if(a->converged()) {
						if( not b->converged()) return true;
						if(a->logLikelihood() < b->logLikelihood()) return true;
						if(a->logLikelihood() == b->logLikelihood() and a->covMatrixValid() and not b->covMatrixValid()) return true;
					} else {
						if((not b->converged()) and a->logLikelihood() < b->logLikelihood()) return true;
					}
					return false;
				});
		_fitResultsInMassbins[massBinCenter].reserve(fitResultPtrs.size());
		for(const auto& fitResultPtr: fitResultPtrs) _fitResultsInMassbins[massBinCenter].push_back(*fitResultPtr);
	}

	// check if the best fit result has the integral and covariance matrix
	for(const auto& massbin_results: _fitResultsInMassbins) {
		const std::vector<fitResult>& results = massbin_results.second;
		if(not results[0].converged()) {
			printErr << "No converged solution found in massbin: " << results[0].massBinCenter() << std::endl;
			throw std::invalid_argument("No converged solution found");
		}
		if (not results[0].hasHessian() or not results[0].covMatrixValid()) {
			printErr << "No solution with valid hessian matrix found in massbin: " << results[0].massBinCenter() << std::endl;
			throw std::invalid_argument("Best solution without valid covariance matrix");
		}
		if (results[0].normIntegralMatrix().nRows()==0) {
			printErr << "Best solution without integral matrix in massbin: " << results[0].massBinCenter() << std::endl;
			throw std::invalid_argument("Best solution without integral matrix");
		}
		const double minLikelihood = std::min_element(results.begin(), results.end(),
				[](const fitResult& a, const fitResult& b)->bool {return a.logLikelihood() < b.logLikelihood();})->logLikelihood();
		if ( minLikelihood != results[0].logLikelihood()) {
			printWarn << "Fit result with best log-likelihood is not converged and will be ignorded in massbin: " << results[0].massBinCenter() << " !" << std::endl;
		}
	}

	// build wave names
	std::set<std::string> waveNamesSet;
	// find all waves
	for(const auto& fitResult: fitresults) {
		for(const auto& waveName: fitResult.waveNames()) {
			waveNamesSet.insert(waveName);
		}
	}
	for(const auto& waveName: waveNamesSet) {
		_metadata.waveNames.push_back(waveName);
	}
}


void
rpwa::multibinPlots::buildDefaultPlots(const std::vector<std::string> waveNamePatterns, const size_t nRefWaves) {
	// build intensity plots
	for (const auto& waveName : _metadata.waveNames) {
		intensitySpectrum(waveName); // triggers generation of intensity spectrum
	}
	std::set<std::string> allWaveNamePatterns = buildAllTotals(_metadata.waveNames);
	allWaveNamePatterns.insert(waveNamePatterns.begin(), waveNamePatterns.end());
	for (const auto& wavenamePattern : allWaveNamePatterns) {
		intensitySpectrumRegEx(wavenamePattern);
	}

	std::vector<std::string> waveNamesOrderedByIntensity = waveNamesSortedByIntensity();
	// build phase plots
	for (const auto& waveNameA : _metadata.waveNames) {
		for(size_t iWaveNameB = 0; iWaveNameB < std::min(nRefWaves,waveNamesOrderedByIntensity.size()); ++iWaveNameB){
			const std::string& waveNameB = waveNamesOrderedByIntensity[iWaveNameB];
			phaseSpectrum(waveNameA, waveNameB); // triggers generation of phase spectrum
		}
	}

	// build plots containing fit information
	negLogLikeSpectrum();
}


bool
rpwa::multibinPlots::buildMultibinSummedPlots(const std::vector<const multibinPlots*>& plotsInMultibins) {
	if (_initialized) {
		printErr<< "Cannot use initialized multibinPlots object to build multibin-summed plots" << std::endl;
		return false;
	}
	_initialized = true;
	if(plotsInMultibins.size() == 0) return false;

	_metadata.labels = plotsInMultibins[0]->labels();
	_metadata.descriptions = plotsInMultibins[0]->descriptions();
	for(const auto& plots: plotsInMultibins){
		if (plots->labels().size() != 1){
			printErr << "Cannot build multibin-summed plots from multibin plots with more than one solution" << std::endl;
			return false;
		}
		if (plots->labels()[0] != _metadata.labels[0]) {
			printErr << "Cannot build multibin-summed plots from multibin plots with different solution labels" << std::endl;
			return false;
		}
	}

	// build combined wave-name list and intensity-names list
	std::set<std::string> waveNames;
	std::set<std::string> intensityNames;
	for(const auto& plotsInMultibin: plotsInMultibins) {
		waveNames.insert(plotsInMultibin->waveNames().begin(), plotsInMultibin->waveNames().end());
		for(const auto& name_plot: plotsInMultibin->intensitySpectra()) intensityNames.insert(name_plot.first);
	}
	for(const auto& waveName: waveNames) _metadata.waveNames.push_back(waveName);

	// build summed intensity spectra
	for(const auto& intensityName: intensityNames) {
		std::vector<const componentPlot*> plots;
		plots.reserve(plotsInMultibins.size());
		for(const auto& plotsInMultibin: plotsInMultibins) { // build list with all plots from the multibins that have this intensity name
			if (plotsInMultibin->hasIntensitySpectrum(intensityName)) {
				plots.push_back(plotsInMultibin->intensitySpectra().at(intensityName));
			}
		}
		std::set<std::string> componentNames;
		for(const auto& plot: plots) { // build list of all component names
			const std::set<std::string> componentsInPlot = plot->getComponentNames();
			componentNames.insert(componentsInPlot.begin(), componentsInPlot.end());
		}
		for(const auto& plot: plots) { // check that all multibins have the same components
			if(plot->getComponentNames() != componentNames) {
				printErr << "Can not build multibin-summed spectra from multibin plots with different component sets" << std::endl;
				return false;
			}
		}

		componentPlot* summedPlot = static_cast<componentPlot*>(new componentPlot::baseType(intensityName.c_str(), intensityName.c_str()));
		for(const auto& component: componentNames) {
			const double* X = plots[0]->getComponent(component)->GetX();
			const double* XErr = plots[0]->getComponent(component)->GetEX();
			const int nPoints = plots[0]->getComponent(component)->GetN();
			componentPlot::plotType* summedComponentPlot = new componentPlot::plotType(nPoints);
			for(int iPoint = 0; iPoint < nPoints; ++iPoint) {
				const double x = X[iPoint];
				double y = 0.0;
				double yErrQrt = 0.0;
				for(const auto& plot: plots) {
					if(plot->getComponents(component).size() > 1) {
						printErr << "Cannot build multibin-summed spectra from multibin plots with more than one plot per component" << std::endl;
						return false;
					}

					const componentPlot::plotType* componentPlot = plot->getComponent(component);
					if(nPoints != componentPlot->GetN() or x != componentPlot->GetX()[iPoint]) {
						printErr << "Cannot build multibin-summed spectra from ultibin plots with different mass binning" << std::endl;
						return false;
					}

					y += componentPlot->GetY()[iPoint];
					yErrQrt += componentPlot->GetEY()[iPoint] * componentPlot->GetEY()[iPoint];
				}
				summedComponentPlot->SetPoint(iPoint, x, y);
				summedComponentPlot->SetPointError(iPoint, XErr[iPoint], std::sqrt(yErrQrt));
			}
			summedPlot->addForComponent(component, label(), summedComponentPlot);
		}
		_intensities[intensityName] = summedPlot;
	}

	// build dummy plots for some stuff
	_negLogLikeSpectrum = static_cast<componentPlot*>(new componentPlot::baseType("negLogLikeSpectrum", "neg. log(likelihood)"));

	return true;
}

bool rpwa::multibinPlots::mergePlotsInto(const multibinPlots& other){

	for( const auto& label: other.labels()){
		if(std::find(_metadata.labels.begin(), _metadata.labels.end(), label) != _metadata.labels.end()){
			printErr << "Label '" << label << "' appears in both plot collections" << std::endl;
			return false;
		}
	}
	// merge intensities
	for(const auto& name_plot: other.intensitySpectra()){
		const std::string& name = name_plot.first;
		const componentPlot* plot = name_plot.second;
		std::map<std::string, componentPlot*>::iterator it = _intensities.find(name);
		if (it == _intensities.end()){ // no intensity with this name, store a copy
			_intensities[name] = plot->copy();
		} else {
			_intensities[name]->merge(plot);
		}
	}

	// merge phases
	for (const auto& elemAB : other._phases) {
		const std::string& nameA = elemAB.first;
		std::map<std::string, std::map<std::string,componentPlot*>>::iterator itAB = _phases.find(nameA);
		for (const auto& elemB : elemAB.second) {
			const std::string& nameB = elemB.first;
			const componentPlot* plot = elemB.second;
			if (itAB == _phases.end()) { // no phases of waveA found
				_phases[nameA][nameB] = plot->copy();
			} else { // phases of waveA found
				std::map<std::string, componentPlot*>::iterator itB = _phases[nameA].find(nameB);
				if (itB == _phases[nameA].end()) { // no phases of waveA w.r.t. waveB found
					_phases[nameA][nameB] = plot->copy();
				} else {
					_phases[nameA][nameB]->merge(plot);
				}
			}
		}
	}

	// merge other plots
	for(const auto& name_plot: other._additionalPlots){
		const std::string& name = name_plot.first;
		const componentPlot* plot = name_plot.second;
		std::map<std::string, componentPlot*>::iterator it = _additionalPlots.find(name);
		if (it == _additionalPlots.end()){ // no plot with this name, store a copy
			_additionalPlots[name] = plot->copy();
		} else {
			_additionalPlots[name]->merge(plot);
		}
	}

	// merge additional fit information
	_negLogLikeSpectrum->merge(other.negLogLikeSpectrum());

	for(const auto& desc: other._metadata.descriptions)_metadata.descriptions.push_back(desc);
	for(const auto& label: other._metadata.labels)_metadata.labels.push_back(label);

	std::set<std::string> waveNameSet;
	waveNameSet.insert(_metadata.waveNames.begin(), _metadata.waveNames.end());
	waveNameSet.insert(other._metadata.waveNames.begin(), other._metadata.waveNames.end());
	_metadata.waveNames.clear();
	_metadata.waveNames.reserve(waveNameSet.size());
	for(const auto& waveName: waveNameSet) _metadata.waveNames.push_back(waveName);

	return true;
}

bool
rpwa::multibinPlots::write(TDirectory* directory) {
	// write intensities
	{
		TDirectory* directoryIntensities = directory->mkdir("intensities");
		directoryIntensities->cd();
		for (const auto& wave_stack : _intensities) {
			componentPlot* const & intensitySpectrum = wave_stack.second;
			static_cast<componentPlot::baseType*>(intensitySpectrum)->Write(intensitySpectrum->GetName());
		}
		directoryIntensities->Close();
	}

	// write phases
	{
		TDirectory* directoryPhases = directory->mkdir("phases");
		directoryPhases->cd();
		for (const auto& waveA_phases : _phases) {
			const std::string& waveNameA = waveA_phases.first;
			TDirectory* directoryWaveA = directoryPhases->mkdir(waveNameA.c_str());
			directoryWaveA->cd();
			for (const auto& waveB_phase : waveA_phases.second) {
				componentPlot* const & phaseSpectrum = waveB_phase.second;
				static_cast<componentPlot::baseType*>(phaseSpectrum)->Write(phaseSpectrum->GetName());
			}
			directoryWaveA->Close();
		}
		directoryPhases->Close();
	}

	// write additional plots
	{
		TDirectory* directoryAddPlots = directory->mkdir("additionalPlots");
		directoryAddPlots->cd();
		for (const auto& name_plot : _additionalPlots) {
			static_cast<componentPlot::baseType*>(name_plot.second)->Write(name_plot.first.c_str());
		}
	}

	// write additonal fit information
	{
		TDirectory* directoryAddFitInfo = directory->mkdir("additionalFitInformation");
		directoryAddFitInfo->cd();
		static_cast<componentPlot::baseType*>(_negLogLikeSpectrum)->Write("negLogLikeSpectrum");
	}

	// write fit results
	{
		directory->cd();
		std::string treeName = "fitResults_";
		treeName += directory->GetName();
		TTree* fitResultTree = new TTree(treeName.c_str(), "fitResults");
		std::vector<fitResult> resultsOut;
		fitResultTree->Branch("results", &resultsOut);
		for (const auto& massbin_result : _fitResultsInMassbins) {
			resultsOut.clear();
			resultsOut = massbin_result.second;
			fitResultTree->Fill();
		}
		fitResultTree->Write();
	}

	// write metadata
	{
		directory->cd();
		_metadata.Write("metaData");
	}

	directory->Close();
	return true;
}


bool
rpwa::multibinPlots::load(TDirectory* directory, const bool onlyBest) {
	if (_initialized) {
		printErr<< "Cannot use initialized multibinPlots object to load multibin plots" << std::endl;
		return false;
	}
	_initialized = true;
	// load intensities
	{
		TDirectory* directoryIntensities = dynamic_cast<TDirectory*>(directory->Get("intensities"));
		if (directoryIntensities != nullptr) {
			TList* keys = directoryIntensities->GetListOfKeys();
			for (int i = 0; i < keys->GetEntries(); ++i) {
				componentPlot::baseType* plots = dynamic_cast<componentPlot::baseType*>(dynamic_cast<TKey*>(keys->At(i))->ReadObj());
				if( plots != nullptr) {
					_intensities[plots->GetName()] = static_cast<componentPlot*>(plots);
				} else {
					printErr<< "Could not load intensity plot for wave '" << keys->At(i)->GetName() << "' in multiplots folder '"
					<< directory->GetName() << "'." << std::endl;
				}
			}
		} else {
			printErr<< "Could not find intensity folder in multibinplots folder '"
			<< directory->GetName() << "'." << std::endl;
			return false;
		}
	}

	// load phases
	{
		TDirectory* directoryPhases = dynamic_cast<TDirectory*>(directory->Get("phases"));
		if (directoryPhases != nullptr) {
			TList* keysA = directoryPhases->GetListOfKeys();
			for (int i = 0; i < keysA->GetEntries(); ++i) {
				const std::string waveNameA = keysA->At(i)->GetName();
				TDirectory* directoryWaveA = dynamic_cast<TDirectory*>(dynamic_cast<TKey*>(keysA->At(i))->ReadObj());
				if(directoryWaveA != nullptr) {
					TList* keysB = directoryWaveA->GetListOfKeys();
					for(int j = 0; j < keysB->GetEntries(); ++j) {
						const std::string waveNameAB = keysB->At(j)->GetName(); // <waveNameA>__<waveNameB>
						const std::string waveNameB = waveNameAB.substr(waveNameA.size() + 2);
						componentPlot::baseType* plots = dynamic_cast<componentPlot::baseType*>(dynamic_cast<TKey*>(keysB->At(j))->ReadObj());
						if( plots != nullptr) {
							_phases[waveNameA][waveNameB] = static_cast<componentPlot*>(plots);
						} else {
							printErr<< "Could not laod phase plot for waves '" << waveNameA << "' and '" << waveNameB << "' in multiplots folder '"
							<< directory->GetName() << "'." << std::endl;
						}
					}
				} else {
					printErr<< "Found object other than directory in phase folder for wave '" << waveNameA << "' in multibinplots folder '"
					<< directory->GetName() << "'." << std::endl;
				}
			}
		} else {
			printErr<< "Could not find phase folder in multibinplots folder '"
			<< directory->GetName() << "'." << std::endl;
			return false;
		}
	}

	// load additional plots
	{
		TDirectory* directoryAddPlots = dynamic_cast<TDirectory*>(directory->Get("additionalPlots"));
		if (directoryAddPlots!= nullptr) {
			TList* keys = directoryAddPlots->GetListOfKeys();
			for (int i = 0; i < keys->GetEntries(); ++i) {
				componentPlot::baseType* plots = dynamic_cast<componentPlot::baseType*>(dynamic_cast<TKey*>(keys->At(i))->ReadObj());
				if( plots != nullptr) {
					_additionalPlots[plots->GetName()] = static_cast<componentPlot*>(plots);
				} else {
					printErr<< "Could not load additional plot '" << keys->At(i)->GetName() << "' in multiplots folder '"
					<< directory->GetName() << "'." << std::endl;
				}
			}
		} else {
			printErr<< "Could not find additional plots folder in multibinplots folder '"
			<< directory->GetName() << "'." << std::endl;
			return false;
		}
	}


	// load fit results
	{
		std::string treeName = "fitResults_";
		treeName += directory->GetName();
		TTree* fitResultTree = dynamic_cast<TTree*>(directory->Get(treeName.c_str()));
		if(fitResultTree != nullptr) {
			std::vector<fitResult> resultsIn;
			std::vector<fitResult>* resultsInPtr = &resultsIn;
			fitResultTree->SetBranchAddress("results", &resultsInPtr);
			for(int i = 0; i < fitResultTree->GetEntries(); ++i) {
				fitResultTree->GetEntry(i);
				if(resultsIn.size() > 0 ) {
					const double massBinCenter = resultsIn[0].massBinCenter();
					if(onlyBest) {
						_fitResultsInMassbins[massBinCenter].push_back(resultsIn[0]);
					} else {
						_fitResultsInMassbins[massBinCenter] = resultsIn;
					}
				} else {
					printErr << "Massbin '" << i << "' in fitresult tree is empty in multibinplots folder '"
					<< directory->GetName() << "'." << std::endl;
					return false;
				}
			}
		} else {
			printErr<< "Could not find fitresults tree in multibinplots folder '"
			<< directory->GetName() << "'." << std::endl;
			return false;
		}
	}

	// load additional fit information
	{
		TDirectory* directoryAddFitInfo = dynamic_cast<TDirectory*>(directory->Get("additionalFitInformation"));
		if (directoryAddFitInfo!= nullptr) {
			TKey* keyNegLogLike = directoryAddFitInfo->FindKey("negLogLikeSpectrum");
			if (keyNegLogLike != nullptr){
				componentPlot::baseType* p = dynamic_cast<componentPlot::baseType*>(keyNegLogLike->ReadObj());
				if (p != nullptr){
					_negLogLikeSpectrum = static_cast<componentPlot*>(keyNegLogLike->ReadObj());
				} else {
					printErr << "Cannot read neg. log-likelihood plot in multibinplots folder '"
					<< directory->GetName() << "'." << std::endl;
					return false;
				}
			}else{
				printWarn << "Cannot find neg. log-likelihood plot in multibinplots folder '"
					<< directory->GetName() << "'. Try to build from fit results." << std::endl;
				if(negLogLikeSpectrum() == nullptr) return false;
			}
		} else {
			printErr<< "Could not find additional fit information folder in multibinplots folder '"
			<< directory->GetName() << "'." << std::endl;
			return false;
		}
	}

	// load meta data
	{
		multibinplotsMetadata* metaData = dynamic_cast<multibinplotsMetadata*>(directory->Get("metaData"));
		if(metaData != nullptr) {
			_metadata = *metaData;
		} else {
			printErr<< "Could not find meta data in multibinplots folder '"
			<< directory->GetName() << "'." << std::endl;
			return false;
		}
	}

	return true;
}


double
rpwa::multibinPlots::calcIntensityIntegralRegEx(const std::string& waveNamePattern, const double xmin, const double xmax) {
	double totIntensity = 0.0;
	if (waveNamePattern.size() > 0) {
		componentPlot* plots = intensitySpectrumRegEx(waveNamePattern);
		if (plots != nullptr) {
			rpwa::componentPlot::plotType* plot = plots->getMassindepComponents()[_metadata.labels[0]];
			if (plot != nullptr) {
				totIntensity = 0.0;
				for (int i = 0; i < plot->GetN(); ++i) {
					const double xLower = plot->GetX()[i] - plot->GetEX()[i];
					const double xUpper = plot->GetX()[i] + plot->GetEX()[i];
					if ((std::isnan(xmin) or xUpper > xmin) and (std::isnan(xmax) or xLower < xmax)) {
						totIntensity += plot->GetY()[i];
					}
				}
			} else {
				printErr<< "Can not calculate total intensity for wave '" << waveNamePattern
				<< "'. Cannot find plot in multiplots!" << std::endl;
			}
		} else {
			printErr<< "Can not calculate total intensity for wave '" << waveNamePattern
			<< "' (cannot find waveName '" << waveNamePattern << "' )!" << std::endl;
		}
	} else { // sum of all waves
		totIntensity = 0.0;
		for (const auto& wave : waveNames()) {
			totIntensity += calcIntensityIntegral(wave, xmin, xmax);
		}
	}
	return totIntensity;
}

componentPlot* rpwa::multibinPlots::negLogLikeSpectrum()
{
	if (_negLogLikeSpectrum == nullptr){
		if (_fitResultsInMassbins.size() > 0) { // we wave loaded fit results
			componentPlot::plotType* p = nullptr;
			p = new componentPlot::plotType(_fitResultsInMassbins.size());
			int i = 0;
			for (const auto& mass_results : _fitResultsInMassbins) {
				const fitResult& result = mass_results.second[0];
				const double negLogLike = result.logLikelihood();
				const double x = result.multibinCenter().at("mass");
				p->SetPoint(i, x, negLogLike);
				i++;
			}
			_negLogLikeSpectrum = static_cast<componentPlot*>(new componentPlot::baseType("negLogLikeSpectrum", "-log(likelihood)"));
			_negLogLikeSpectrum->addForComponent(plotComponents::massIndependent, label(), p);
		} else {
			printErr<< "Can not generate likelihood stack because fit results are not stored!" << std::endl;
		}
	}
	return _negLogLikeSpectrum;
}

componentPlot*
rpwa::multibinPlots::_intensitySpectrum(const std::string& waveNamePattern) {
	if (_intensities.find(waveNamePattern) != _intensities.end()) {
		return _intensities.at(waveNamePattern);
	} else {
		componentPlot* plots = nullptr;
		componentPlot::plotType* p = genIntensitySpectrumRegEx(waveNamePattern);
		if (p != nullptr) { // if we can generate the intensity spectrum
			plots = static_cast<componentPlot*>(new componentPlot::baseType(waveNamePattern.c_str(), waveNamePattern.c_str()));
			plots->addForComponent(plotComponents::massIndependent, label(), p);
			_intensities[waveNamePattern] = plots;
		}
		return plots;
	}
}


componentPlot*
rpwa::multibinPlots::phaseSpectrum(const std::string& waveNameA,
                                                  const std::string& waveNameB) {
	std::map<std::string, std::map<std::string, componentPlot*>>::const_iterator itA = _phases.find(waveNameA);
	std::map<std::string, componentPlot*>::const_iterator itB;
	if (itA != _phases.end() and (itB = itA->second.find(waveNameB)) != itA->second.end()) { // wave already exists
		return itB->second;
	} else {
		if(hasMultipleSolutions()){
			printErr << "Cannot build phase spectra for multibin plots with more than one solution" << std::endl;
			return nullptr;
		}
		componentPlot* plots = nullptr;
		componentPlot::plotType* p = genPhaseSpectrum(waveNameA, waveNameB);
		if (p != nullptr) { // if we can generate the phase spectrum
			const std::string spectrumName = waveNameA + "__" + waveNameB;
			plots = static_cast<componentPlot*>(new componentPlot::baseType(spectrumName.c_str(), spectrumName.c_str()));
			plots->addForComponent(plotComponents::massIndependent, label(), p);
			_phases[waveNameA][waveNameB] = plots;
		}
		return plots;
	}
}


componentPlot::plotType*
rpwa::multibinPlots::genIntensitySpectrumRegEx(const std::string& waveNamePattern) {
	componentPlot::plotType* p = nullptr;
	if (_fitResultsInMassbins.size() > 0) { // we wave loaded fit results
		p = new componentPlot::plotType(_fitResultsInMassbins.size());
		int i = 0;
		for (const auto& mass_results : _fitResultsInMassbins) {
			const fitResult& result = mass_results.second[0];
			const double I = result.intensity(waveNamePattern);
			const double Ierr = result.hasHessian() ? result.intensityErr(waveNamePattern) : 0.0;
			const double x = result.multibinCenter().at("mass");
			const double xerr = 0.5 * (result.multibinBoundaries().at("mass").second - result.multibinBoundaries().at("mass").first);
			p->SetPoint(i, x, I);
			p->SetPointError(i, xerr, Ierr);
			i++;
		}
	} else {
		printErr<< "Can not generate intensity stack because fit results are not stored!" << std::endl;
	}
	return p;
}


rpwa::componentPlot::plotType*
rpwa::multibinPlots::genPhaseSpectrum(const std::string& waveNameA, const std::string& waveNameB) {
	componentPlot::plotType* p = nullptr;
	if (_fitResultsInMassbins.size() > 0) { // we wave loaded fit results
		p = new componentPlot::plotType(_fitResultsInMassbins.size());
		int i = 0;
		for (const auto& mass_results : _fitResultsInMassbins) {
			const fitResult& result = mass_results.second[0];
			const double phase = result.phase(waveNameA, waveNameB);
			const double phaseerr = result.hasHessian() ? result.phaseErr(waveNameA, waveNameB) : 0.0;
			const double x = result.multibinCenter().at("mass");
			const double xerr = 0.5 * (result.multibinBoundaries().at("mass").second - result.multibinBoundaries().at("mass").first);
			p->SetPoint(i, x, phase);
			p->SetPointError(i, xerr, phaseerr);
			i++;
		}
	} else {
		printErr<< "Can not generate phase stack because fit results not stored!" << std::endl;
	}
	return p;
}


bool
rpwa::multibinPlots::addAdditionalPlot(componentPlot::plotType* plot, const std::string& name) {
	componentPlot* cPlot = new componentPlot();
	cPlot->addForComponent(plot->GetName(), label(), plot);
	return addAdditionalPlot(cPlot, name);
}


bool
rpwa::multibinPlots::addAdditionalPlot(componentPlot* plot, const std::string& name) {
	if (_additionalPlots.find(name) == _additionalPlots.end()) {
		plot->SetName(name.c_str());
		_additionalPlots[name] = plot;
		return true;
	}
	return false;
}


componentPlot*
rpwa::multibinPlots::getAdditionalPlot(const std::string& name) {
	std::map<std::string, componentPlot*>::iterator it = _additionalPlots.find(name);
	if (it != _additionalPlots.end()) {
		return it->second;
	}
	return nullptr;
}


std::vector<std::string>
rpwa::multibinPlots::waveNamesSortedByIntensity()
{
	std::vector<std::string> waveNames = this->waveNames();
	// build map of intensities
	std::map<std::string, double> intensitys;
	for(const auto& waveName: waveNames) intensitys[waveName] = calcIntensityIntegral(waveName);
	std::sort(waveNames.begin(), waveNames.end(),
	          [&intensitys](const std::string& a, const std::string& b) -> bool {return intensitys[a] > intensitys[b];});
	return waveNames;
}


std::map<std::string, rpwa::componentPlot::plotType*>
rpwa::componentPlot::getComponents(const std::string& componentName) const {
	std::map<std::string, rpwa::componentPlot::plotType*> histsOfComponentName;
	TList* allPlots = GetListOfGraphs();
	for (int i = 0; i < allPlots->GetEntries(); ++i) {
		rpwa::componentPlot::plotType* h = static_cast<rpwa::componentPlot::plotType*>(allPlots->At(i));
		if (std::string(h->GetName()).find(componentName) != std::string::npos) {
			// the graph belongs to the requested component name
			histsOfComponentName[std::string(h->GetTitle())] = h;
		}
	}
	return histsOfComponentName;
}


std::map<std::string, rpwa::componentPlot::plotType*>
rpwa::componentPlot::getComponents(const int componentType) const {
	return getComponents(plotComponents::names.at(componentType));
}


void
rpwa::componentPlot::addForComponent(const std::string& componentName, const std::string& label, rpwa::componentPlot::plotType* p, Option_t* chopt) {
	p->SetName(componentName.c_str());
	p->SetTitle(label.c_str());
	TMultiGraph::Add(p, chopt);
}


rpwa::componentPlot::plotType*
rpwa::componentPlot::getComponent(const std::string& componentName) const {
	rpwa::componentPlot::plotType* hist = nullptr;
	const std::map<std::string, rpwa::componentPlot::plotType*> all = getComponents(componentName);
	if (all.size() == 1) {
		hist = all.begin()->second;
	} else {
		printErr<< "componentPlot has more than one plot per component" << std::endl;
	}
	return hist;
}


rpwa::componentPlot::plotType*
rpwa::componentPlot::getComponent(const int componentType) const {
	return getComponent(plotComponents::names.at(componentType));
}


std::set<std::string>
rpwa::componentPlot::getComponentNames() const {
	std::set<std::string> names;
	for (int i = 0; i < GetListOfGraphs()->GetEntries(); ++i) {
		names.insert(GetListOfGraphs()->At(i)->GetName());
	}
	return names;
}


void
rpwa::componentPlot::addForComponent(const int componentType, const std::string& label, rpwa::componentPlot::plotType* p, Option_t* chopt) {
	componentPlot::addForComponent(plotComponents::names.at(componentType), label, p, chopt);
}


void
rpwa::componentPlot::merge(const componentPlot* other){
	TList* otherGraphs = other->GetListOfGraphs();
	if (otherGraphs != nullptr){ // if the other componentPlot has graphs
		for( int i = 0; i < otherGraphs->GetEntries(); ++i){
			plotType* g = dynamic_cast<plotType*>(otherGraphs->At(i));
			if (g == nullptr){
				printErr << "Cannot merge componentPlots" << std::endl;
				throw std::invalid_argument("Cannot merge componentPlots");
			}
			plotType* gNew = dynamic_cast<plotType*>(g->Clone());
			gNew->SetName(g->GetName());
			gNew->SetTitle(g->GetTitle());
			Add(gNew);
		}
	}
}


componentPlot*
rpwa::componentPlot::copy() const {
	componentPlot* newPlot = static_cast<componentPlot*>(new componentPlot::baseType());
	newPlot->SetName(this->GetName());
	newPlot->SetTitle(this->GetTitle());
	newPlot->merge(this); // generates a copy of each plot in this object
	return newPlot;
}


void
rpwa::plottingtools::shiftPhaseSpectrumInRange(componentPlot::plotType* plot, const double rangeStart) {
	const double rangeEnd = rangeStart + 360.0;

	for (int i = 0; i < plot->GetN(); ++i) {
		double y = plot->GetY()[i];
		while (y < rangeStart)
			y += 360.0;
		while (y >= rangeEnd)
			y -= 360.0;
		plot->SetPoint(i, plot->GetX()[i], y);
	}
}


double
rpwa::plottingtools::makePhaseContinousWithinRange(componentPlot::plotType* plot, const int nTrails) {

	double bestRangeStart = -180.0;
	double bestRangeDiscontinuity = std::numeric_limits<double>::max();
	for (int iTrail = 0; iTrail < nTrails; ++iTrail) {
		const double rangeStart = -220.0 + 260.0 / (nTrails - 1) * iTrail;
		plottingtools::shiftPhaseSpectrumInRange(plot, rangeStart);
		const double discontinuity = phaseSpectrumDiscontinuity(plot);
		if (discontinuity != 0 and discontinuity < bestRangeDiscontinuity) {
			bestRangeDiscontinuity = discontinuity;
			bestRangeStart = rangeStart;
		}
	}
	plottingtools::shiftPhaseSpectrumInRange(plot, bestRangeStart);
	return bestRangeStart;
}


namespace {
	double
	phaseSpectrumDiscontinuity(componentPlot::plotType* plot) {
		double discontinuity = 0.0;
		const double* y = plot->GetY();
		const double* yErr = plot->GetEY();
		for (int i = 0; i < plot->GetN() - 1; ++i) {
			if (y[i] != 0 and y[i + 1] != 0)
				if (yErr[i] != 0 and yErr[i + 1] != 0)
					discontinuity += (y[i] - y[i + 1]) * (y[i] - y[i + 1]) / yErr[i] / yErr[i + 1];
		}
		return discontinuity;
	}


	std::set<std::string>
	buildAllTotals(const std::vector<std::string>& waveNames) {

		// create set of all totals
		std::set<std::string> allTotals;
		allTotals.insert(".*");

		std::string IG;
		bool firstWave(true);
		bool oldNaming(true);
		for (const auto& waveName : waveNames) {
			if (waveName != "flat") {
				if (firstWave) {
					if (waveName[0] == '[')
						oldNaming = false;
				} else {
					if (not oldNaming)
						assert(waveName[0] == '[');
					else
						assert(waveName[0] != '[');
				}

				// we can only handle fits were the quantum numbers IG are the same for all waves
				const size_t length = oldNaming ? 2 : 4;
				if (firstWave) {
					IG = waveName.substr(0, length);
				} else {
					assert(IG == waveName.substr(0, length));
				}

				firstWave = false;
			}

			const std::set<std::string> newTotals = extractTotalPatternsFromWavename(waveName, oldNaming);
			allTotals.insert(newTotals.begin(), newTotals.end());

		}
		return allTotals;
	}


	std::set<std::string>
	extractTotalPatternsFromWavename(std::string name, const bool oldNaming) {
		std::set<std::string> patterns;

		if (name == "flat") {
			patterns.insert("flat");
			// everything except flat wave for plot showing total intensity minus
			// intensity of the flat wave
			patterns.insert("(?!flat$).*");
			return patterns;
		}

		name = rpwa::escapeRegExpSpecialChar(name);

		if (not oldNaming) {
			// new naming convention "[IG,JPC,ME]=..."
			// * the C might not be present
			// assuming: * '[' and ']' have been escaped to '\[' and '\]'

			// from [IG,JPC,MR]=[...] create
			const size_t splitQn1 = name.find(',');
			const size_t posP = std::min(name.find('+', splitQn1 + 1), name.find('-', splitQn1 + 1));
			const size_t splitQn2 = name.find(',', posP + 1);
			const size_t splitDecay = name.find('=');
			std::string dots = "";
			for (size_t i = posP; i < splitQn2; ++i)
				dots += ".";
			const std::string jpcPattern = "\\d+" + dots;
			// 1.) IG JPC MR
			patterns.insert(name.substr(0, splitDecay) + "=\\[.*\\]");
			// 2.) IG JPC M
			patterns.insert(name.substr(0, splitQn2 + 2) + "." + name.substr(splitQn2 + 4, splitDecay - splitQn2 - 4) + "=\\[.*\\]");
			// 3.) IG JPC R
			patterns.insert(name.substr(0, splitQn2 + 1) + "\\d+" + name.substr(splitQn2 + 2, splitDecay - splitQn2 - 2) + "=\\[.*\\]");
			// 4.) IG JPC
			patterns.insert(name.substr(0, splitQn2 + 1) + "\\d+." + name.substr(splitDecay - 2, 2) + "=\\[.*\\]");
			// 5.) IG MR
			patterns.insert(name.substr(0, splitQn1 + 1) + jpcPattern + name.substr(splitQn2, splitDecay - splitQn2) + "=\\[.*\\]");
			// 6.) IG M
			patterns.insert(name.substr(0, splitQn1 + 1) + jpcPattern + name.substr(splitQn2, 2) + "." + name.substr(splitDecay - 2, 2) + "=\\[.*\\]");
			// 7.) IG R
			patterns.insert(
			        name.substr(0, splitQn1 + 1) + jpcPattern + name.substr(splitQn2, 1) + "\\d+" + name.substr(splitQn2 + 2, splitDecay - splitQn2 - 2)
			                + "=\\[.*\\]");

		} else {
			name = rpwa::unescapeRegExpSpecialChar(name);

			// old naming convention "IGJPCME..."
			// from IG JPC MR create
			// 1.) IG JPC MR
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 7)) + ".*");
			// 2.) IG JPC M
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 6)) + "." + ".*");
			// 3.) IG JPC R
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 5)) + "." + rpwa::escapeRegExpSpecialChar(name.substr(6, 1)) + ".*");
			// 4.) IG JPC
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 5)) + ".." + ".*");
			// 5.) IG MR
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 2)) + "..." + rpwa::escapeRegExpSpecialChar(name.substr(5, 2)) + ".*");
			// 6.) IG M
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 2)) + "..." + rpwa::escapeRegExpSpecialChar(name.substr(5, 1)) + "." + ".*");
			// 7.) IG R
			patterns.insert(rpwa::escapeRegExpSpecialChar(name.substr(0, 2)) + "...." + rpwa::escapeRegExpSpecialChar(name.substr(6, 1)) + ".*");
		}
		return patterns;
	}
}

