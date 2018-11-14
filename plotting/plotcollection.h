#ifndef PLOTCOLLECTION_H
#define PLOTCOLLECTION_H

#include <map>
#include <string>
#include <vector>
#include <set>

#include "TObject.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH2D.h"
#include "fitResult.h"

namespace rpwa {

	class plotComponents {
	public:

		enum IDs{
			massIndependent
		};
		const static std::vector<std::string> names;

		plotComponents() {}

		~plotComponents() {}

	};

	class componentPlot: public TMultiGraph{
	public:
		typedef TGraphErrors plotType;
		typedef TMultiGraph baseType;
		componentPlot(){}
		componentPlot(const std::string& name):
			TMultiGraph(name.c_str(), name.c_str()){}
		componentPlot(const std::string& name, const std::string& title):
			TMultiGraph(name.c_str(), title.c_str()){}

		virtual ~componentPlot(){}

		/***
		 * @param componentName component name according to plotcomponents::componentNames
		 * @return label, plot for all plots that correspond to the given component name
		 */
		std::map<std::string, plotType*> getComponents(const std::string& componentName)const;

		/***
		 * @param componentType component type according to plotcomponents::components
		 * @return label, plot for all plots that correspond to the given component name
		 */
		std::map<std::string, plotType*> getComponents(const int componentType)const;

		/***
		 * @return label, plot for all plots that correspond to the mass-independent component
		 */
		std::map<std::string, plotType*> getMassindepComponents()const{return getComponents(plotComponents::massIndependent);}

		/***
		 * @param componentName component name according to plotcomponents::componentNames
		 * @return the first plot that correspond to the given component name (nulltpr if not found)
		 */
		plotType* getComponent(const std::string& componentName)const;

		/***
		 * @param componentType component type according to plotcomponents::components
		 * @return the first plot that correspond to the given component name (nulltpr if not found)
		 */
		plotType* getComponent(const int componentType)const;

		/***
		 * @return the first plot that correspond to the mass-independent component (nulltpr if not found)
		 */
		plotType* getMassindepComponent()const{return getComponent(plotComponents::massIndependent);}

		/***
		 * @return set of all component names present in the componentPlot object
		 */
		std::set<std::string> getComponentNames()const;

		/***
		 * @param componentName component name of the component the graph represents
		 * @param p plot to add to the multigraph
		 * @param chopt Options
		 */
		void addForComponent(const std::string& componentName, const std::string& label, plotType* p, Option_t* chopt = "");

		/***
		 * @param componentType component type according to plotcomponents::components
		 * @param p plot to add to the multigraph
		 * @param chopt Options
		 */
		void addForComponent(const int componentType, const std::string& label, plotType* p, Option_t* chopt = "");

		/***
		 * Merge the plots of the other component plot into this one.
		 *
		 * A copy of each plot is stored in this object
		 */
		void merge(const componentPlot* other);

		/***
		 * Generate a full copy, including copy of the individual plots, of this componentPlot
		 * @return pointer of type 'componentPlot*' to the new object of type 'componentPlot::baseType'
		 */
		componentPlot* copy() const;

		ClassDef(componentPlot, 1)
	};


	/***
	 *  Class to handle meta data for multibinplots
	 */
	class multibinplotsMetadata: public TObject {
	public:
		std::vector<std::string> descriptions;
		std::vector<std::string> labels;
		std::vector<std::string> waveNames;

	ClassDef(multibinplotsMetadata, 1);
	};

	class multibinPlots {
	public:
		multibinPlots();
		/***
		 *
		 * @param fitresults List of fit results (can have multiple results per mass bin)
		 * @param soltionDescription Description for the soltion the fit results belong to
		 */
		multibinPlots(const std::vector<rpwa::fitResult>& fitresults, const std::string& label, const std::string& description);

		~multibinPlots() {}

		/***
		 * Build all intensity spectra and phase plots.
		 * Requires fit results to be loaded.
		 * Also builds spin-totals of I^G J^PC M^e and subsets
		 * @param waveNamePatterns additional wave-name patterns for which generate intensity spectra will be generated
		 * @param nRefWaves Rel. phases are builded for all waves w.r.t. the nRefWaves largest waves
		 */
		void buildDefaultPlots(const std::vector<std::string> waveNamePatterns = {}, const size_t nRefWaves = 4);

		/***
		 * Build the multibin-summed spectra from the list of plots in the individual multibins.
		 * The following summed plots will be generated:
		 *     - incoh. summed intensity spectra of waves and all spin-totals in all multi-bins
		 * @param plotsInMultibins list of plots in the individual multibins
		 * @return true if buldint the plots was successful
		 */
		bool buildMultibinSummedPlots(const std::vector<const multibinPlots*>& plotsInMultibins);

		/***
		 * Merge the plots of the other multibinPlots object into this one
		 * Fit resutls will not be merged -> only the fit results of the first label stay in the object
		 * @return true if merging was successful
		 */
		bool mergePlotsInto(const multibinPlots& other);

		/***
		 * Store fit results and plots to the given root directory
		 * @return true if write was successful
		 */
		bool write(TDirectory* directory);

		/***
		 * Load fit results and plots from the given root directory
		 * @param onlyBest if true, load only the best fit result ber mass bin
		 * @return true if load was successful
		 */
		bool load(TDirectory* directory, const bool onlyBest=false);

		/**
		 * Override the internally stored intensities
		 * @param intensities Intensities, where the key is the wave name
		 */
		void setIntensitySpectra(const std::map<std::string, componentPlot*>& intensities) {_intensities = intensities;}

		/**
		 * Override the internally stored phases
		 * @param phases Phase spectra, where the key is the wave name
		 */
		void setPhaseSpectra(const std::map<std::string, std::map<std::string,componentPlot*>>& phases) {_phases = phases;}

		/**
		 * Set the (first) label to the given value
		 */
		void setLabel(const std::string& label){ if (_metadata.labels.size() > 0) _metadata.labels[0] = label; else _metadata.labels.push_back(label);}

		/**
		 * Set the (first) description to the given value
		 */
		void setDescription(const std::string& desc){ if (_metadata.descriptions.size() > 0) _metadata.descriptions[0] = desc; else _metadata.descriptions.push_back(desc);}

		/***
		 * Get intensity spectrum as multigraph.
		 * If the intensity spectrum does not jet exist, it generates its.
		 * @param waveNamePattern Any wave-name pattern (sam as in fitResult intensity)
		 * @return
		 */
		componentPlot* intensitySpectrumRegEx(const std::string& waveNamePattern){return _intensitySpectrum(waveNamePattern);}

		/***
		 * Get intensity spectrum as multigraph.
		 * If the intensity spectrum does not jet exist, it generates its.
		 * @param waveName any wave name (not regex pattern)
		 * @return
		 */
		componentPlot* intensitySpectrum(const std::string& waveName){return _intensitySpectrum(rpwa::escapeRegExpSpecialChar(waveName));}

		const std::map<std::string, componentPlot*>& intensitySpectra() const {return _intensities;}

		/***
		 * Checks wether the multibinPlots object has an intensity spectrum with the given name
		 * @param name wave-name for waves or wave-name pattern for spin-totals
		 * @return
		 */
		bool hasIntensitySpectrum(const std::string& name)const {return _intensities.find(name) != _intensities.end();}

		/***
		 * Get phase spectrum as multigraph.
		 * If the intensity spectrum does not jet exist, it generates its.
		 * @param waveNameA any wave name (not regex pattern)
		 * @param waveNameB any wave name (not regex pattern)
		 * @return
		 */
		componentPlot* phaseSpectrum(const std::string& waveNameA, const std::string& waveNameB);

		/***
		 * Generate an intensity spectrum
		 * in each bin
		 * @param waveNamePattern Any wave-name pattern (same as in fitResults intensity)
		 * @return Intensity spectrum of the given wave-name pattern
		 */
		rpwa::componentPlot::plotType* genIntensitySpectrumRegEx(const std::string& waveNamePattern);

		/***
		 * Generate an intensity spectrum
		 * in each bin
		 * @param waveName any wave name (not a regex pattern)
		 * @return Intensity spectrum of the given wave-name pattern
		 */
		rpwa::componentPlot::plotType* genIntensitySpectrum(const std::string& waveName){return genIntensitySpectrumRegEx(rpwa::escapeRegExpSpecialChar(waveName));}

		/***
		 * Generate a phase spectrum for [waveNameA] - [waveNameB]
		 * in each bin
		 * @param waveNameA Phase of this wave
		 * @param waveNameB relative to phase of this wave
		 * @return Phase spectrum of the given waves
		 */
		rpwa::componentPlot::plotType* genPhaseSpectrum(const std::string& waveNameA, const std::string& waveNameB);

		/***
		 * Add an additional plot to this collection
		 * @return true if the plot was added, false if the name already exists
		 */
		bool addAdditionalPlot(componentPlot* plot, const std::string& name);

		/***
		 * Add an additional plot to this collection
		 * @return true if the plot was added, false if the name already exists
		 */
		bool addAdditionalPlot(componentPlot::plotType* plot, const std::string& name);

		/***
		 * Get additional plot
		 * @return return nullptr if the plot does not exist
		 */
		componentPlot* getAdditionalPlot(const std::string& name);

		/***
		 * Return
		 * @return map of mass-bin centers and fit results per mass bin ordered convergence and neg-log-likelihood
		 */
		std::map<double, std::vector<fitResult> >& fitResultsInMassbins(){return _fitResultsInMassbins;}
		const std::map<double, std::vector<fitResult> >& fitResultsInMassbins() const {return _fitResultsInMassbins;}

		std::vector<double> massBinCenters() const;
		std::vector<boundaryType> massBinBoundaries() const;

		/***
		 * @return ordered list of wave names
		 */
		const std::vector<std::string>& waveNames() const {return _metadata.waveNames;}

		/**
		 * @return list of wave names ordered by intensity (largest intensity first)
		 */
		std::vector<std::string> waveNamesSortedByIntensity();

		/***
		 * @return descriptions of the multibinPlots object
		 */
		const std::vector<std::string>& descriptions() const {return _metadata.descriptions;}

		/***
		 * @return first description of the multibinPlots object
		 */
		const std::string& description() const {return _metadata.descriptions.at(0);}

		/***
		 * @return labels of the multibinPlots object
		 */
		const std::vector<std::string>& labels() const {return _metadata.labels;}

		/***
		 * @return first label of the multibinPlots object
		 */
		const std::string& label() const {return _metadata.labels.at(0);}

		/***
		 * @return true if the multibinPlots object has more than one solution
		 */
		bool hasMultipleSolutions() const {return _metadata.labels.size() > 1;}

		/***
		 * @return total intensity for the given wave name pattern.
		 *         If no name given, calculate the total intensity of all waves.
		 *         If xmin or xmax given, calculate the intensity within this range.
		 *         A x-bin is included in the total intensity if a part of it is within [xmin, xmax]
		 *         (upper bound > xmin && lower bound < xmax). Thus, the integrated range is >= [xmin, xmax].
		 *         If multiple solutions are stored, the first one is used to calculated the integral.
		 */
		double calcIntensityIntegralRegEx(const std::string& waveNamePattern, const double xmin = nan(""), const double xmax = nan(""));

		/***
		 * @return total intensity for the given wave name.
		 *         If no name given, calculate the total intensity of all waves.
		 *         If xmin or xmax given, calculate the intensity within this range.
		 *         A x-bin is included in the total intensity if a part of it is within [xmin, xmax]
		 *         (upper bound > xmin && lower bound < xmax). Thus, the integrated range is >= [xmin, xmax].
		 *         If multiple solutions are stored, the first one is used to calculated the integral.
		 */
		double calcIntensityIntegral(const std::string& waveName, const double xmin = nan(""), const double xmax = nan("")) {
			return calcIntensityIntegralRegEx(rpwa::escapeRegExpSpecialChar(waveName), xmin, xmax);
		}

		componentPlot* negLogLikeSpectrum() const {return _negLogLikeSpectrum;}
		componentPlot* negLogLikeSpectrum();
		TH2D* negLogLikeDistribution(const int nBinsNegLogLike = 100);
		/***
		 * @param massBinCenter mass bin center
		 * @return List of negative log likelihoods of all fit attempts in the given mass bin
		 */
		std::vector<double> negLogLikelihoods(const double massBinCenter) const;



	private:
		bool _initialized;
		std::map<std::string, componentPlot*> _intensities;
		// indices are [<waveA>][<waveB>}
		// shows phi(waveA) - phi(waveB)
		std::map<std::string, std::map<std::string, componentPlot*> > _phases;
		std::map<std::string, componentPlot*> _additionalPlots;
		// indices are [<massbin-center>][<fitattempt-id>], where
		// <fitattempt-id> is the id of the fit attempt
		// the fit attempts are ordered by first convergence and then by neg. log-likelihood
		std::map< double, std::vector< fitResult > > _fitResultsInMassbins;
		rpwa::multibinplotsMetadata _metadata;
		componentPlot* _negLogLikeSpectrum;

		/***
		 * Get intensity spectrum as multigraph.
		 * If the intensity spectrum does not jet exist, generates its.
		 * @param waveNamePattern Any wave-name pattern (sam as in fitResult intensity)
		 * @param escapedWaveNamePattern wave-name pattern used when calling fit-results intensity function
		 * @return
		 */
		componentPlot* _intensitySpectrum(const std::string& waveNamePattern);



	};

	namespace plottingtools {
		/***
		 * Centers the given phase plot in the range of [rangeStart, rangeStart+360).
		 * @param plot plot of phase spectrum that will be centered
		 * @param rangeStart start of the range in degree
		 */
		void shiftPhaseSpectrumInRange(componentPlot::plotType* plot, const double rangeStart);


		/***
		 * Tries to avoid phase jumps when plotting the phases within +-180 deg, by shifting the points individually by n*360.
		 *
		 * This is done by trying different ranges for the phase spectrum from [-220,140) until [40,400).
		 * Then the discontinuity of phase spectrum of each range is measured as the sum of quadratic distances of point i and point (i+1).
		 * Finally the phase spectrum is shifted in the range with the smallest discontinuity.
		 * This works only with pahse spectra with error bars.
		 * @param plot plot of phase spectrum that should be made more continuous
		 * @param nTrails number of different ranges tried
		 * @return range start of the range [rangeStart, rangeStart+360) that is used to make the plot more continuous.
		 */
		double makePhaseContinousWithinRange(componentPlot::plotType* plot, const int nTrails = 50);
	}


}// namespace rpwa

#endif
