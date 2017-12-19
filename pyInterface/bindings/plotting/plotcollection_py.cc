#include "plotcollection_py.h"

#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <map>
#include <vector>
#include <stdexcept>

#include <TTree.h>
#include <TPython.h>
#include "../../../plotting/plotcollection.h"
#include "fitResult.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"
#include "reportingUtils.hpp"

namespace bp = boost::python;

namespace {
	PyObject* componentPlot_root(rpwa::componentPlot& self) {
		return TPython::ObjectProxy_FromVoidPtr(&self, self.ClassName(), false);
	}


	bp::dict
	multigraph_getComponentsHandler(rpwa::componentPlot& self, const std::string& componentName) {
		bp::dict components;
		for (const auto& comp : self.getComponents(componentName)) {
			PyObject* obj = static_cast<PyObject*>(TPython::ObjectProxy_FromVoidPtr(comp.second, comp.second->ClassName(), false));
			// use borrowed such that the obj will not be destroyed when object goes out of scope
			bp::object object(bp::handle<>(bp::borrowed(obj)));
			components[comp.first] = object;
		}
		return components;
	}


	bp::dict
	componentPlot_getMassindepComponents(rpwa::componentPlot& self) {
		return multigraph_getComponentsHandler(self, rpwa::plotComponents::names[rpwa::plotComponents::massIndependent]);
	}


	bp::dict
	componentPlot_getComponents(rpwa::componentPlot& self, bp::object& componentIdentifier) {
		bp::extract<int> componentType(componentIdentifier);
		if (componentType.check()) {
			return multigraph_getComponentsHandler(self, rpwa::plotComponents::names.at(static_cast<int>(componentType)));
		} else {
			bp::extract<std::string> componentName(componentIdentifier);
			if (componentName.check()) {
				return multigraph_getComponentsHandler(self, componentName);
			} else {
				printErr<< "multigraph::getComponents can only be called with integer or string" << std::endl;
			}
		}
		return bp::dict();
	}


	bp::object
	componentPlot_getComponent(rpwa::componentPlot& self, bp::object& componentIdentifier) {
		bp::dict components = componentPlot_getComponents(self, componentIdentifier);
		if (bp::len(components) == 1) {
			return components.values()[0];
		} else {
			printErr << "Componentplot has more than one label" << std::endl;
			return bp::object();
		}
	}


	bp::object
	componentPlot_getMassindepComponent(rpwa::componentPlot& self) {
		bp::str name(rpwa::plotComponents::names[rpwa::plotComponents::massIndependent]);
		return componentPlot_getComponent(self, name);
	}


	bp::list
	componentPlot_getComponentNames(rpwa::componentPlot& self) {
		bp::list names;
		for (const auto& name : self.getComponentNames())
			names.append(name);
		return names;
	}


	void
	componentPlot_addForComponent(rpwa::componentPlot& self, bp::object& typeIdOrName, const std::string& label, PyObject* pyPlot) {
		rpwa::componentPlot::plotType* plot = rpwa::py::convertFromPy<rpwa::componentPlot::plotType*>(pyPlot);
		bp::extract<int> extractId(typeIdOrName);
		if (extractId.check()) {
			const int id = extractId;
			self.addForComponent(id, label, plot);
		} else {
			bp::extract<std::string> extractName(typeIdOrName);
			if (extractName.check()) {
				const std::string& name = extractName;
				self.addForComponent(name, label, plot);
			} else {
				printErr<< "Cannot add the given component to the componentPlot, because the component type cannot be converted" << std::endl;
				throw std::exception();
			}
		}
	}


	boost::shared_ptr<rpwa::multibinPlots>
	multibinPlots_constructor(const bp::list& fitResults, const bp::str& label, const bp::str& desc) {

		std::vector<rpwa::fitResult> cFitResults;
		for (int i = 0; i < bp::len(fitResults); ++i) {
			cFitResults.push_back(bp::extract<rpwa::fitResult>(fitResults[i]));

		}
		return boost::shared_ptr<rpwa::multibinPlots>(new rpwa::multibinPlots(cFitResults,
		        bp::extract<std::string>(label),
		        bp::extract<std::string>(desc)));
	}


	void
	multibinPlots_buildDefaultPlots(rpwa::multibinPlots& self, bp::list waveNamePatterns) {
		std::vector<std::string> waveNamePatternsC;
		rpwa::py::convertBPObjectToVector(waveNamePatterns, waveNamePatternsC);
		self.buildDefaultPlots(waveNamePatternsC);
	}


	bool
	multibinPlots_buildMultibinSummedPlots(rpwa::multibinPlots& self, bp::list& plotsInMultibins) {
		std::vector<const rpwa::multibinPlots*> plotsInMultibinsC;
		for (int i = 0; i < bp::len(plotsInMultibins); ++i) {
			plotsInMultibinsC.push_back(bp::extract<rpwa::multibinPlots*>(plotsInMultibins[i]));
		}
		return self.buildMultibinSummedPlots(plotsInMultibinsC);
	}


	void
	multibinPlots_setIntensities(rpwa::multibinPlots& self, bp::dict& pyIntensities){
		std::map<std::string, rpwa::componentPlot*> intensities;
		bp::list pyWaveNames = pyIntensities.keys();
		for(int i = 0; i < bp::len(pyWaveNames); ++i){
			const std::string waveName = bp::extract<std::string>(pyWaveNames[i]);
			rpwa::componentPlot* plot;
			bp::extract<rpwa::componentPlot*> componentsplot(pyIntensities[pyWaveNames[i]]);
			if(componentsplot.check()){
				plot = componentsplot();
			} else {
				bp::object o = pyIntensities[pyWaveNames[i]];
				plot = static_cast<rpwa::componentPlot*>(rpwa::py::convertFromPy<rpwa::componentPlot::baseType*>(o.ptr()));
				if (plot == nullptr){ // it was not a componentPlot::baseType
					printErr << "Cannot convert plot from Python to Cpp" << std::endl;
					throw std::exception();
				}
			}
			intensities[waveName] = plot;
		}
		self.setIntensitySpectra(intensities);
	}


	int
	multibinPlots_write(rpwa::multibinPlots& self, PyObject* directory) {
		TDirectory* cDirectory = rpwa::py::convertFromPy<TDirectory*>(directory);
		return self.write(cDirectory) ? 0 : 1;
	}


	bool
	multibinPlots_load(rpwa::multibinPlots& self, PyObject* directory, const int onlyBest) {
		TDirectory* cDirectory = rpwa::py::convertFromPy<TDirectory*>(directory);

		return self.load(cDirectory, onlyBest > 0);
	}


	bp::dict
	multibinPlots_fitResultsInMassbins(rpwa::multibinPlots& self) {
		std::map<double, std::vector<rpwa::fitResult> >& fitResultsInMassbinsC = self.fitResultsInMassbins();
		bp::dict fitResultsInMassbins;
		for (const auto& mass_fitResultsC : fitResultsInMassbinsC) {
			bp::list fitResults;
			for (const auto& fitResultC : mass_fitResultsC.second) {
				fitResults.append(&fitResultC);
			}
			fitResultsInMassbins[mass_fitResultsC.first] = fitResults;
		}
		return fitResultsInMassbins;
	}


	bp::list
	multibinPlots_waveNames(rpwa::multibinPlots& self) {
		bp::list waveNames;
		for (const auto& waveNameC : self.waveNames()) {
			waveNames.append(bp::str(waveNameC));
		}
		return waveNames;
	}


	bp::list
	multibinPlots_descriptions(rpwa::multibinPlots& self) {
		bp::list descriptions;
		for (const auto& desc: self.descriptions()){
			descriptions.append(bp::str(desc));
		}
		return descriptions;
	}


	bp::list
	multibinPlots_labels(rpwa::multibinPlots& self) {
		bp::list labels;
		for (const auto& label: self.labels()){
			labels.append(bp::str(label));
		}
		return labels;
	}

	bool
	multibinPlots_addAdditionalPlot(rpwa::multibinPlots& self, PyObject* plot, const bp::str& name) {
		TObject* objRoot = (TObject*) (TPython::ObjectProxy_AsVoidPtr(plot));
		TMultiGraph* mg = dynamic_cast<TMultiGraph*>(objRoot);
		TGraphErrors* g = dynamic_cast<TGraphErrors*>(objRoot);
		if (mg != nullptr) {
			return self.addAdditionalPlot(static_cast<rpwa::componentPlot*>(mg), bp::extract<std::string>(name));
		} else if (g != nullptr) {
			return self.addAdditionalPlot(g, bp::extract<std::string>(name));
		}

		printErr<< "Can not add plot of type '" << objRoot->ClassName() << "' to a multibinplot object." << std::endl;
		return false;
	}

	rpwa::componentPlot* multibinPlots_negLogLikeSpectrum(rpwa::multibinPlots& self){
		return self.negLogLikeSpectrum();
	}


	PyObject*
	multibinPlots_getAdditionalPlot(rpwa::multibinPlots& self, const bp::str& name) {
		TMultiGraph* mp = self.getAdditionalPlot(bp::extract<std::string>(name));
		return TPython::ObjectProxy_FromVoidPtr(mp, mp->ClassName(), false);
	}


	void
	plottingtools_shiftPhaseSpectrumInRange(PyObject* plot, const double rangeStart) {
		rpwa::plottingtools::shiftPhaseSpectrumInRange(rpwa::py::convertFromPy<rpwa::componentPlot::plotType*>(plot), rangeStart);
	}


	double
	plottingtools_makePhaseContinousWithinRange(PyObject* plot, const int nTrails) {
		return rpwa::plottingtools::makePhaseContinousWithinRange(rpwa::py::convertFromPy<rpwa::componentPlot::plotType*>(plot), nTrails);
	}

}


void
rpwa::py::exportPlotcollection() {
	{ // scope for componentsScope to apply to the enum
		bp::scope componentsScope = bp::class_<rpwa::plotComponents>("plotComponents")
		        .def_readonly("names", &rpwa::plotComponents::names)
		        ;
		bp::enum_<rpwa::plotComponents::IDs>("plotComponentIDs")
		        .value("massIndependent", rpwa::plotComponents::massIndependent)
		        .export_values()
		;
	}

	bp::class_<rpwa::componentPlot>("componentPlot")
	        .def("root", &componentPlot_root)
	        .def("getMassindepComponents", &componentPlot_getMassindepComponents)
	        .def("getComponents", &componentPlot_getComponents)
	        .def("getMassindepComponent", &componentPlot_getMassindepComponent)
	        .def("getComponent", &componentPlot_getComponent)
	        .def("getComponentNames", &componentPlot_getComponentNames)
	        .def("addForComponent", &componentPlot_addForComponent)
	        ;
	bp::class_<rpwa::multibinPlots, boost::shared_ptr<rpwa::multibinPlots>>("multibinPlots")
	        .def(bp::init<const rpwa::multibinPlots&>())
	        .def("__init__", bp::make_constructor(&multibinPlots_constructor))
	        .def("buildDefaultPlots", &multibinPlots_buildDefaultPlots, bp::arg("waveNamePatterns") = bp::list())
	        .def("buildMultibinSummedPlots", &multibinPlots_buildMultibinSummedPlots)
	        .def("setIntensitySpectra", &multibinPlots_setIntensities)
	        .def("setLabel", &rpwa::multibinPlots::setLabel)
	        .def("setDescription", &rpwa::multibinPlots::setDescription)
	        .def("mergePlotsInto", &rpwa::multibinPlots::mergePlotsInto)
	        .def("intensitySpectrum", &rpwa::multibinPlots::intensitySpectrum, bp::return_internal_reference<>())
	        .def("intensitySpectrumRegEx", &rpwa::multibinPlots::intensitySpectrumRegEx, bp::return_internal_reference<>())
	        .def("phaseSpectrum", &rpwa::multibinPlots::phaseSpectrum, bp::return_internal_reference<>())
	        .def("write", &multibinPlots_write)
	        .def("load", &multibinPlots_load, (bp::args("directory"), bp::args("onlyBest") = false))
	        .def("fitResultsInMassbins", &multibinPlots_fitResultsInMassbins)
	        .def("waveNames", &multibinPlots_waveNames)
	        .def("getAdditionalPlot", &multibinPlots_getAdditionalPlot)
	        .def("addAdditionalPlot", &multibinPlots_addAdditionalPlot)
	        .def("descriptions", &multibinPlots_descriptions)
	        .def("labels", &multibinPlots_labels)
	        .def("calcIntensityIntegral", &rpwa::multibinPlots::calcIntensityIntegral,
	             (bp::args("waveName"), bp::args("xmin") = nan(""), bp::args("xmax") = nan("")))
	        .def("calcIntensityIntegralRegEx", &rpwa::multibinPlots::calcIntensityIntegralRegEx,
	             (bp::args("waveNamePattern"), bp::args("xmin") = nan(""), bp::args("xmax") = nan("")))
	        .def("negLogLikeSpectrum", &multibinPlots_negLogLikeSpectrum, bp::return_internal_reference<>())
	        ;
	bp::def("shiftPhaseSpectrumInRange", &plottingtools_shiftPhaseSpectrumInRange);
	bp::def("makePhaseContinousWithinRange", &plottingtools_makePhaseContinousWithinRange, (bp::args("plot"), bp::arg("nTrails") = 50));

}
