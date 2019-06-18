#include "calcIntegral_py.h"

#include <boost/python.hpp>
#include "stlContainers_py.h"

#include "calcIntegral.h"

namespace bp = boost::python;


namespace {

	bp::tuple
		calcIntegralOnTheFly( const rpwa::eventMetadata& eventMeta,
		                      bp::list& pyAmplitudes,
		                      bp::dict& pyMultibinBoundaries,
		                      long int minEvent,
		                      long int maxEvent,
		                      const long int treeCacheSize)
	{
		rpwa::multibinBoundariesType multibinBoundaries = rpwa::py::convertMultibinBoundariesFromPy(pyMultibinBoundaries);
		std::vector<rpwa::isobarAmplitudePtr> amplitudes;
		rpwa::py::convertBPObjectToVector(pyAmplitudes, amplitudes);
		auto result = rpwa::hli::calcIntegralOnTheFly(eventMeta, amplitudes, multibinBoundaries, minEvent, maxEvent, treeCacheSize);
		bp::list pyHashers;
		for(const auto& hasher: result.second) pyHashers.append(hasher);
		return bp::make_tuple(result.first, pyHashers);
	}

}


void rpwa::py::exportCalcIntegral()
{

	bp::def(
		"calcIntegralOnTheFly"
		, &::calcIntegralOnTheFly
		, (bp::arg("eventMeta"),
		   bp::arg("amplitudes"),
		   bp::arg("pyMultibinBoundaries"),
		   bp::arg("minEvent") = 0,
		   bp::arg("maxEvent") = -1,
		   bp::arg("treeCacheSize") = 25000000)
	);

}
