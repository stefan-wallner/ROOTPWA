#ifndef HLI_CALCINTEGRAL_H
#define HLI_CALCINTEGRAL_H

#include <utility>
#include <vector>
#include <string>

#include <ampIntegralMatrix.h>
#include <eventMetadata.h>
#include <hashCalculator.h>
#include <isobarAmplitude.h>

namespace rpwa {
	namespace hli {

		/**
		 * Calculates the integral matrix for the given amplitudes
		 * by calculating the amplitudes on the fly for the events.
		 * given in eventMeta.
		 * @param amplitudes Vector of amplitudes to calculate
		 * @param eventMeta Input events
		 * @param waveNames Vector of wave names for amplitudes
		 * @param minEvent First event to process
		 * @param maxEvent Last event to process (exclusive). If -1, all events from minEvent until the end are processed
		 * @param multibinBoundaries Boundaries for on-the-fly binning
		 * @param treeCacheSize Size for the input tree cache
		 * @return Integral matrix and the hasher for each amplitude
		 */
		std::pair<rpwa::ampIntegralMatrix, std::vector<rpwa::hashCalculator>>
		calcIntegralOnTheFly(const rpwa::eventMetadata& eventMeta,
		                     const std::vector<rpwa::isobarAmplitudePtr>& amplitudes,
		                     rpwa::multibinBoundariesType& multibinBoundaries,
		                     const long int minEvent = 0,
		                     long int maxEvent = -1,
		                     const long int treeCacheSize = 25000000
		                     );

	}
}

#endif // HLI_CALCINTEGRAL_H
