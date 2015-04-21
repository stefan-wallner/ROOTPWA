///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert, Boris Grube
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
//
// Description:
//      Implementation of class pwaLikelihood
//      see pwaLikelihood.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//      Boris Grube          Universe Cluster Munich
//
//
//-----------------------------------------------------------


#include <iomanip>
#include <fstream>
#include <complex>
#include <cassert>
#include <limits>

#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "complexMatrix.h"
#include "reportingUtils.hpp"
#include "conversionUtils.hpp"
#include "fileUtils.hpp"
#ifdef USE_CUDA
#include "complex.cuh"
#include "likelihoodInterface.cuh"
#endif
#include "amplitudeTreeLeaf.h"
#include "pwaLikelihood.h"


// #define USE_FDF


using namespace std;
using namespace rpwa;
using namespace boost;
using namespace boost::accumulators;
namespace bt = boost::tuples;


template<typename complexT> const double pwaLikelihood<complexT>::_cauchyWidth = 0.5;
template<typename complexT> bool pwaLikelihood<complexT>::_debug = true;


namespace {

	double cauchyFunction(const double x, const double gamma)
	{
		if(x < 0.) {
			printWarn << "got negative argument." << endl;
		}
		return 1. / (1. + ((x*x) / (gamma*gamma)));
	}

	double cauchyFunctionDerivative(const double x, const double gamma)
	{
		if(x < 0.) {
			printWarn << "got negative argument." << endl;
		}
		return (-2. * x * gamma*gamma) / ((gamma*gamma + x*x) * (gamma*gamma + x*x));
	}

	double cauchyFunctionSecondDerivative(const double x, const double gamma)
	{
		if(x < 0.) {
			printWarn << "got negative argument." << endl;
		}
		const double x2 = x*x;
		const double gamma2 = gamma*gamma;
		return (6. * x2 * gamma2 - 2 * gamma2*gamma2) / ((gamma2 + x2) * (gamma2 + x2) + (gamma2 + x2));
	}

}


template<typename complexT>
pwaLikelihood<complexT>::pwaLikelihood()
	: _nmbEvents        (0),
	  _rank             (1),
	  _nmbWaves         (0),
	  _nmbWavesReflMax  (0),
	  _nmbPars          (0),
#ifdef USE_CUDA
	  _cudaEnabled      (false),
#endif
	  _useNormalizedAmps(true),
	  _priorType        (FLAT),
	  _numbAccEvents    (0)
{
	_nmbWavesRefl[0] = 0;
	_nmbWavesRefl[1] = 0;
	resetFuncCallInfo();
#ifdef USE_FDF
	printInfo << "using FdF() to calculate likelihood" << endl;
#else
	printInfo << "using DoEval() to calculate likelihood" << endl;
#endif
}


template<typename complexT>
pwaLikelihood<complexT>::~pwaLikelihood()
{
	clear();
}


template<typename complexT>
void
pwaLikelihood<complexT>::FdF
(const double* par,             // parameter array; reduced by rank conditions
 double&       funcVal,         // function value
 double*       gradient) const  // array of derivatives
{
	++(_funcCallInfo[FDF].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// copy arguments into parameter cache
	for (unsigned int i = 0; i < _nmbPars; ++i)
		_parCache[i] = par[i];

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// create array of likelihood derivatives w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	value_type                                     derivativeFlat = 0;
	boost::array<typename ampsArrayType::index, 3> derivShape     = {{ _rank, 2, _nmbWavesReflMax }};
	ampsArrayType                                  derivatives(derivShape);

	// loop over events and calculate real-data term of log likelihood
	// as well as derivatives with respect to parameters
	TStopwatch timer;
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > logLikelihoodAcc;
	accumulator_set<value_type, stats<tag::sum(compensated)> > derivativeFlatAcc;
	multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
		derivativesAcc(derivShape);
	for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
		accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
		ampsArrayType derivative(derivShape);  // likelihood derivatives for this event
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
				}
				const complexT ampProdSum = sum(ampProdAcc);
				likelihoodAcc(norm(ampProdSum));
				// set derivative term that is independent on derivative wave index
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][iWave] = ampProdSum;
			}
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iEvt][iRefl][iWave]);
		}  // end loop over rank
		likelihoodAcc   (prodAmpFlat2            );
		logLikelihoodAcc(-log(sum(likelihoodAcc)));
		// incorporate factor 2 / sigma
		const value_type factor = 2. / sum(likelihoodAcc);
		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivativesAcc[iRank][iRefl][iWave](-factor * derivative[iRank][iRefl][iWave]);
		derivativeFlatAcc(-factor * prodAmpFlat);
	}  // end loop over events
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				derivatives[iRank][iRefl][iWave] = sum(derivativesAcc[iRank][iRefl][iWave]);
	derivativeFlat = sum(derivativeFlatAcc);
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[FDF].funcTime(timer.RealTime());

	// compute normalization term of log likelihood and normalize derivatives w.r.t. parameters
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > normFactorAcc;
	const value_type nmbEvt      = (_useNormalizedAmps) ? 1 : _nmbEvents;
	const value_type twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				accumulator_set<complexT, stats<tag::sum(compensated)> > normFactorDerivAcc;
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					normFactorAcc(real((prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave]))
					                   * I));
					normFactorDerivAcc(prodAmps[iRank][iRefl][jWave] * conj(I));
				}
				derivatives[iRank][iRefl][iWave] += sum(normFactorDerivAcc) * twiceNmbEvt;  // account for 2 * nmbEvents
			}
	// take care of flat wave
	normFactorAcc(prodAmpFlat2 * _totAcc);
	derivativeFlat += prodAmpFlat * twiceNmbEvt * _totAcc;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[FDF].normTime(timer.RealTime());

	double priorValue = 0.;
	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);
						const double cauchyFunctionValue = cauchyFunction(r, _cauchyWidth);
						const double factor = (1./(r*cauchyFunctionValue)) * cauchyFunctionDerivative(r, _cauchyWidth);
						complexT derivative(factor*prodAmps[iRank][iRefl][iWave].real(), factor*prodAmps[iRank][iRefl][iWave].imag());
						priorValue -= log(cauchyFunctionValue);
						derivatives[iRank][iRefl][iWave] -= derivative;
					}
				}
			}
			break;
	}

	// sort derivative results into output array and cache
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, _derivCache.data());

	// calculate log likelihood value
	funcVal = sum(logLikelihoodAcc) + nmbEvt * sum(normFactorAcc) + priorValue;

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[FDF].totalTime(timerTot.RealTime());

	if (_debug)
		printDebug << "raw log likelihood = "        << maxPrecisionAlign(sum(logLikelihoodAcc)) << ", "
		           << "normalization = "             << maxPrecisionAlign(sum(normFactorAcc)   ) << ", "
		           << "prior = "                     << maxPrecisionAlign(priorValue           ) << ", "
		           << "normalized log likelihood = " << maxPrecisionAlign(funcVal              ) << endl;
}


template<typename complexT>
double
pwaLikelihood<complexT>::DoEval(const double* par) const
{
	++(_funcCallInfo[DOEVAL].nmbCalls);

#ifdef USE_FDF

	// call FdF
	double logLikelihood;
	double gradient[_nmbPars];
	FdF(par, logLikelihood, gradient);
	return logLikelihood;

#else  // USE_FDF

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// loop over events and calculate real-data term of log likelihood
	TStopwatch timer;
	timer.Start();
	value_type logLikelihood = 0;
#ifdef USE_CUDA
	if (_cudaEnabled) {
		logLikelihood = cuda::likelihoodInterface<cuda::complex<value_type> >::logLikelihood
			(reinterpret_cast<cuda::complex<value_type>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, _rank);
	} else
#endif
	{
		accumulator_set<value_type, stats<tag::sum(compensated)> > logLikelihoodAcc;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
					}
					const complexT ampProdSum = sum(ampProdAcc);
					likelihoodAcc(norm(ampProdSum));
				}
			}  // end loop over rank
			likelihoodAcc   (prodAmpFlat2            );
			logLikelihoodAcc(-log(sum(likelihoodAcc)));
		}  // end loop over events
		logLikelihood = sum(logLikelihoodAcc);
	}
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[DOEVAL].funcTime(timer.RealTime());

	// compute normalization term of log likelihood
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > normFactorAcc;
	const value_type nmbEvt = (_useNormalizedAmps) ? 1 : _nmbEvents;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					normFactorAcc(real((prodAmps[iRank][iRefl][iWave] * conj(prodAmps[iRank][iRefl][jWave]))
					                   * I));
				}
			}
	// take care of flat wave
	normFactorAcc(prodAmpFlat2 * _totAcc);
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[DOEVAL].normTime(timer.RealTime());

	double priorValue = 0.;
	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);
						const double cauchyFunctionValue = cauchyFunction(r, _cauchyWidth);
						priorValue -= log(cauchyFunctionValue);
					}
				}
			}
			break;
	}

	// calculate log likelihood value
	const double funcVal = logLikelihood + nmbEvt * sum(normFactorAcc) + priorValue;

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[DOEVAL].totalTime(timerTot.RealTime());

	if (_debug)
		printDebug << "raw log likelihood = "        << maxPrecisionAlign(logLikelihood     ) << ", "
		           << "normalization = "             << maxPrecisionAlign(sum(normFactorAcc)) << ", "
		           << "prior = "                     << maxPrecisionAlign(priorValue        ) << ", "
		           << "normalized log likelihood = " << maxPrecisionAlign(funcVal           ) << endl;

	return funcVal;

#endif  // USE_FDF
}


template<typename complexT>
double
pwaLikelihood<complexT>::DoDerivative(const double* par,
                                      unsigned int  derivativeIndex) const
{
	++(_funcCallInfo[DODERIVATIVE].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// check whether parameter is in cache
	bool samePar = true;
	for (unsigned int i = 0; i < _nmbPars; ++i)
		if (_parCache[i] != par[i]) {
			samePar = false;
			break;
		}
	timerTot.Stop();
	_funcCallInfo[DODERIVATIVE].totalTime(timerTot.RealTime());
	if (samePar) {
		//cout << "using cached derivative! " << endl;
		return _derivCache[derivativeIndex];
	}
	// call FdF
	double logLikelihood;
	double gradient[_nmbPars];
	FdF(par, logLikelihood, gradient);
	return gradient[derivativeIndex];
}


// calculate derivatives with respect to parameters
template<typename complexT>
void
pwaLikelihood<complexT>::Gradient
(const double* par,             // parameter array; reduced by rank conditions
 double*       gradient) const  // array of derivatives
{
	++(_funcCallInfo[GRADIENT].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

#ifdef USE_FDF

	// check whether parameter is in cache
	bool samePar = true;
	for (unsigned int i = 0; i < _nmbPars; ++i)
		if (_parCache[i] != par[i]) {
			samePar = false;
			break;
		}
	timerTot.Stop();
	_funcCallInfo[GRADIENT].totalTime(timerTot.RealTime());
	if (samePar) {
		for (unsigned int i = 0; i < _nmbPars ; ++i)
			gradient[i] = _derivCache[i];
		return;
	}
	// call FdF
	double logLikelihood;
	FdF(par, logLikelihood, gradient);

#else  // USE_FDF

	// copy arguments into parameter cache
	for (unsigned int i = 0; i < _nmbPars; ++i)
		_parCache[i] = par[i];

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);

	// create array of likelihood derivatives w.r.t. real and imaginary
	// parts of the production amplitudes
	// !NOTE! although stored as and constructed from complex values,
	// the dL themselves are _not_ well defined complex numbers!
	value_type                                     derivativeFlat = 0;
	boost::array<typename ampsArrayType::index, 3> derivShape     = {{ _rank, 2, _nmbWavesReflMax }};
	ampsArrayType                                  derivatives(derivShape);

	// loop over events and calculate derivatives with respect to parameters
	TStopwatch timer;
	timer.Start();
#ifdef USE_CUDA
	if (_cudaEnabled) {
		cuda::likelihoodInterface<cuda::complex<value_type> >::logLikelihoodDeriv
			(reinterpret_cast<cuda::complex<value_type>*>(prodAmps.data()),
			 prodAmps.num_elements(), prodAmpFlat, _rank,
			 reinterpret_cast<cuda::complex<value_type>*>(derivatives.data()),
			 derivativeFlat);
	} else
#endif
	{
		accumulator_set<value_type, stats<tag::sum(compensated)> > derivativeFlatAcc;
		multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
			derivativesAcc(derivShape);
		const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;
		for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
			accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
			ampsArrayType derivative(derivShape);  // likelihood derivatives for this event
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
					}
					const complexT ampProdSum = sum(ampProdAcc);
					likelihoodAcc(norm(ampProdSum));
					// set derivative term that is independent on derivative wave index
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
						// amplitude sums for current rank and for waves with same reflectivity
						derivative[iRank][iRefl][iWave] = ampProdSum;
				}
				// loop again over waves for current rank and multiply with complex conjugate
				// of decay amplitude of the wave with the derivative wave index
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
						derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iEvt][iRefl][iWave]);
			}  // end loop over rank
			likelihoodAcc(prodAmpFlat2);
			// incorporate factor 2 / sigma
			const value_type factor = 2. / sum(likelihoodAcc);
			for (unsigned int iRank = 0; iRank < _rank; ++iRank)
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
						derivativesAcc[iRank][iRefl][iWave](-factor * derivative[iRank][iRefl][iWave]);
			derivativeFlatAcc(-factor * prodAmpFlat);
		}  // end loop over events
		for (unsigned int iRank = 0; iRank < _rank; ++iRank)
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivatives[iRank][iRefl][iWave] = sum(derivativesAcc[iRank][iRefl][iWave]);
		derivativeFlat = sum(derivativeFlatAcc);
	}
	// log time needed for likelihood calculation
	timer.Stop();
	_funcCallInfo[GRADIENT].funcTime(timer.RealTime());

	// normalize derivatives w.r.t. parameters
	timer.Start();
	const value_type nmbEvt      = (_useNormalizedAmps) ? 1 : _nmbEvents;
	const value_type twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				accumulator_set<complexT, stats<tag::sum(compensated)> > normFactorDerivAcc;
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {  // inner loop over waves with same reflectivity
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					normFactorDerivAcc(prodAmps[iRank][iRefl][jWave] * conj(I));
				}
				derivatives[iRank][iRefl][iWave] += sum(normFactorDerivAcc) * twiceNmbEvt;  // account for 2 * nmbEvents
			}
	// take care of flat wave
	derivativeFlat += prodAmpFlat * twiceNmbEvt * _totAcc;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[GRADIENT].normTime(timer.RealTime());

	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);
						const double cauchyFunctionValue = cauchyFunction(r, _cauchyWidth);
						const double factor = (1./(r*cauchyFunctionValue)) * cauchyFunctionDerivative(r, _cauchyWidth);
						complexT derivative(factor*prodAmps[iRank][iRefl][iWave].real(), factor*prodAmps[iRank][iRefl][iWave].imag());
						derivatives[iRank][iRefl][iWave] -= derivative;
					}
				}
			}
			break;
	}

	// sort derivative results into output array and cache
	copyToParArray(derivatives, derivativeFlat, gradient);
	copyToParArray(derivatives, derivativeFlat, _derivCache.data());

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[GRADIENT].totalTime(timerTot.RealTime());

#endif  // USE_FDF
}


// calculate Hessian with respect to parameters
template<typename complexT>
TMatrixT<double>
pwaLikelihood<complexT>::Hessian
(const double* par) const  // parameter array; reduced by rank conditions
{
	++(_funcCallInfo[HESSIAN].nmbCalls);

	// timer for total time
	TStopwatch timerTot;
	timerTot.Start();

	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);
	const value_type prodAmpFlat2 = prodAmpFlat * prodAmpFlat;

	// create array to store likelihood hessian w.r.t. real and imaginary
	// parts of the production amplitudes
	value_type                                     hessianFlat  = 0;       // term in Hessian matrix where we derive twice w.r.t. the flat wave
	boost::array<typename ampsArrayType::index, 3> derivShape   = {{ _rank, 2, _nmbWavesReflMax }};
	ampsArrayType                                  flatTerms(derivShape);  // array for terms where we first derive w.r.t to
	                                                                       // a non-flat term and then w.r.t. the flat wave
	boost::array<typename ampsArrayType::index, 7> hessianShape = {{ _rank, 2, _nmbWavesReflMax, _rank, 2, _nmbWavesReflMax, 3 }};
	boost::multi_array<value_type, 7>              hessian(hessianShape);  // array to store components for Hessian matrix

	// loop over events and calculate second derivatives with respect to
	// parameters for the raw likelihood part
	TStopwatch timer;
	timer.Start();
	accumulator_set<value_type, stats<tag::sum(compensated)> > hessianFlatAcc;
	multi_array<accumulator_set<complexT, stats<tag::sum(compensated)> >, 3>
		flatTermsAcc(derivShape);
	multi_array<accumulator_set<value_type, stats<tag::sum(compensated)> >, 7>
		hessianAcc(hessianShape);
	for (unsigned int iEvt = 0; iEvt < _nmbEvents; ++iEvt) {
		accumulator_set<value_type, stats<tag::sum(compensated)> > likelihoodAcc;
		ampsArrayType derivative(derivShape);  // likelihood derivatives for this event
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
				accumulator_set<complexT, stats<tag::sum(compensated)> > ampProdAcc;
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					ampProdAcc(prodAmps[iRank][iRefl][iWave] * _decayAmps[iEvt][iRefl][iWave]);
				}
				const complexT ampProdSum = sum(ampProdAcc);
				likelihoodAcc(norm(ampProdSum));
				// set derivative term that is independent on derivative wave index
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					// amplitude sums for current rank and for waves with same reflectivity
					derivative[iRank][iRefl][iWave] = ampProdSum;
			}
			// loop again over waves for current rank and multiply with complex conjugate
			// of decay amplitude of the wave with the derivative wave index
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					derivative[iRank][iRefl][iWave] *= conj(_decayAmps[iEvt][iRefl][iWave]);
		}  // end loop over rank
		likelihoodAcc(prodAmpFlat2);
		// incorporate factor 2 / sigma
		const value_type factor  = 2. / sum(likelihoodAcc);
		const value_type factor2 = factor*factor;
		for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
					for (unsigned int jRank = 0; jRank < _rank; ++jRank) {
						for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
							for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
								// last array index 0 indicates derivative w.r.t. real part of first prodAmp and real part of the second prodAmp
								hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][0](factor2 * derivative[jRank][jRefl][jWave].real() * derivative[iRank][iRefl][iWave].real());
								// last array index 1 indicates derivative w.r.t. real part of first prodAmp and imaginary part of the second prodAmp
								hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][1](factor2 * derivative[jRank][jRefl][jWave].imag() * derivative[iRank][iRefl][iWave].real());
								// last array index 2 indicates derivative w.r.t. imaginary part of first prodAmp and imaginary part of the second prodAmp
								hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][2](factor2 * derivative[jRank][jRefl][jWave].imag() * derivative[iRank][iRefl][iWave].imag());
								if(iRank == jRank and iRefl == jRefl) {
									const complexT uPrime = conj(_decayAmps[iEvt][jRefl][jWave]) * _decayAmps[iEvt][iRefl][iWave];
									hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][0](-factor * uPrime.real());
									hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][1](-factor * uPrime.imag());
									hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][2](-factor * uPrime.real());
								}
							}
						}
					}
					flatTermsAcc[iRank][iRefl][iWave](factor2 * prodAmpFlat * derivative[iRank][iRefl][iWave]);  // calculate terms where we first derive w.r.t. real/imag part
					                                                                                             // of a prodAmp and then w.r.t. the flat wave
				}
			}
		}
		hessianFlatAcc(factor2 * prodAmpFlat2 - factor);
	}  // end loop over events
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				for (unsigned int jRank = 0; jRank < _rank; ++jRank) {
					for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
						for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave){
							hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][0] = sum(hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][0]);
							hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][1] = sum(hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][1]);
							hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][2] = sum(hessianAcc[iRank][iRefl][iWave][jRank][jRefl][jWave][2]);
						}
					}
				}
				flatTerms[iRank][iRefl][iWave] = sum(flatTermsAcc[iRank][iRefl][iWave]);
			}
		}
	}
	hessianFlat = sum(hessianFlatAcc);
	// log time needed for calculation of second derivatives of raw likelhood part
	timer.Stop();
	_funcCallInfo[HESSIAN].funcTime(timer.RealTime());

	// normalize second derivatives w.r.t. parameters
	timer.Start();
	const value_type nmbEvt      = (_useNormalizedAmps) ? 1 : _nmbEvents;
	const value_type twiceNmbEvt = 2 * nmbEvt;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[iRefl]; ++jWave) {
					const complexT I = _accMatrix[iRefl][iWave][iRefl][jWave];
					hessian[iRank][iRefl][iWave][iRank][iRefl][jWave][0] += I.real() * twiceNmbEvt;
					hessian[iRank][iRefl][iWave][iRank][iRefl][jWave][1] += I.imag() * twiceNmbEvt;
					hessian[iRank][iRefl][iWave][iRank][iRefl][jWave][2] += I.real() * twiceNmbEvt;
				}
	hessianFlat += twiceNmbEvt * _totAcc;
	// log time needed for normalization
	timer.Stop();
	_funcCallInfo[HESSIAN].normTime(timer.RealTime());

	switch(_priorType)
	{
		case FLAT:
			break;
		case HALF_CAUCHY:
			for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
				for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
					for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
						const double r = abs(prodAmps[iRank][iRefl][iWave]);

						const double cauchyFunctionValue                 = cauchyFunction                (r, _cauchyWidth);
						const double cauchyFunctionDerivativeValue       = cauchyFunctionDerivative      (r, _cauchyWidth);
						const double cauchyFunctionSecondDerivativeValue = cauchyFunctionSecondDerivative(r, _cauchyWidth);

						const double a1 = cauchyFunctionSecondDerivativeValue / (cauchyFunctionValue * r*r);
						const double a2 = cauchyFunctionDerivativeValue*cauchyFunctionDerivativeValue / (cauchyFunctionValue*cauchyFunctionValue * r*r);
						const double a3 = cauchyFunctionDerivativeValue / (cauchyFunctionValue * r);
						const double a4 = cauchyFunctionDerivativeValue / (cauchyFunctionValue * r*r*r);
						const double a124 = a1 - a2 - a4;

						const double re = prodAmps[iRank][iRefl][iWave].real();
						const double im = prodAmps[iRank][iRefl][iWave].imag();

						hessian[iRank][iRefl][iWave][iRank][iRefl][iWave][0] -= (a124) * re*re + a3;
						hessian[iRank][iRefl][iWave][iRank][iRefl][iWave][1] -= (a124) * re*im;
						hessian[iRank][iRefl][iWave][iRank][iRefl][iWave][2] -= (a124) * im*im + a3;
					}
				}
			}
			break;
	}

	TMatrixT<double> hessianMatrix(_nmbPars, _nmbPars);
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				bt::tuple<int, int> parIndices1 = _prodAmpToFuncParMap[iRank][iRefl][iWave];
				const int r1 = get<0>(parIndices1);
				const int i1 = get<1>(parIndices1);

				if (r1 < 0)
					continue;
				if (_parFixed[r1])
					continue;
				assert(i1 < 0 || not _parFixed[i1]);

				for (unsigned int jRank = 0; jRank < _rank; ++jRank) {
					for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
						for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
							bt::tuple<int, int> parIndices2 = _prodAmpToFuncParMap[jRank][jRefl][jWave];
							const int r2 = get<0>(parIndices2);
							const int i2 = get<1>(parIndices2);

							if (r2 < 0)
								continue;
							if (_parFixed[r2])
								continue;
							assert(i2 < 0 || not _parFixed[i2]);

							hessianMatrix[r1][r2] = hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][0];
							if (i2 >= 0) // real/imaginary derivative
								hessianMatrix[r1][i2] = hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][1];
							if (i1 >= 0) // imaginary/real derivative
								hessianMatrix[i1][r2] = hessian[jRank][jRefl][jWave][iRank][iRefl][iWave][1];
							if (i1 >= 0 && i2 >= 0) // imaginary/imaginary derivative
								hessianMatrix[i1][i2] = hessian[iRank][iRefl][iWave][jRank][jRefl][jWave][2];
						}
					}
				}

				// second derivatives w.r.t. flat wave
				hessianMatrix[r1][_nmbPars - 1] = flatTerms[iRank][iRefl][iWave].real();
				hessianMatrix[_nmbPars - 1][r1] = flatTerms[iRank][iRefl][iWave].real();
				if (i1 >= 0) {
					hessianMatrix[i1][_nmbPars - 1] = flatTerms[iRank][iRefl][iWave].imag();
					hessianMatrix[_nmbPars - 1][i1] = flatTerms[iRank][iRefl][iWave].imag();
				}
			}
		}
	}
	hessianMatrix[_nmbPars - 1][_nmbPars - 1] = hessianFlat;  // enter flat/flat term

	// log total consumed time
	timerTot.Stop();
	_funcCallInfo[HESSIAN].totalTime(timerTot.RealTime());

	return hessianMatrix;
}


/// calculates covariance matrix of function at point defined by par
template<typename complexT>
TMatrixT<double>
pwaLikelihood<complexT>::CovarianceMatrix(const double* par) const
{
	const TMatrixT<double> hessian = Hessian(par);
	return CovarianceMatrix(hessian);
}


/// turns hessian into covariance matrix
template<typename complexT>
TMatrixT<double>
pwaLikelihood<complexT>::CovarianceMatrix(const TMatrixT<double>& hessian) const
{
	// reduce Hessian matrix by removing the rows and columns corresponding
	// to fixed parameters
	TMatrixT<double> hessianRed(_nmbPars - _nmbParsFixed, _nmbPars - _nmbParsFixed);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < _nmbPars; ++i) {
			if (_parFixed[i]) {
				iSkip++;
				continue;
			}

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < _nmbPars; ++j) {
				if (_parFixed[j]) {
					jSkip++;
					continue;
				}

				hessianRed[i - iSkip][j - jSkip] = hessian[i][j];
			}
		}
	}

	// invert Hessian to get covariance matrix
	TMatrixT<double> covarianceRed(TMatrixT<double>::kInverted, hessianRed);

	// blow up covariance matrix by adding the rows and columns corresponding
	// to fixed parameters
	TMatrixT<double> covariance(_nmbPars, _nmbPars);
	{
		unsigned int iSkip = 0;
		for (unsigned int i = 0; i < _nmbPars; ++i) {
			if (_parFixed[i]) {
				iSkip++;
				continue;
			}

			unsigned int jSkip = 0;
			for (unsigned int j = 0; j < _nmbPars; ++j) {
				if (_parFixed[j]) {
					jSkip++;
					continue;
				}

				covariance[i][j] = covarianceRed[i - iSkip][j - jSkip];
			}
		}
	}

	return covariance;
}


template<typename complexT>
std::vector<double>
pwaLikelihood<complexT>::CorrectParamSigns(const double* par) const
{
	// build complex production amplitudes from function parameters taking into account rank restrictions
	value_type    prodAmpFlat;
	ampsArrayType prodAmps;
	copyFromParArray(par, prodAmps, prodAmpFlat);

	if (prodAmpFlat < 0) {
		prodAmpFlat *= -1;  // flip sign of flat wave if it is negative
	}

	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {  // incoherent sum over ranks
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {  // incoherent sum over reflectivities
			if (prodAmps[iRank][iRefl][iRank].real() < 0 and prodAmps[iRank][iRefl][iRank].imag() == 0) {
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {  // coherent sum over waves
					prodAmps[iRank][iRefl][iWave] *= -1;  // flip sign of coherent waves if anchorwave is negative
				}
			}
		}
	}

	std::vector<double> returnParams(_nmbPars);
	copyToParArray(prodAmps, prodAmpFlat, returnParams.data());
	return returnParams;
}


template<typename complexT>
unsigned int
pwaLikelihood<complexT>::nmbWaves(const int reflectivity) const
{
	if (reflectivity == 0)
		return _nmbWaves;
	else if (reflectivity > 0)
		return _nmbWavesRefl[1];  // positive reflectivity
	else
		return _nmbWavesRefl[0];  // negative reflectivity
}


template<typename complexT>
void
#ifdef USE_CUDA
pwaLikelihood<complexT>::enableCuda(const bool enableCuda)
{
	_cudaEnabled = enableCuda;
}
#else
pwaLikelihood<complexT>::enableCuda(const bool) { }
#endif


template<typename complexT>
bool
pwaLikelihood<complexT>::cudaEnabled() const
{
#ifdef USE_CUDA
	return _cudaEnabled;
#else
	return false;
#endif
}


template<typename complexT>
void
pwaLikelihood<complexT>::init(const unsigned int rank,
                              const double       massBinCenter,
                              const std::string& waveListFileName,
                              const std::string& normIntFileName,
                              const std::string& accIntFileName,
                              const std::string& ampDirName,
                              const unsigned int numbAccEvents,
                              const bool         useRootAmps)
{
	_numbAccEvents = numbAccEvents;
	readWaveList(waveListFileName);
	buildParDataStruct(rank, massBinCenter);
	readIntegrals(normIntFileName, accIntFileName);
	readDecayAmplitudes(ampDirName, useRootAmps);
#ifdef USE_CUDA
	if (_cudaEnabled)
		cuda::likelihoodInterface<cuda::complex<value_type> >::init
			(reinterpret_cast<cuda::complex<value_type>*>(_decayAmps.data()),
			 _decayAmps.num_elements(), _nmbEvents, _nmbWavesRefl, true);
#endif
}


template<typename complexT>
void
pwaLikelihood<complexT>::readWaveList(const string& waveListFileName)
{
	printInfo << "reading amplitude names and thresholds from wave list file "
	          << "'" << waveListFileName << "'." << endl;
	ifstream waveListFile(waveListFileName.c_str());
	if (not waveListFile) {
		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
		throw;
	}
	vector<string>       waveNames     [2];
	vector<double>       waveThresholds[2];
	unsigned int         countWave = 0;
	unsigned int         lineNmb   = 0;
	string               line;
	while (getline(waveListFile, line)) {
		if (line[0] == '#')  // comments start with #
			continue;
		stringstream lineStream;
		lineStream.str(line);
		string waveName;
		if (lineStream >> waveName) {
			double threshold;
			// !!! it would be safer to make the threshold value in the wave list file mandatory
			if (not (lineStream >> threshold))
				threshold = 0;
			if (_debug)
				printDebug << "reading line " << setw(3) << lineNmb + 1 << ": " << waveName<< ", "
				           << "threshold = " << setw(4) << threshold << " MeV/c^2" << endl;
			if (getReflectivity(waveName) > 0) {
				++_nmbWavesRefl[1];  // positive reflectivity
				waveNames     [1].push_back(waveName);
				waveThresholds[1].push_back(threshold);
			} else {
				++_nmbWavesRefl[0];  // negative reflectivity
				waveNames     [0].push_back(waveName);
				waveThresholds[0].push_back(threshold);
			}
			++countWave;
		} else
			printWarn << "cannot parse line '" << line << "' in wave list file "
			          << "'" << waveListFileName << "'" << endl;
		++lineNmb;
	}
	waveListFile.close();
	printInfo << "read " << lineNmb << " lines from wave list file "
	          << "'" << waveListFileName << "'" << endl;
	_nmbWaves        = _nmbWavesRefl[0] + _nmbWavesRefl[1];
	_nmbWavesReflMax = max(_nmbWavesRefl[0], _nmbWavesRefl[1]);
	_waveNames.resize      (extents[2][_nmbWavesReflMax]);
	_waveThresholds.resize (extents[2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			_waveNames      [iRefl][iWave] = waveNames     [iRefl][iWave];
			_waveThresholds [iRefl][iWave] = waveThresholds[iRefl][iWave];
		}
}


template<typename complexT>
void
pwaLikelihood<complexT>::buildParDataStruct(const unsigned int rank,
                                            const double       massBinCenter)
{
	if ((_nmbWavesRefl[0] + _nmbWavesRefl[1] == 0) or (_waveThresholds.size() == 0)) {
		printErr << "no wave info. was readWaveList() executed successfully? Aborting...";
		throw;
	}
	_rank = rank;
	// calculate dimension of function taking into account rank restrictions and flat wave
	_nmbPars = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		const int nmbProdAmpsPos  = _nmbWavesRefl[1] - iRank;  // number non-zero production amplitudes in this rank with positive reflectivity
		int       nmbProdAmpsPosC = nmbProdAmpsPos - 1;        // number of complex-valued production amplitudes in this rank with positive reflectivity
		if (nmbProdAmpsPosC < 0)
			nmbProdAmpsPosC = 0;
		const int nmbProdAmpsNeg  = _nmbWavesRefl[0] - iRank;  // number non-zero production amplitudes in this rank with negative reflectivity
		int       nmbProdAmpsNegC = nmbProdAmpsNeg - 1;        // number of complex-valued production amplitudes in this rank with negative reflectivity
		if (nmbProdAmpsNegC < 0)
			nmbProdAmpsNegC = 0;
		_nmbPars += 2 * (nmbProdAmpsPosC + nmbProdAmpsNegC);  // 2 parameters for each complex production amplitude
		// 1 real production amplitude per rank and reflectivity
		if (nmbProdAmpsPos > 0)
			++_nmbPars;
		if (nmbProdAmpsNeg > 0)
			++_nmbPars;
	}
	_nmbPars += 1;  // additonal flat wave
	printInfo << "dimension of likelihood function is " << _nmbPars << "." << endl;
	_nmbParsFixed = 0;
	_parNames.resize     (_nmbPars, "");
	_parThresholds.resize(_nmbPars, 0);
	_parFixed.resize     (_nmbPars, false);
	_parCache.resize     (_nmbPars, 0);
	_derivCache.resize   (_nmbPars, 0);
	_prodAmpToFuncParMap.resize(extents[_rank][2][_nmbWavesReflMax]);
	// build parameter names
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				ostringstream parName;
				if (iWave < iRank)  // production amplitude is zero
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = bt::make_tuple(-1, -1);
				else if (iWave == iRank) {  // production amplitude is real
					parName << "V" << iRank << "_" << _waveNames[iRefl][iWave] << "_RE";
					_parNames     [parIndex]                  = parName.str();
					_parThresholds[parIndex]                  = _waveThresholds[iRefl][iWave];
					if (_parThresholds[parIndex] != 0 && _parThresholds[parIndex] >= massBinCenter) {
						_parFixed[parIndex] = true;
						++_nmbParsFixed;
					}
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = bt::make_tuple(parIndex, -1);
					++parIndex;
				} else {  // production amplitude is complex
					parName << "V" << iRank << "_" << _waveNames[iRefl][iWave];
					_parNames     [parIndex]                  = parName.str() + "_RE";
					_parThresholds[parIndex]                  = _waveThresholds[iRefl][iWave];
					if (_parThresholds[parIndex] != 0 && _parThresholds[parIndex] >= massBinCenter) {
						_parFixed[parIndex] = true;
						++_nmbParsFixed;
					}
					_parNames     [parIndex + 1]              = parName.str() + "_IM";
					_parThresholds[parIndex + 1]              = _waveThresholds[iRefl][iWave];
					if (_parThresholds[parIndex + 1] != 0 && _parThresholds[parIndex + 1] >= massBinCenter) {
						_parFixed[parIndex + 1] = true;
						++_nmbParsFixed;
					}
					_prodAmpToFuncParMap[iRank][iRefl][iWave] = bt::make_tuple(parIndex, parIndex + 1);
					parIndex += 2;
				}
			}
	// flat wave
	_parNames     [parIndex] = "V_flat";
	_parThresholds[parIndex] = 0;
}


// returns integral matrix reordered according to _waveNames array
template<typename complexT>
void
pwaLikelihood<complexT>::reorderIntegralMatrix(const ampIntegralMatrix& integral,
                                               normMatrixArrayType&     reorderedMatrix) const
{
	// create reordered matrix
	reorderedMatrix.resize(extents[2][_nmbWavesReflMax][2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
					const complex<double> val = integral.element(_waveNames[iRefl][iWave],
					                                             _waveNames[jRefl][jWave]);
					reorderedMatrix[iRefl][iWave][jRefl][jWave] = complexT(val.real(), val.imag());
				}
}


template<typename complexT>
void
pwaLikelihood<complexT>::readIntegrals
(const string& normIntFileName,   // name of file with normalization integrals
 const string& accIntFileName,    // name of file with acceptance integrals
 const string& integralTKeyName)  // name of TKey which stores integral in .root file
{
	printInfo << "loading normalization integral from '" << normIntFileName << "'" << endl;
	const string normIntFileExt  = extensionFromPath(normIntFileName);
	if (normIntFileExt == "root") {
		TFile* intFile  = TFile::Open(normIntFileName.c_str(), "READ");
		if (not intFile or intFile->IsZombie()) {
			printErr << "could not open normalization integral file '" << normIntFileName << "'. "
			         << "Aborting..." << endl;
			throw;
		}
		ampIntegralMatrix* integral = 0;
		intFile->GetObject(integralTKeyName.c_str(), integral);
		if (not integral) {
			printErr << "cannot find integral object in TKey '" << integralTKeyName << "' in file "
			         << "'" << normIntFileName << "'. Aborting..." << endl;
			throw;
		}
		reorderIntegralMatrix(*integral, _normMatrix);
		intFile->Close();
	} else if(normIntFileExt == "int") {
			ampIntegralMatrix integral;
			integral.readAscii(normIntFileName);
			reorderIntegralMatrix(integral, _normMatrix);
	} else {
		printErr << "unknown file type '" << normIntFileName << "'. "
		         << "only .int and .root files are supported. Aborting..." << endl;
		throw;
	}

	printInfo << "loading acceptance integral from '" << accIntFileName << "'" << endl;
	const string accIntFileExt  = extensionFromPath(accIntFileName);
	if (accIntFileExt == "root") {
		TFile* intFile  = TFile::Open(accIntFileName.c_str(), "READ");
		if (not intFile or intFile->IsZombie()) {
			printErr << "could not open normalization integral file '" << accIntFileName << "'. "
			         << "Aborting..." << endl;
			throw;
		}
		ampIntegralMatrix* integral = 0;
		intFile->GetObject(integralTKeyName.c_str(), integral);
		if (not integral) {
			printErr << "cannot find integral object in TKey '" << integralTKeyName << "' in file "
			         << "'" << accIntFileName << "'. Aborting..." << endl;
			throw;
		}
		if (_numbAccEvents != 0) {
			_totAcc = ((double)integral->nmbEvents()) / (double)_numbAccEvents;
			printInfo << "total acceptance in this bin: " << _totAcc << endl;
			integral->setNmbEvents(_numbAccEvents);
		} else
			_totAcc = 1;
		reorderIntegralMatrix(*integral, _accMatrix);
		intFile->Close();
	} else if (accIntFileExt == "int") {
			ampIntegralMatrix integral;
			integral.readAscii(accIntFileName);
			if (_numbAccEvents != 0) {
				_totAcc = ((double)integral.nmbEvents()) / (double)_numbAccEvents;
				printInfo << "total acceptance in this bin: " << _totAcc << endl;
				integral.setNmbEvents(_numbAccEvents);
			} else {
				_totAcc = 1;
			}
			reorderIntegralMatrix(integral, _accMatrix);
	} else {
		printErr << "unknown file type '" << accIntFileName << "'. "
		         << "only .int and .root files are supported. exiting." << endl;
		throw;
	}
}


template<typename complexT>
void
pwaLikelihood<complexT>::readDecayAmplitudes(const string& ampDirName,
                                             const bool    useRootAmps,
                                             const string& ampLeafName)
{
	// check that normalization integrals are loaded
	if (_normMatrix.num_elements() == 0) {
		printErr << "normalization integrals have to be loaded before loading the amplitudes. "
		         << "Aborting..." << endl;
		throw;
	}
	clear();

	printInfo << "loading amplitude data from " << ((useRootAmps) ? ".root" : ".amp")
	          << " files" << endl;
	// loop over amplitudes and read in data
	bool         firstWave = true;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			// get normalization
			const complexT   normInt = _normMatrix[iRefl][iWave][iRefl][iWave];
			vector<complexT> amps;
			if (!firstWave)  // number of events is known except for first wave that is read in
				amps.reserve(_nmbEvents);
			// read decay amplitudes
			string ampFilePath = ampDirName + "/" + _waveNames[iRefl][iWave];
			if (useRootAmps) {
				ampFilePath = changeFileExtension(ampFilePath, ".root");
				printInfo << "loading amplitude data from '" << ampFilePath << "'" << endl;
				// open amplitude file
				TFile* ampFile = TFile::Open(ampFilePath.c_str(), "READ");
				if (not ampFile or ampFile->IsZombie()) {
					printWarn << "cannot open amplitude file '" << ampFilePath << "'. skipping." << endl;
					continue;
				}
				// find amplitude tree
				TTree*       ampTree     = 0;
				const string ampTreeName = changeFileExtension(_waveNames[iRefl][iWave], ".amp");
				ampFile->GetObject(ampTreeName.c_str(), ampTree);
				if (not ampTree) {
					printWarn << "cannot find tree '" << ampTreeName << "' in file "
					          << "'" << ampFilePath << "'. skipping." << endl;
					continue;
				}
				// connect tree leaf
				amplitudeTreeLeaf* ampTreeLeaf = 0;
				ampTree->SetBranchAddress(ampLeafName.c_str(), &ampTreeLeaf);
				for (long int eventIndex = 0; eventIndex < ampTree->GetEntriesFast(); ++eventIndex) {
					ampTree->GetEntry(eventIndex);
					if (!ampTreeLeaf) {
						printWarn << "null pointer to amplitude leaf for event " << eventIndex << ". "
						          << "skipping." << endl;
						continue;
					}
					assert(ampTreeLeaf->nmbIncohSubAmps() == 1);
					complexT amp(ampTreeLeaf->incohSubAmp(0).real(), ampTreeLeaf->incohSubAmp(0).imag());
					if (_useNormalizedAmps)         // normalize data, if option is switched on
						amp /= sqrt(normInt.real());  // rescale decay amplitude
					amps.push_back(amp);
				}
			} else {
				printInfo << "loading amplitude data from '" << ampFilePath << "'" << endl;
				ifstream ampFile(ampFilePath.c_str());
				if (not ampFile) {
					printErr << "cannot open amplitude file '" << ampFilePath << "'. Aborting..." << endl;
					throw;
				}
				complexT amp;
				while (ampFile.read((char*)&amp, sizeof(complexT))) {
				  if (_useNormalizedAmps) {        // normalize data, if option is switched on
					        value_type mynorm=normInt.real();
					        if(mynorm==0)mynorm=1;
						amp /= sqrt(mynorm);  // rescale decay amplitude
				  }
				  amps.push_back(amp);
				}
			}

			unsigned int nmbEvents = amps.size();
			if (firstWave) {
				_nmbEvents = nmbEvents;
				_decayAmps.resize(extents[_nmbEvents][2][_nmbWavesReflMax]);
				firstWave = false;
			} else {
				if (nmbEvents != _nmbEvents)
					printWarn << "size mismatch in amplitude files: this file contains " << nmbEvents
					          << " events, previous file had " << _nmbEvents << " events." << endl;
				// make sure not to store more events than _decayAmps can keep
				nmbEvents = std::min(nmbEvents, _nmbEvents);
			}

			// copy decay amplitudes into array that is indexed [event index][reflectivity][wave index]
			// this index scheme ensures a more linear memory access pattern in the likelihood function
			for (unsigned int iEvt = 0; iEvt < nmbEvents; ++iEvt)
				_decayAmps[iEvt][iRefl][iWave] = amps[iEvt];
			if (_debug)
				printDebug << "read " << nmbEvents << " events from file "
				           << "'" << _waveNames[iRefl][iWave] << "'" << endl;
		}
	printInfo << "loaded decay amplitudes for " << _nmbEvents << " events into memory" << endl;

	// save phase space integrals
	_phaseSpaceIntegral.resize(extents[2][_nmbWavesReflMax]);
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			_phaseSpaceIntegral[iRefl][iWave] = sqrt(_normMatrix[iRefl][iWave][iRefl][iWave].real());

	// rescale integrals, if necessary
	if (_useNormalizedAmps) {
		// matrices _normMatrix and _accMatrix are already normalized to number of Monte Carlo events
		printInfo << "rescaling integrals" << endl;
		// rescale normalization and acceptance integrals
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				value_type norm_i = sqrt(_normMatrix[iRefl][iWave][iRefl][iWave].real());
				if(norm_i==0)norm_i=1;
				for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
					for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
						value_type norm_j = sqrt(_normMatrix[jRefl][jWave][jRefl][jWave].real());
						// protect against empty amplitudes (which will not contribute anyway)
						if(norm_j==0)norm_j=1;
						if ((iRefl != jRefl) or (iWave != jWave))
							// set diagonal terms later so that norm_i,j stay unaffected
							_normMatrix[iRefl][iWave][jRefl][jWave] /= norm_i * norm_j;
						_accMatrix[iRefl][iWave][jRefl][jWave] /= norm_i * norm_j;
					}
			}
		// set diagonal elements of normalization matrix
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
				_normMatrix[iRefl][iWave][iRefl][iWave] = 1;  // diagonal term
		if (_debug) {
			printDebug << "normalized integral matrices" << endl;
			for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
				for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
					for (unsigned int jRefl = 0; jRefl < 2; ++jRefl)
						for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
							cout << "    normalization matrix [" << sign((int)iRefl * 2 - 1) << ", "
							     << setw(3) << iWave << ", " << sign((int)jRefl * 2 - 1) << ", "
							     << setw(3) << jWave << "] = "
							     << "("  << maxPrecisionAlign(_normMatrix[iRefl][iWave][jRefl][jWave].real())
							     << ", " << maxPrecisionAlign(_normMatrix[iRefl][iWave][jRefl][jWave].imag())
							     << "), acceptance matrix [" << setw(3) << iWave << ", "
							     << setw(3) << jWave << "] = "
							     << "("  << maxPrecisionAlign(_accMatrix[iRefl][iWave][jRefl][jWave].real())
							     << ", " << maxPrecisionAlign(_accMatrix[iRefl][iWave][jRefl][jWave].imag()) << ")"
							     << endl;
						}
		}
	}  // _useNormalizedAmps
}


template<typename complexT>
void
pwaLikelihood<complexT>::getIntegralMatrices(complexMatrix&  normMatrix,
                                             complexMatrix&  accMatrix,
                                             vector<double>& phaseSpaceIntegral) const
{
	phaseSpaceIntegral.clear();
	phaseSpaceIntegral.resize(_nmbWaves + 1, 0);
	unsigned int iIndex = 0;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
			phaseSpaceIntegral[iIndex] = _phaseSpaceIntegral[iRefl][iWave];
			unsigned int jIndex = 0;
			for (unsigned int jRefl = 0; jRefl < 2; ++jRefl) {
				for (unsigned int jWave = 0; jWave < _nmbWavesRefl[jRefl]; ++jWave) {
					const complexT     normVal = _normMatrix[iRefl][iWave][jRefl][jWave];
					const complexT     accVal  = _accMatrix [iRefl][iWave][jRefl][jWave];
					normMatrix.set(iIndex, jIndex, complex<double>(normVal.real(), normVal.imag()));
					accMatrix.set (iIndex, jIndex, complex<double>(accVal.real(),  accVal.imag() ));
					++jIndex;
				}
			}
			++iIndex;
		}
	}
	// set unused entries to 0
	for (unsigned int i = 0; i < normMatrix.nCols(); ++i) {
		normMatrix.set(_nmbWaves, i, complexT(0., 0.));
		normMatrix.set(i, _nmbWaves, complexT(0., 0.));
		accMatrix.set(_nmbWaves, i, complexT(0., 0.));
		accMatrix.set(i, _nmbWaves, complexT(0., 0.));
	}
	// add flat
	normMatrix.set(_nmbWaves, _nmbWaves, complexT(1,0));
	accMatrix.set (_nmbWaves, _nmbWaves, complexT(_totAcc,0));
	phaseSpaceIntegral[_nmbWaves] = 1;
}


// builds complex numbers from parameters
// maps real and imaginary part of amplitudes to error matrix
// for both rank restrictions are taken into account
template<typename complexT>
void
pwaLikelihood<complexT>::buildProdAmpArrays(const double*             inPar,
                                            vector<complex<double> >& prodAmps,
                                            vector<pair<int, int> >&  parIndices,
                                            vector<string>&           prodAmpNames,
                                            const bool                withFlat) const
{
	prodAmps.clear();
	parIndices.clear();
	prodAmpNames.clear();
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank) {
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl) {
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				double re, im;
				if (iWave < iRank)  // zero production amplitude
					continue;
				else if (iWave == iRank) {  // real production amplitude
					parIndices.push_back(make_pair(parIndex, -1));
					re = inPar[parIndex];
					im = 0;
					++parIndex;
				} else {  // complex production amplitude
					parIndices.push_back(make_pair(parIndex, parIndex + 1));
					re = inPar[parIndex];
					im = inPar[parIndex + 1];
					parIndex += 2;
				}
				prodAmps.push_back(complex<double>(re, im));
				stringstream prodAmpName;
				prodAmpName << "V" << iRank << "_" << _waveNames[iRefl][iWave];
				prodAmpNames.push_back(prodAmpName.str());
			}
		}
	}
	if (withFlat) {
		prodAmps.push_back(complex<double>(inPar[parIndex], 0));
		parIndices.push_back(make_pair(parIndex, -1));
		prodAmpNames.push_back("V_flat");
	}
}


template<typename complexT>
void
pwaLikelihood<complexT>::clear()
{
	_decayAmps.resize(extents[0][0][0]);
}


// depends on naming convention for waves!!!
// VR_IGJPCMEIso....
template<typename complexT>
int
pwaLikelihood<complexT>::getReflectivity(const TString& waveName)
{
	int refl = 0;
	unsigned int reflIndex = 6;  // position of reflectivity in wave
	// check whether it is parameter or wave name
	if (waveName[0] == 'V')
		reflIndex = 9;
	if (waveName[reflIndex] == '-')
		refl= -1;
	else if (waveName[reflIndex] == '+')
		refl= +1;
	else {
		printErr << "cannot parse parameter/wave name '" << waveName << "'. "
		         << "cannot not determine reflectivity. Aborting..." << endl;
		throw;
	}
	if (_debug)
		printDebug << "extracted reflectivity = " << refl << " from parameter name "
		           << "'" << waveName << "' (char position " << reflIndex << ")" << endl;
	return refl;
}


// copy values from array that corresponds to the function parameters
// to structure that corresponds to the complex production amplitudes
// taking into account rank restrictions
template<typename complexT>
void
pwaLikelihood<complexT>::copyFromParArray
(const double*  inPar,             // input parameter array
 ampsArrayType& outVal,            // array of complex output values [rank][reflectivity][wave index]
 value_type&    outFlatVal) const  // output value corresponding to flat wave
{
	// group parameters by rank and wave index only; keeps the order defined in wave list
	outVal.resize(extents[_rank][2][max(_nmbWavesRefl[0], _nmbWavesRefl[1])]);
	unsigned int parIndex = 0;
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				value_type re, im;
				if (iWave < iRank)          // production amplitude is zero
					re = im = 0;
				else if (iWave == iRank) {  // production amplitude is real
					re = inPar[parIndex];
					im = 0;
					++parIndex;
				} else {                    // production amplitude is complex
					re = inPar[parIndex];
					im = inPar[parIndex + 1];
					parIndex += 2;
				}
				outVal[iRank][iRefl][iWave] = complexT(re, im);
			}
	outFlatVal = inPar[parIndex];
}


// copy values from structure that corresponds to complex
// production amplitudes to array that corresponds to function
// parameters taking into account rank restrictions
template<typename complexT>
void
pwaLikelihood<complexT>::copyToParArray
(const ampsArrayType& inVal,         // values corresponding to production amplitudes
 const value_type     inFlatVal,     // value corresponding to flat wave
 double*              outPar) const  // output parameter array
{
	for (unsigned int iRank = 0; iRank < _rank; ++iRank)
		for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
			for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave) {
				bt::tuple<int, int> parIndices = _prodAmpToFuncParMap[iRank][iRefl][iWave];
				if (bt::get<0>(parIndices) >= 0)  // real part
					outPar[bt::get<0>(parIndices)] = inVal[iRank][iRefl][iWave].real();
				if (bt::get<1>(parIndices) >= 0)  // imaginary part
					outPar[bt::get<1>(parIndices)] = inVal[iRank][iRefl][iWave].imag();
			}
	outPar[_nmbPars - 1] = inFlatVal;
}


template<typename complexT>
ostream&
pwaLikelihood<complexT>::print(ostream& out) const
{
	out << "pwaLikelihood parameters:" << endl
	    << "number of events ........................ " << _nmbEvents         << endl
	    << "rank .................................... " << _rank              << endl
	    << "number of waves ......................... " << _nmbWaves          << endl
	    << "number of positive reflectivity waves ... " << _nmbWavesRefl[1]   << endl
	    << "number of negative reflectivity waves ... " << _nmbWavesRefl[0]   << endl
	    << "number of function parameters ........... " << _nmbPars           << endl
	    << "number of fixed function parameters ..... " << _nmbParsFixed      << endl
	    << "print debug messages .................... " << _debug             << endl
#ifdef USE_CUDA
	    << "use CUDA kernels ........................ " << _cudaEnabled       << endl
#endif
	    << "use normalized amplitudes ............... " << _useNormalizedAmps << endl
	    << "list of waves: " << endl;
	for (unsigned int iRefl = 0; iRefl < 2; ++iRefl)
		for (unsigned int iWave = 0; iWave < _nmbWavesRefl[iRefl]; ++iWave)
			out << "        [" << setw(2) << sign((int)iRefl * 2 - 1) << " " << setw(3) << iWave << "] "
			    << _waveNames[iRefl][iWave] << "    threshold = "
			    << _waveThresholds[iRefl][iWave] << " MeV/c^2" << endl;
	out << "list of function parameters: " << endl;
	for (unsigned int iPar = 0; iPar < _nmbPars; ++iPar)
		out << "        [" << setw(3) << iPar << "] " << _parNames[iPar] << "    "
		    << "threshold = " << _parThresholds[iPar] << " MeV/c^2" << endl;
	return out;
}


template<typename complexT>
ostream&
pwaLikelihood<complexT>::printFuncInfo(ostream& out) const
{
	const string funcNames[NMB_FUNCTIONCALLENUM] = {"FdF", "Gradient", "DoEval", "DoDerivative", "Hessian"};
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i)
		if (_funcCallInfo[i].nmbCalls > 0)
			out << "    " << _funcCallInfo[i].nmbCalls
			    << " calls to pwaLikelihood<complexT>::" << funcNames[i] << "()" << endl
			    << "    time spent in pwaLikelihood<complexT>::" << funcNames[i] << "(): " << endl
			    << "        total time ................. " << sum(_funcCallInfo[i].totalTime) << " sec" << endl
			    << "        return value calculation ... " << sum(_funcCallInfo[i].funcTime ) << " sec" << endl
			    << "        normalization .............. " << sum(_funcCallInfo[i].normTime ) << " sec" << endl;
	return out;
}


template<typename complexT>
void
pwaLikelihood<complexT>::resetFuncCallInfo() const
{
	for (unsigned int i = 0; i < NMB_FUNCTIONCALLENUM; ++i) {
		_funcCallInfo[i].nmbCalls  = 0;
		// boost accumulators do not have a reset function
		_funcCallInfo[i].funcTime  = typename functionCallInfo::timeAccType();
		_funcCallInfo[i].normTime  = typename functionCallInfo::timeAccType();
		_funcCallInfo[i].totalTime = typename functionCallInfo::timeAccType();
	}
}


// explicit specializations
template class rpwa::pwaLikelihood<complex<float > >;
template class rpwa::pwaLikelihood<complex<double> >;