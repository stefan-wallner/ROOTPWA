//-----------------------------------------------------------
//
// Description:
//    (BW) Component of mass dependent fit
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITCOMPONENTS_HH
#define MASSDEPFITCOMPONENTS_HH

#include <map>

#include <boost/multi_array.hpp>

#include <Math/Interpolator.h>

namespace libconfig {
	class Setting;
}
namespace ROOT {
	namespace Math {
		class Minimizer;
	}
}

namespace rpwa {

	namespace massDepFit {

		class pwachannel {

		public:

			pwachannel(const std::string& waveName,
			           std::complex<double> coupling,
			           const std::vector<double>& massBinCenters,
			           const std::vector<double>& phaseSpace);
			pwachannel(const rpwa::massDepFit::pwachannel& ch);
			~pwachannel();
    
			const std::string& getWaveName() const { return _waveName; }

			std::complex<double> getCoupling() const { return _coupling; }
			void setCoupling(std::complex<double> coupling) { _coupling = coupling; }

			double getPhaseSpace(const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::complex<double> getCouplingPhaseSpace(const double mass,
			                                           const size_t idxMass = std::numeric_limits<size_t>::max()) const;
    
		private:

			std::string _waveName;

			std::complex<double> _coupling;

			std::vector<double> _massBinCenters;
			std::vector<double> _phaseSpace;
			ROOT::Math::Interpolator* _interpolator;

		};

		class component {

		public:

			component(const std::string& name,
			          const size_t nrParameters);
			virtual ~component() {};

			const std::string& getName() const { return _name; }

			virtual bool init(const libconfig::Setting* configComponent,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual bool update(const libconfig::Setting* configComponent,
			                    const ROOT::Math::Minimizer* minimizer,
			                    const bool debug) const;

			const size_t getNrChannels() const { return _channels.size(); }
			const std::vector<pwachannel>& getChannels() const { return _channels; }
			const pwachannel& getChannel(const size_t i) const { return _channels[i]; } 
			const std::string& getChannelWaveName(const size_t i) const { return _channels[i].getWaveName(); }

			void getCouplings(double* par) const;
			void setCouplings(const double* par);

			size_t getNrParameters() const { return _nrParameters; }
			virtual void getParameters(double* par) const;
			virtual void setParameters(const double* par);

			virtual double getParameter(const size_t idx) const;
			virtual bool getParameterFixed(const size_t idx) const;
			virtual double getParameterLimitLower(const size_t idx) const;
			virtual bool getParameterLimitedLower(const size_t idx) const;
			virtual double getParameterLimitUpper(const size_t idx) const;
			virtual bool getParameterLimitedUpper(const size_t idx) const;
			virtual const std::string& getParameterName(const size_t idx) const;
			virtual double getParameterStep(const size_t idx) const;

			virtual std::complex<double> val(const double m) const = 0;

			virtual std::ostream& print(std::ostream& out) const;

		private:

			const std::string _name;

			std::vector<pwachannel> _channels;

			const size_t _nrParameters;

		protected:

			std::vector<double> _parameters;
			std::vector<bool> _parametersFixed;
			std::vector<double> _parametersLimitLower;
			std::vector<bool> _parametersLimitedLower;
			std::vector<double> _parametersLimitUpper;
			std::vector<bool> _parametersLimitedUpper;
			std::vector<std::string> _parametersName;
			std::vector<double> _parametersStep;

		};

		class pwacomponent : public component {

		public:

			pwacomponent(const std::string& name);
    
			virtual bool init(const libconfig::Setting* configComponent,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const double m) const;

			std::ostream& print(std::ostream& out) const;

		private:

			bool _constWidth;

		};

		class pwabkg : public component {

		public:

			pwabkg(const std::string& name);
    
			virtual bool init(const libconfig::Setting* configComponent,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const double m) const;

		private:

			double _m1;
			double _m2;
    
		};

	} // end namespace massDepFit

} // end namespace rpwa


inline
double
rpwa::massDepFit::pwachannel::getPhaseSpace(const double mass,
                                            const size_t idxMass) const
{
	if(idxMass != std::numeric_limits<size_t>::max()) {
		return _phaseSpace[idxMass];
	}
  
	return _interpolator->Eval(mass);
}


inline
std::complex<double>
rpwa::massDepFit::pwachannel::getCouplingPhaseSpace(const double mass,
                                                    const size_t idxMass) const
{
	return _coupling * getPhaseSpace(mass, idxMass);
}


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::component& component)
{
	return component.print(out);
}


#endif // MASSDEPFITCOMPONENTS_HH