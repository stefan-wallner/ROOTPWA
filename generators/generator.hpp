#ifndef GENERATOR_HH_
#define GENERATOR_HH_

#include "generatorParameters.hpp"
#include "generatorPickerFunctions.h"
#include "primaryVertexGen.h"

namespace rpwa {

	class generator {

	  public:

		generator()
			: _pickerFunction(NULL),
			  _primaryVertexGen(NULL) { }

		virtual ~generator() {
			delete _pickerFunction;
			delete _primaryVertexGen;
		};

		virtual void setBeam(const rpwa::Beam& beam) { _beam = beam; };
		virtual void setTarget(const rpwa::Target& target) { _target = target; };
		virtual void setTPrimeAndMassPicker(const rpwa::massAndTPrimePicker& pickerFunction) {
			(*_pickerFunction) = pickerFunction;
		}
		virtual void setPrimaryVertexGenerator(rpwa::primaryVertexGen* primaryVertexGen) {
			_primaryVertexGen = primaryVertexGen;
		}


	  protected:

		rpwa::Beam _beam;
		rpwa::Target _target;
		rpwa::massAndTPrimePicker* _pickerFunction;
		rpwa::primaryVertexGen* _primaryVertexGen;

	};

}

#endif
