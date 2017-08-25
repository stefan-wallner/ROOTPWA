#include "ampIntegralMatrixMetadata_py.h"

#include <boost/python.hpp>

#include <TDirectory.h>

#include "ampIntegralMatrix.h"
#include "ampIntegralMatrixMetadata.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {
	int ampIntegralMatrixMetadata_write(rpwa::ampIntegralMatrixMetadata& self, const char* name = 0) {
		return self.Write(name);
	}

	bp::list ampIntegralMatrixMetadata_getKeyFileContents(const rpwa::ampIntegralMatrixMetadata& self) {
		return bp::list(self.getKeyFileContents());
	}

	bp::list ampIntegralMatrixMetadata_getAmplitudeHashes(const rpwa::ampIntegralMatrixMetadata& self) {
		return bp::list(self.getAmplitudeHashes());
	}

	bp::dict ampIntegralMatrixMetadata_multibinBoundaries(const rpwa::ampIntegralMatrixMetadata& self) {
		bp::dict retval;
		rpwa::multibinBoundariesType multibinBoundaries = self.multibinBoundaries();
		for (rpwa::multibinBoundariesType::const_iterator iter = multibinBoundaries.begin(); iter != multibinBoundaries.end(); ++iter) {
			retval[iter->first] = bp::make_tuple(iter->second.first, iter->second.second);
		}
		return retval;
	}

	void ampIntegralMatrixMetadata_setMultibinBoundaries(rpwa::ampIntegralMatrixMetadata& self, bp::dict pyMultibinBoundaries) {
		rpwa::multibinBoundariesType multibinBoundaries;
		const bp::list keys = pyMultibinBoundaries.keys();
		for (int i = 0; i < bp::len(keys); ++i) {
			std::string binningVar = bp::extract<std::string>(keys[i]);
			double lowerBound      = bp::extract<double>(pyMultibinBoundaries[binningVar][0]);
			double upperBound      = bp::extract<double>(pyMultibinBoundaries[binningVar][1]);
			multibinBoundaries[binningVar] = rpwa::boundaryType(lowerBound, upperBound);
		}
		self.setMultibinBoundaries(multibinBoundaries);
	}

	const rpwa::ampIntegralMatrixMetadata* ampIntegralMatrixeMetadata_readIntegralFile(PyObject* pyInputFile,
	                                                                                   const bool& quiet = false)
	{
		TFile* inputFile = rpwa::py::convertFromPy<TFile*>(pyInputFile);
		if(not inputFile) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for inputFile when executing rpwa::ampIntegralMatrixMetadata::readIntegralFile()");
			bp::throw_error_already_set();
		}
		return rpwa::ampIntegralMatrixMetadata::readIntegralFile(inputFile, quiet);
	}

	bool ampIntegralMatrixMetadata_setAmpIntegralMatrix_1(rpwa::ampIntegralMatrixMetadata& self, PyObject* integralMatrixPy) {
		rpwa::ampIntegralMatrix* integralMatrix = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix*>(integralMatrixPy);
		if(not integralMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for ampIntegralMatrix when executing rpwa::ampIntegralMatrixMetadata::setAmpIntegralMatrix()");
			bp::throw_error_already_set();
		}
		return self.setAmpIntegralMatrix(integralMatrix);
	}

	bool ampIntegralMatrixMetadata_setAmpIntegralMatrix_2(rpwa::ampIntegralMatrixMetadata& self, rpwa::ampIntegralMatrix* integralMatrix) {
		return self.setAmpIntegralMatrix(integralMatrix);
	}

	bool ampIntegralMatrixMetadata_writeToFile(rpwa::ampIntegralMatrixMetadata& self, PyObject* pyOutputFile) {
		TFile* outputFile = rpwa::py::convertFromPy<TFile*>(pyOutputFile);
		if(not outputFile) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for outputFile when executing rpwa::ampIntegralMatrixMetadata::writeToFile()");
			bp::throw_error_already_set();
		}
		return self.writeToFile(outputFile);
	}
}

void rpwa::py::exportAmpIntegralMatrixMetadata() {

	bp::class_<rpwa::ampIntegralMatrixMetadata>("ampIntegralMatrixMetadata")
		.def(bp::init<const  rpwa::ampIntegralMatrixMetadata&>())

		.def(bp::self_ns::str(bp::self))

		.def("Write", &::ampIntegralMatrixMetadata_write, bp::arg("name")=0)

		.def("getKeyFileContents", &::ampIntegralMatrixMetadata_getKeyFileContents)
		.def("getAmplitudeHashes", &::ampIntegralMatrixMetadata_getAmplitudeHashes)
		.def("readIntegralFile"
		     , &::ampIntegralMatrixeMetadata_readIntegralFile
		     , (bp::arg("inputFile"), bp::arg("quiet")=false)
		     , bp::return_value_policy<bp::reference_existing_object>()
		)
		.staticmethod("readIntegralFile")
		.def("getAmpIntegralMatrix", &rpwa::ampIntegralMatrixMetadata::getAmpIntegralMatrix, bp::return_value_policy<bp::reference_existing_object>())

		.def("writeToFile", &::ampIntegralMatrixMetadata_writeToFile, bp::arg("outputFile"))

		.def("setAmpIntegralMatrix", &::ampIntegralMatrixMetadata_setAmpIntegralMatrix_1, bp::with_custodian_and_ward<1,2>(), bp::arg("integralMatrix"))
		.def("setAmpIntegralMatrix", &::ampIntegralMatrixMetadata_setAmpIntegralMatrix_2, bp::with_custodian_and_ward<1,2>(), bp::arg("integralMatrix"))

		.def("contentHash", &rpwa::ampIntegralMatrixMetadata::contentHash, bp::return_value_policy<bp::return_by_value>())
		.def("rootpwaGitHash", &rpwa::ampIntegralMatrixMetadata::rootpwaGitHash, bp::return_value_policy<bp::return_by_value>())
		.def("mergeIntegralMatrix", &rpwa::ampIntegralMatrixMetadata::mergeIntegralMatrix, bp::arg("secondMatrix"))

		.def("multibinBoundaries", &ampIntegralMatrixMetadata_multibinBoundaries)
		.def("addEventMetadata", &ampIntegralMatrixMetadata::addEventMetadata, bp::arg("eventMetadata"))
		.def("addAmplitudeHash", &ampIntegralMatrixMetadata::addAmplitudeHash, bp::arg("amplitudehash"))
		.def("hasAmplitudeHash", &ampIntegralMatrixMetadata::hasAmplitudeHash, bp::arg("amplitudehash"))
		.def("setHash",  &ampIntegralMatrixMetadata::setHash)
		.def("recalculateHash", &rpwa::ampIntegralMatrixMetadata::recalculateHash)
		.def("setGitHash",  &ampIntegralMatrixMetadata::setGitHash)

		.def("setMultibinBoundaries", &::ampIntegralMatrixMetadata_setMultibinBoundaries)

		.def("addKeyFileContent" , &rpwa::ampIntegralMatrixMetadata::addKeyFileContent)
		.def("hasKeyFileContent" , &rpwa::ampIntegralMatrixMetadata::hasKeyFileContent)
		.def_readonly("integralObjectName", &rpwa::ampIntegralMatrixMetadata::objectNameInFile);
}
