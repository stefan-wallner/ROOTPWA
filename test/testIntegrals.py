#!/usr/bin/env python

import argparse
import sys
import os
import subprocess as sp
import glob
import numpy as np

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Test on the fly integration"
	                                )

	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	fileManager = pyRootPwa.loadFileManager("fileManager.pkl")
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	binId = 0
	nBunches = 3
	printInfo("Calculate integrals in bin {0} in {1} bunches".format(binId, nBunches))
	for i in range(nBunches):
		sp.check_call("${{ROOTPWA}}/build/bin/calcIntegrals -e generated -b {b} --on-the-fly --split-events-N {en} --split-events-bunch {eb}".format(
	                  b = binId, en = nBunches, eb = i), shell=True)

	printInfo("Merge integral files")
	origIntegralFile = fileManager.getIntegralFilePath(fileManager.binList[binId], pyRootPwa.core.eventMetadata.GENERATED)
	if not os.path.exists(origIntegralFile):
		printErr("Integral file '{0}' does not exist".format(origIntegralFile))
		sys.exit(1)
	mergedFilePath = origIntegralFile + ".merged"
	bunchFilePaths = glob.glob(origIntegralFile + ".???")
	sp.check_call("${{ROOTPWA}}/build/bin/mergeIntegrals --same-multibin '{i}' '{o}' ".format(
	               i = "' '".join(bunchFilePaths), o = mergedFilePath), shell=True)

	printInfo("Compare integral files")
	origFile = ROOT.TFile.Open(origIntegralFile, "READ")
	origMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(origFile)
	origMatrix = origMeta.getAmpIntegralMatrix()

	mergedFile = ROOT.TFile.Open(mergedFilePath, "READ")
	mergedMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(mergedFile)
	mergedMatrix = mergedMeta.getAmpIntegralMatrix()

	if origMatrix.nmbWaves() != mergedMatrix.nmbWaves():
		printErr("Different number of waves")
		sys.exit(100)

	if origMatrix.nmbEvents() != mergedMatrix.nmbEvents():
		printErr("Different number of events")
		sys.exit(101)

	for i in range(origMatrix.nmbWaves()):
		if origMatrix.waveName(i) != mergedMatrix.waveName(i):
			print "Different wave ordering / wave names"
			sys.exit(102)

	deviations = np.empty((origMatrix.nmbWaves(), origMatrix.nmbWaves(),2))
	for i in range(origMatrix.nmbWaves()):
		for j in range(origMatrix.nmbWaves()):
			deviations[i,j,0] = abs((mergedMatrix.element(i,j) - origMatrix.element(i,j)).real) / abs(origMatrix.element(i,j).real)
			if mergedMatrix.element(i,j).imag != 0.0:
				deviations[i,j,1] = abs((mergedMatrix.element(i,j) - origMatrix.element(i,j)).imag) / abs(origMatrix.element(i,j).imag)
			else:
				deviations[i,j,1] = abs((mergedMatrix.element(i,j) - origMatrix.element(i,j)).imag)
	deviations = deviations.reshape(deviations.shape[0]*deviations.shape[1]*2)
	if np.max(deviations) > 1e-12:
		printErr("Maximal deviation = {0}".format(np.max(deviations)))
		printErr("Mean deviation    = {0}".format(np.mean(deviations)))
		sys.exit(1)
	else:
		printSucc("Maximal deviation = {0}".format(np.max(deviations)))
		printSucc("Mean deviation    = {0}".format(np.mean(deviations)))
		sys.exit(0)




