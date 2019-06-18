#!/usr/bin/env python

import argparse
import sys
import os
import itertools
import subprocess

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="merge fit results"
	                                )
	parser.add_argument("-q", "--quiet", action="store_true")
	parser.add_argument("--same-multibin", dest="sameMultibin", action="store_true",
	                    help="The fit results, which should be merged, must be from the same multibin.")
	parser.add_argument("--not-strip", dest="strip", action="store_false",
	                    help="The integral and covariance matrices are striped from each fit result except the best one and the best converged one in each multibin")

	parser.add_argument("inputFileNames", type=str, metavar="inFileNames", nargs='*', help="path to input files")
	parser.add_argument("outputFileName", type=str, metavar="outFileName", help="path to output file")

	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	if os.path.exists(args.outputFileName):
		printErr("Output file '{0}' exists!".format(args.outputFileName))
		sys.exit(1)
	if not args.strip:
		subprocess.check_call("hadd '{o}' '{i}'".format(o = args.outputFileName, i="' '".join(args.inputFileNames)), shell=True)
		sys.exit(0)

	if not args.quiet:
		printInfo("Load fit results")
	fitResultsInBins = pyRootPwa.utils.getFitResultsFromFiles(args.inputFileNames, stripMatricesFromFurtherAttempts=args.strip, quiet=args.quiet)

	if args.sameMultibin and len(fitResultsInBins) > 1:
		printErr("Fund fit results in the following multibins:")
		for b in sorted(fitResultsInBins.keys()):
			printInfo("\t{0}".format(b))
		sys.exit(1)

	fitResults = itertools.chain.from_iterable([v for v in fitResultsInBins.values()])
	nResults = sum([len(res) for res in fitResultsInBins.values()])

	if os.path.exists(args.outputFileName): # check again as it takes some time to load the results
		printErr("Output file '{0}' exists!".format(args.outputFileName))
		sys.exit(1)
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "CREATE")
	if not outputFile: # Do this up here. Without the output file, nothing else makes sense
		pyRootPwa.utils.printErr("could not open output file. Aborting...")
		sys.exit(1)
	if not args.quiet:
		printInfo("writing result(s) to '" + args.outputFileName + "'")
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	fitResult = pyRootPwa.core.fitResult()
	tree = outputFile.Get(valTreeName)
	if not tree:
		if not args.quiet:
			printInfo("file '" + args.outputFileName + "' is empty. "
			        + "creating new tree '" + valTreeName + "' for PWA result.")
		tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
		if not fitResult.branch(tree, valBranchName):
			printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
			sys.exit(1)
	else:
		fitResult.setBranchAddress(tree, valBranchName)

	if not args.quiet:
		progressBar = pyRootPwa.utils.progressBar(0, nResults)
		progressBar.start()
	for i,result in enumerate(fitResults):
		fitResult.fill(result)
		tree.Fill()
		if not args.quiet:
			progressBar.update(i)

	nmbBytes = tree.Write()
	outputFile.Close()
	if nmbBytes == 0:
		printErr("problems writing fit result to TKey 'fitResult' "
		       + "in file '" + args.outputFileName + "'")
		sys.exit(1)
	else:
		if not args.quiet:
			printSucc("wrote fit result to TKey 'fitResult' "
			        + "in file '" + args.outputFileName + "'")
