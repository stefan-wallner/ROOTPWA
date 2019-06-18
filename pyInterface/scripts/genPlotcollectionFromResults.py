#!/usr/bin/env python

import argparse
import atexit
import multiprocessing
import sys
import tempfile
import os
import glob
import pyRootPwa
import pyRootPwa.utils

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="""
	                                                Generate a plotcollection collection from the given
	                                                fit resutls and stores it to the given
	                                                output file.
	                                             """
	                                )

	parser.add_argument("fitResults", type=str, metavar="fitResults", nargs='+', help="Files with fit results or folder containing .root files with fit resuls.")
	parser.add_argument("-o", "--output", type=str, metavar="output", dest="output", default="plotcollection.root",
	                    help="path to output file (default: '%(default)s')")
	parser.add_argument("-d", "--description", type=str, metavar="description", dest="description", default="" ,
	                    help="Description stored with the plotcollection collection")
	parser.add_argument("--label", type=str, metavar="label", dest="label", default="" ,
	                    help="Label of the result. If not given, a label is generated automatically.")
	parser.add_argument("--parallel", action="store_true", help="Build multibins in parallel")
	parser.add_argument("--n-ref-waves", dest='nRefWaves', type=int, default=4, help="Number of reference waves for which phases are generated.")
	parser.add_argument("--no-totals", dest='buildAllTotals', action="store_false", default=True, help="Do not build wave-totals.")
	parser.add_argument("--clearNames", action="store_true", help="Clear prod.-amp and wave names from further results. (Ugly hack. Should be used with care)")

	args = parser.parse_args()

	if not args.description:
		pyRootPwa.utils.printErr("No description given")

	if os.path.exists(args.output):
		pyRootPwa.utils.printErr("The output file '{0}' exist.".format(args.output))
		yesNo = raw_input("Override? [y/N]:").lower()
		if yesNo == 'y' or yesNo == "yes":
			os.remove(args.output)
		else:
			sys.exit(1)

	fitResults = []
	for path in args.fitResults:
		if os.path.isfile(path):
			fitResults.append(path)
		elif os.path.isdir(path):
			fitResults += glob.glob(os.path.join(path, "*.root"))
		else:
			raise Exception("Cannot find '{0}'!".format(path))

	plotcollection = None

	if not args.parallel:
		plotcollection = pyRootPwa.plotcollection(fitResults, description=args.description, label=args.label if args.label else None,
		                                          clearNames=args.clearNames)
	else:
		label = pyRootPwa.plotcollection.buildLabelFromHash(fitResults, args.description)
		processes = []
		def _worker(fitResultFilename):
			filedescriptor, tmpFilename = tempfile.mkstemp("_genPlotcollectionFromResults.root")
			plotcollection = pyRootPwa.plotcollection([fitResultFilename], description=args.description, label=label)
			os.close(filedescriptor)
			os.remove(tmpFilename)
			statusOk = plotcollection.write(tmpFilename)
			if not statusOk:
				raise Exception("Cannot build plotcollection forom file '{0}'.".format(fitResultFilename))
			return tmpFilename
		pool = multiprocessing.Pool(min(len(fitResults), multiprocessing.cpu_count()*2))
		tmpFilenames = pool.map(_worker, fitResults)
		pool.terminate()
		def _cleanup():
			for filename in tmpFilenames:
				if os.path.exists(filename):
					os.remove(filename)
		atexit.register(_cleanup)

		plotcollections = []
		for tmpFilename in tmpFilenames:
			plotcollection = pyRootPwa.plotcollection()
			plotcollection.load(tmpFilename)
			plotcollections.append(plotcollection)

		master = plotcollections[0]
		for plotcollection in plotcollections[1:]:
			master.mergePlotsInto(plotcollection)
		plotcollection = master

	pyRootPwa.utils.printInfo("Building default plots")
	plotcollection.buildDefaultPlots(nRefWaves=args.nRefWaves, buildAllTotals=args.buildAllTotals)

	pyRootPwa.utils.printInfo("Building multibin summed plots")
	plotcollection.buildMultibinSummedPlots()

	statusOk = plotcollection.write(args.output)

	if statusOk:
		sys.exit(0)
	else:
		sys.exit(2)
