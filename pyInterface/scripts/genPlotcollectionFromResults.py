#!/usr/bin/env python

import argparse
import atexit
import multiprocessing
import sys
import tempfile
import os
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

	parser.add_argument("fitResults", type=str, metavar="fitResults", nargs='+', help="fitResult to get the production amplitudes")
	parser.add_argument("-o", "--output", type=str, metavar="output", dest="output", default="plotcollection.root",
	                    help="path to output file (default: '%(default)s')")
	parser.add_argument("-d", "--description", type=str, metavar="description", dest="description", default="" ,
	                    help="Description stored with the plotcollection collection")
	parser.add_argument("--label", type=str, metavar="label", dest="label", default="" ,
	                    help="Label of the result. If not given, a label is generated automatically.")
	parser.add_argument("--parallel", action="store_true", help="Build multibins in parallel")

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

	if not args.parallel:
		plotcollection = pyRootPwa.plotcollection(args.fitResults, description=args.description, label=args.label if args.label else None)
		plotcollection.buildMultibinSummedPlots()
		statusOk = plotcollection.write(args.output)
	else:
		label = pyRootPwa.plotcollection.buildLabelFromHash(args.fitResults, args.description)
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
		pool = multiprocessing.Pool(min(len(args.fitResults), multiprocessing.cpu_count()*2))
		tmpFilenames = pool.map(_worker, args.fitResults)
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

		master.buildMultibinSummedPlots()

		statusOk = master.write(args.output)

	if statusOk:
		sys.exit(0)
	else:
		sys.exit(2)
