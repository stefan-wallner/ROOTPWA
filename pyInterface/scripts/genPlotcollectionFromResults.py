#!/usr/bin/env python

import argparse
import sys
import os
import pyRootPwa
import pyRootPwa.utils

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="""
	                                                Generate a plot collection from the given
	                                                fit resutls and stores it to the given
	                                                output file.
	                                             """
	                                )

	parser.add_argument("fitResults", type=str, metavar="fitResults", nargs='+', help="fitResult to get the production amplitudes")
	parser.add_argument("-o", "--output", type=str, metavar="output", dest="output", default="plotcollection.root",
	                    help="path to output file (default: '%(default)s')")
	parser.add_argument("-d", "--description", type=str, metavar="description", dest="description", default="" ,
	                    help="Description stored with the plot collection")

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

	plotcollection = pyRootPwa.plotcollection(args.fitResults, description=args.description)
	statusOk = plotcollection.write(args.output)

	if statusOk:
		sys.exit(0)
	else:
		sys.exit(2)
