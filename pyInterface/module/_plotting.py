import gc
import hashlib
import os
from collections import defaultdict
import pyRootPwa.core
import pyRootPwa.utils
from pyRootPwa.utils import printErr
ROOT = pyRootPwa.utils.ROOT



class plotcollection(object):
	def __init__(self, fitResultFilenames = None, description = "", label=None):
		self._multibinPlots = {}
		self._multibinSummedPlots = {} # indices: [<summing-variable>][<multibin>],
		                               # where <multibin> is a multibin in all variables except the summed one
		self._waveNames = []
		self._descriptions = {}
		self._labels = []

		if fitResultFilenames is not None:
			if label is None:
				label = plotcollection.buildLabelFromHash(fitResultFilenames, description)
			self._labels.append(label)
			if fitResultFilenames is not None:
				self._descriptions[self._labels[0]] = description
				self._initFromFitresultsFilenames(fitResultFilenames)


	def _initFromFitresultsFilenames(self, fitResultFilenames):
		'''
		Initialize plots from fit-result files
		'''
		# build bash
		self._initFromFitresults(pyRootPwa.utils.getFitResultsFromFiles(fitResultFilenames,stripMatricesFromFurtherAttempts=True))


	def _initFromFitresults(self, fitresults):
		'''
		Initialize from fitresults.
		@param fitresults: dictionary with multibin as key and list of attempts as value
		'''
		# get list of all multibins except binning in mass
		multibins = sorted(list(set( [b.getSubMultiBin(exception="mass") for b in fitresults.keys()] )))


		pyRootPwa.utils.printInfo("Building multibin plots")
		# first group fit results by multibin
		fitresultsInMultibin = defaultdict(list)
		for multibin in multibins:
			for multibinResults, results in fitresults.iteritems():
				if multibinResults in multibin:
					fitresultsInMultibin[multibin] += results
		del fitresults.values()[:]
		del fitresults
		gc.collect()

		for multibin in multibins:
			self.addMultiBin(multibin, fitresults=fitresultsInMultibin[multibin])
			del fitresults[multibin][:]
			del fitresults[multibin]
			gc.collect()

		pyRootPwa.utils.printInfo("Building default plots")
		for _,multibinPlots in self.iterMultibinPlots():
			multibinPlots.buildDefaultPlots()


	def buildMultibinSummedPlots(self):
		if not self._multibinPlots:
			return

		variables = self._multibinPlots.keys()[0].variables()
		for variable in variables:
			multibinsInOtherVariables = sorted(list(set(b.getSubMultiBin(exception=variable) for b in self.multibins())))
			for multibinInOtherVariables in multibinsInOtherVariables:
				multibinsToSumOver = [b for b in self.multibins() if b in multibinInOtherVariables]
				multibinPlotsToSumOver = [ self.multibinPlots(mb) for mb in multibinsToSumOver ]
				multibinSummedPlots = pyRootPwa.core.multibinPlots()
				if multibinSummedPlots.buildMultibinSummedPlots(multibinPlotsToSumOver):
					if variable not in self._multibinSummedPlots:
						self._multibinSummedPlots[variable] = {}
					self._multibinSummedPlots[variable][multibinInOtherVariables] = multibinSummedPlots


	def mergePlotsInto(self, other):
		'''
		Merges the plots of other into the current collection

		This collection will contain the multibins of both collection
		If a multibin exists in both collections, the multibin will be merged
		@param other: other plot collection that will be merged into this one
		@type other: plotcollection
		@return: True if merging was successful, False otherwise
		'''
		pyRootPwa.utils.printInfo("Merge plot collections")
		# merge multibins
		for multibin in other.multibins():
			if multibin in self._multibinPlots:
				statusOK = self.multibinPlots(multibin).mergePlotsInto(other.multibinPlots(multibin))
				if not statusOK:
					return False
			else:
				self.addMultiBin(multibin, multibinplots = other.multibinPlots(multibin))

		# merge multibin-summed plots
		for variable, multibin, multibinplots in other.iterMultibinSummedPlots():
			if variable in self._multibinSummedPlots:
				if multibin in self._multibinSummedPlots[variable]:
					statusOK = self.multibinSummedPlots(variable, multibin).mergePlotsInto(multibinplots)
					if not statusOK:
						return False
				else:
					self._multibinSummedPlots[variable][multibin] = multibinplots
			else:
				self._multibinSummedPlots[variable] = {multibin: multibinplots}

		for label in other.labels(): # ensures to keep the order
			if label not in self._labels:
				self._labels.append(label)
		self._descriptions.update(other.descriptions())
		self._waveNames = sorted(list(set(  self._waveNames + other.waveNames() )))
		return True

	def addMultiBin(self, multibin, fitresults = None, dirIn = None, multibinplots = None, onlyBest = False):
		'''
		Add a further multibin to this plot collection
		@param fitresults: if given, multibin is greated from fit results
		@param dirIn: if given, multibin is loaded from ROOT-file directory
		@param onlyBest: if dirIn given, load only the best results from dir
		@param multibinplots: if given, the multibinplots object will be added to this plot collection
		'''
		if multibin in self._multibinPlots:
			msg = "Try to add multibin '{0}' to plotcollection, which is already in the plot collection.".format(multibin)
			pyRootPwa.utils.printErr(msg)
			raise Exception(msg)

		if fitresults is not None:
			self._multibinPlots[multibin] = pyRootPwa.core.multibinPlots(fitresults, self._labels[0], self._descriptions[self._labels[0]])
		elif dirIn is not None:
			mbp = pyRootPwa.core.multibinPlots()
			if mbp.load(dirIn, onlyBest):
				self._multibinPlots[multibin] = mbp
			else:
				msg = "Cannot load multibinPlots from file"
				pyRootPwa.utils.printErr(msg)
				raise Exception
		elif multibinplots is not None:
			self._multibinPlots[multibin] = multibinplots

		self._waveNames = sorted(list(set(self._waveNames + self._multibinPlots[multibin].waveNames())))
		for label in self._multibinPlots[multibin].labels(): # ensures to keep the order
			if label not in self._labels:
				self._labels.append(label)
		# append descriptions
		for label, description in zip(self._multibinPlots[multibin].labels(), self._multibinPlots[multibin].descriptions()):
			if label not in self._descriptions:
				self._descriptions[label] = description
			elif self._descriptions[label] != description:
				pyRootPwa.utils.printErr("Found different descriptions for the same label '{0}': '{1}' != '{2}'".format(label,
				                                                                                                        description,
				                                                                                                        self._descriptions[label]))
				raise "Found different descriptions"


	def multibins(self):
		'''
		@return: sorted list of multibins
		'''
		return sorted(self._multibinPlots.keys())


	def waveNames(self):
		'''
		@return: sorted list of wave names of all multibins
		'''
		return self._waveNames


	def descriptions(self):
		'''
		@return: dict of <label>: <description>
		'''
		return self._descriptions


	def labels(self):
		'''
		@return: list of labels for all contained plots
		'''
		return self._labels


	def multibinPlots(self, multibin):
		'''
		@return multibinPlots object of the given multibin
		'''
		if isinstance(multibin, int):
			multibin = self.multibins()[multibin]
		return self._multibinPlots[multibin]


	def multibinSummedPlots(self, variable = None, multibinInOtherVariables = None):
		'''
		@param variable variable in which the summation was performed
		                if None, returen the multibin summed plot in all variables
		@param multibinInOtherVariables multibin in all variables except the summed variable
		                                if None and only one multibin in other variables, return the correspinding multibinPlot object
		@return multibinPlots object of the multibin-summed spectra, summed in the given variable
		'''
		if variable is None:
			if len(self._multibinSummedPlots) == 1:
				return self._multibinSummedPlots.values()[0].values()[0]
			printErr("Need to implement multibinSummed plot over all variables") # TODO
			return None
		if variable in self._multibinSummedPlots:
			if multibinInOtherVariables is None:
				if len(self._multibinSummedPlots[variable]) == 1:
					return self._multibinSummedPlots[variable].values()[0]
				printErr("More than one multibin in other variables than the summed variable. Define the multibin")
				return None
			else:
				if multibinInOtherVariables in self._multibinSummedPlots[variable]:
					return self._multibinSummedPlots[variable][multibinInOtherVariables]
				printErr("Multibin '{0}' not in list of multibins for multibin-summed plots in variable '{1}'".format(multibinInOtherVariables, variable))
				return None
		printErr("Variable '{0}' not in list of summed variables".format(variable))
		return None


	def iterMultibinPlots(self):
		'''
		Iterate over sorted multibins
		@return multibin, multibinPlots
		'''
		for multibin in self.multibins():
			yield multibin, self._multibinPlots[multibin]


	def iterMultibinSummedPlots(self):
		'''
		Iterate over sorted multibin-summed plots
		@return summed-variable, multibin in other variables, multibinPlots
		'''
		for variable in sorted(self._multibinSummedPlots.keys()):
			for multibin in sorted(self._multibinSummedPlots[variable].keys()):
				yield variable, multibin, self._multibinSummedPlots[variable][multibin]


	def write(self, filenameOrRootfile):
		fileOut = None
		if isinstance(filenameOrRootfile, ROOT.TFile):
			fileOut = filenameOrRootfile
		else:
			fileOut = ROOT.TFile(filenameOrRootfile, "CREATE")

		if not fileOut.IsOpen():
			pyRootPwa.utils.printErr("Cannot open output file")
			return False

		pyRootPwa.utils.printInfo("Write plot collection to file '{0}'".format(fileOut.GetName()))


		# multibin plots
		for multibin in self.multibins():
			dirName = "bin__" + multibin.uniqueStr()
			dirOut = fileOut.mkdir(dirName)
			self._multibinPlots[multibin].write(dirOut)

		# multibin-summed plots
		dirSummed = fileOut.mkdir("multibinSummed")
		for variable in sorted(self._multibinSummedPlots.keys()):
			dirSummedVariable = dirSummed.mkdir(variable)
			for multibin in sorted(self._multibinSummedPlots[variable].keys()):
				dirName = "bin__" + multibin.uniqueStr()
				dirOut = dirSummedVariable.mkdir(dirName)
				self._multibinSummedPlots[variable][multibin].write(dirOut)
		fileOut.Close()
		return True


	def load(self, filenameOrRootfile, onlyBest=False):
		'''
		@param onlyBest: If true, load only the best fit result in each bin
		@return True if loading was successful
		'''
		fileIn = None
		if isinstance(filenameOrRootfile, ROOT.TFile):
			fileIn = filenameOrRootfile
		else:
			filename = os.path.expandvars(os.path.expanduser(filenameOrRootfile))
			fileIn = ROOT.TFile(filename, "OPEN")

		if not fileIn or not fileIn.IsOpen():
			printErr("Cannot read from input file")
			return False

		pyRootPwa.utils.printInfo("Load plot collection from file '{0}'".format(fileIn.GetName()))

		for k in fileIn.GetListOfKeys():
			if k.GetName().startswith("bin__"): # load multibin plots
				multibin = pyRootPwa.utils.multiBin.fromUniqueStr(k.GetName()[5:])
				dirIn = k.ReadObj()
				self.addMultiBin(multibin, dirIn = dirIn, onlyBest=onlyBest)
			if k.GetName() == "multibinSummed": # load multibin-summed plots
				dirSummed = k.ReadObj()
				for kVariable in dirSummed.GetListOfKeys():
					variable = kVariable.GetName()
					dirVariable = kVariable.ReadObj()
					for kMultibin in dirVariable.GetListOfKeys():
						multibin = pyRootPwa.utils.multiBin.fromUniqueStr(kMultibin.GetName()[5:])
						dirIn = kMultibin.ReadObj()
						mbsp = pyRootPwa.core.multibinPlots()
						mbsp.load(dirIn, onlyBest)
						if variable not in self._multibinSummedPlots:
							self._multibinSummedPlots[variable] = {}
						self._multibinSummedPlots[variable][multibin] = mbsp
		fileIn.Close()
		return True

	@classmethod
	def loadFromFile(cls, filenameOrRootfile, onlyBest=False):
		'''
		@param onlyBest: If true, load only the best fit result in each bin
		@rtype: plotcollection
		'''
		collection = cls()
		collection.load(filenameOrRootfile, onlyBest=onlyBest)
		return collection

	@classmethod
	def buildLabelFromHash(cls, fitResultFilenames, description):
		label = hashlib.sha1()
		label.update(description)
		for filename in fitResultFilenames:
			label.update(os.path.realpath(filename))
			label.update(str(os.path.getmtime(filename)))
		return label.hexdigest()[0:7]
