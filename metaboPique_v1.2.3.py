#!/usr/bin/env python3

import sys, math, os, csv, gc, re, warnings
from multiprocessing import Process, Queue

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

import numpy as np

from scipy.optimize import curve_fit, OptimizeWarning
from scipy.stats import norm
from scipy.stats import spearmanr

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

warnings.simplefilter("ignore", OptimizeWarning)

leVersion = "v1.2.3"

if len(sys.argv) > 1:
	allMSdataPaths = sys.argv[1:]
else:

	import tkinter
	from tkinter import filedialog

	root = tkinter.Tk()
	root.withdraw()  # use to hide tkinter window

	#filez = tkFileDialog.askopenfilenames(parent=root,title='Choose a file')
	#print root.tk.splitlist(filez)

	currdir = os.getcwd()
	tempdir = filedialog.askdirectory(parent=root, initialdir=currdir, title='MetaboPique %s: Select analysis directory' % leVersion)
	if len(tempdir) > 0:
		#print "You chose %s" % tempdir
		allMSdataPaths = [tempdir]

allFiles = []
for allMSdataPath in allMSdataPaths:
	if os.path.exists(allMSdataPath):
		allMSdata = next(iter(os.walk(allMSdataPath)))
		for shortFilename in allMSdata[2]:
			truPath = os.path.join(allMSdataPath, shortFilename)
			allFiles.append([truPath, shortFilename])
	else:
		print("Data path %s is incorrect!!" % allMSdataPath)
		sys.exit()

print("\n<--------------------------------------------------------- THE WINDOW SHOULD BE THIS WIDE FOR OPTIMAL RUNNING ----------------------------------------------------------------------->\n")

print("------ MetaboPique %s ------\n" % leVersion)

topfolderName = input(" > New folder name for file output? ")
if os.path.exists(topfolderName): print("Folder already exists! Files may be overwritten!")
else: os.makedirs(topfolderName)

parameterDict = {"Version": "MetaboPique %s" % leVersion}

SMPcount = input(" > Multiprocessing thread count [8 default]? ")
SMPcount = 8 if SMPcount == "" else int(SMPcount)

polaritySet = input(" > Ionization mode: (1) Positive (2) negative [1 default]? ")
polaritySet = "+" if ((polaritySet == '') | (polaritySet == '1')) else "-"

print("\n  >>> Chromatogram selection parameters <<<")
precursor_ppm = input(" > MS1 scanning precursor m/z tolerance (ppm) [20 default]? ")
precursor_ppm = 20.0 if precursor_ppm == "" else float(precursor_ppm)

gaussianErr = input(" > Maximum threshold for MS1 chromatographic peak Gaussian error [0.45 default]? ")
gaussianErr = 0.45 if gaussianErr == "" else float(gaussianErr)

minPoints = input(" > Minimum data points contained within calculated peak width [4 default]? ")
minPoints = 4 if minPoints == "" else int(minPoints)

#minIsoTrace = input(" > Filter chromatographic peaks by adduct trace presence [Y/n]? ")
#minIsoTrace = True if ((minIsoTrace == "") | (minIsoTrace.lower() == "y")) else False

parameterDict["polaritySet"] = polaritySet; parameterDict["precursor_ppm"] = precursor_ppm; parameterDict["gaussianErr"] = gaussianErr; parameterDict["minPoints"] = minPoints;# parameterDict["minIsoTrace"] = minIsoTrace

print("\n  >>> MS2 filtering parameters <<<")

maxPrecursorMZ = input(" > Maximum MS2 precursor m/z value [1200 default]? ")
maxPrecursorMZ = 1200.0 if maxPrecursorMZ == "" else float(maxPrecursorMZ)

minFragmented = input(" > Minimum precursor fragmentation [default 0.1]? ")
minFragmented = 0.1 if minFragmented == '' else float(minFragmented)

selectCount = input(" > Maximum MS2 acquisition retention time distance from peak center (standard deviations) [3 default]? ")
selectCount = 3 if selectCount == "" else int(selectCount)

ms2minAb = input(" > Minimum MS2 fragment peak abundance threshold for inclusion (as a function of the base peak) [0.005 default]? ")
ms2minAb = 0.005 if ms2minAb == "" else float(ms2minAb)

parameterDict["minFragmented"] = minFragmented; parameterDict["maxPrecursorMZ"] = maxPrecursorMZ; parameterDict["selectCount"] = selectCount; parameterDict["ms2minAb"] = ms2minAb

print("\n  >>> Chromatographic clustering parameters <<<")
maxRTwindow = input(" > Maximum retention time window for clustering (seconds) [default 60]? ")
maxRTwindow = 60 if maxRTwindow == "" else int(maxRTwindow)

simThresh = input(" > Cosine similarity threshold for MS2 clustering validation [default 0.7]? ")
simThresh = 0.7 if simThresh == "" else float(simThresh)

parameterDict["maxRTwindow"] = maxRTwindow; parameterDict["simThresh"] = simThresh

r = robjects.r
importr('mzR')
#r.sink("/dev/null")
#r.options(warn=-1)
slot = r.slot


theAdducts = dict([('H', 1.00728), ('K', 38.96316), ('Na', 22.98922), ('Li', 7.01546), ('NH4', 18.03383), ('ACN', 42.03382), ('IsoProp', 61.06534), ('DMSO', 78.01394), ('CH3OH', 34.04077), ('Ca', 19.47738),
				('Cl', 34.96940), ('H2O', 18.010565), ('FA', 46.00548), ('Br', 78.91889), ('HAc', 60.02113), ('TFA', 113.99287), ('M', 0.0), ('H2CO2', 46.00548), ('C2HF3O2', 113.99287)])

theUnit = re.compile(r'(\-|\+)?(\d)?([A-Z][a-zA-Z0-9_]*)')  # this deconstructs the adduct annotation

# positive mode adducts
posAdducts = [(1, '[M+H]1+'), (1, '[M+Na]1+'), (1, '[M+NH4]1+'), (1, '[M+2Na-H]1+'), (1, '[2M+H]1+'), (1, '[2M+Na]1+'), (1, '[2M+NH4]1+'), (1, '[M+K]1+'),
				(2, '[M+2H]2+'), (2, '[M+H+Na]2+'), (2, '[M+H+K]2+'), (2, '[M+H+NH4]2+'), (2, '[M+2Na]2+'), (2, '[M+2NH4]2+'),
				(1, '[M+Li]1+'), (1, '[M+2K-H]1+'), (3, '[M+H+2Na]3+'), (3, '[M+2H+Na]3+')]

# negative mode adducts
negAdducts = [(-1, '[M-H]1-'), (-1, '[M+Na-2H]1-'), (-1, '[M+Cl]1-'), (-1, '[M+K-2H]1-'), (-1, '[2M-H]1-'), (-1, '[2M-H]1-'),
				(-2, '[M-2H]2-'), (-2, '[M+H2CO2-2H]2-'),
				(-1, '[M+H2CO2-H]1-'), (-1, '[M+C2HF3O2-H]1-'), (-1, '[2M+H2CO2-H]1-'), (-1, '[2M+C2HF3O2-H]1-'), (-1, '[3M-H]1-')]

'''
testSetAdducts = [(1, 1.00728, 1, '[M+H]1+'), (1, 22.98922, 1, '[M+Na]1+'), (1, 18.03383, 1, '[M+NH4]1+'), (1, 38.96316, 1, '[M+K]1+'),
					(2, 2*1.00728, 1, '[M+2H]2+'), (2, 1.00728+22.98922, 1, '[M+H+Na]2+'), (2, 1.00728+38.96316, 1, '[M+H+K]2+'),
					(1, 1.00728, 2, '[2M+H]1+'), (1, 22.98922, 2, '[2M+Na]1+'), (1, 38.96316, 2, '[2M+K]1+')]
'''
'''
testSetAdducts = [(1, 1.00728, 1, '[M+H]1+'), (1, 22.98922, 1, '[M+Na]1+'), (1, 18.03383, 1, '[M+NH4]1+'), (1, 38.96316, 1, '[M+K]1+'), (1, 44.971160, 1, '[M+2Na-H]1+'), (1, 42.033823, 1, '[M+ACN+H]1+'),
					(2, 2*1.00728, 1, '[M+2H]2+'), (2, 1.00728+22.98922, 1, '[M+H+Na]2+'), (2, 1.00728+38.96316, 1, '[M+H+K]2+'), (2, 1.00728+18.03383, 1, '[M+H+NH4]2+'),
					(1, 1.00728, 2, '[2M+H]1+'), (1, 22.98922, 2, '[2M+Na]1+'), (1, 38.96316, 2, '[2M+K]1+'), (1, 18.03383, 2, '[2M+NH4]1+')]
'''


postSetAdducts = [(1, 1.00728, 1, '[M+H]1+'), (1, 22.98922, 1, '[M+Na]1+'), (1, 18.03383, 1, '[M+NH4]1+'), (1, 38.96316, 1, '[M+K]1+'), (1, -17.0032881, 1, '[M-H2O+H]1+'),
					(2, 2*1.00728, 1, '[M+2H]2+'), (2, 1.00728+22.98922, 1, '[M+H+Na]2+'), (2, 1.00728+38.96316, 1, '[M+H+K]2+'),
					(3, 3*1.00728, 1, '[M+3H]3+'), (3, 2*1.00728+22.98922, 1, '[M+2H+Na]3+'),
					(1, 1.00728, 2, '[2M+H]1+'), (1, 22.98922, 2, '[2M+Na]1+'), (1, 38.96316, 2, '[2M+K]1+')]

negSetAdducts = [(1, -1.00728, 1, '[M-H]1-'), (1, -19.01839, 1, '[M-H2O-H]1-'), (1, 20.974666, 1, '[M+Na-2H]1-'), (1, 34.969402, 1, '[M+Cl]1-'), (1, 36.948606, 1, '[M+K-2H]1-'),
					(2, -2*1.00728, 1, '[M-2H]2-'),
					(3, -3*1.00728, 1, '[M-3H]3-'),
					(1, -1.00728, 2, '[2M-H]1-')]


if polaritySet == "+":
	testSetAdducts = postSetAdducts
else:
	testSetAdducts = negSetAdducts

massesToSearch = dict()

for aAdduct in testSetAdducts:
	aAdductName = aAdduct[3]

	for bAdduct in testSetAdducts:
		if aAdduct != bAdduct:
			bAdductName = bAdduct[3]

			pairName = (aAdductName, bAdductName)
			#				   multiply chg, subtract adduct, divide M | multiply M, add adduct, divide chg
			orderOfOperations = [aAdduct[0], aAdduct[1], aAdduct[2], bAdduct[2], bAdduct[1], bAdduct[0]]

			massesToSearch[pairName] = list(orderOfOperations)



def gauss(x, *p):
	A, mu, sigma = p
	return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def gaussSimp(x, *p):
	global isomean, isosd
	A = p[0]
	return A*np.exp(-(x-isomean)**2/(2.*isosd**2))


def ms1tracer(MS1data, precursorMZofInterest, precursorRIofInterest, theI, sortedAdductKeys, adductRangeSet, storedData, direction):

	bridgeTooFarLimit = 1000

	maxScan = max([scanNum for scanNum in MS1data])

	movingIndex = 1
	retentionTime = precursorRIofInterest
	threeStrikes = 0
	while (threeStrikes < 3):

		if theI+direction*movingIndex in MS1data:
			retentionTime = MS1data[theI+direction*movingIndex][2]
			ms1peakData = MS1data[theI+direction*movingIndex][0]
			theSize = MS1data[theI+direction*movingIndex][1]

			checkDivisor = 90.0 if theSize > 500 else 40.0

			theKset = dict()

			#----------- this is for all the adduct MZs + the precursor ------------------------------
			previousK = ""
			for adductMass in sortedAdductKeys:
				minThresholdMZ = (1-30.0/10**6)*adductMass
				maxThresholdMZ = (30.0/10**6+1)*adductMass

				for numerator in range(int(checkDivisor)-1, -1, -1):
					theIndex = int(round(theSize*(numerator/checkDivisor)))
					theMZval = ms1peakData[theIndex]


					if minThresholdMZ > theMZval:

						alreadyInRange = False

						if previousK != "":
							searchMinMZ = theKset[previousK][0]
							searchMaxMZ = theKset[previousK][1]

							if (theMZval >= searchMinMZ) & (theMZval <= searchMaxMZ):  # it is within the bounds of an existing search range!
								if maxThresholdMZ > searchMaxMZ:
									theKset[previousK][1] = maxThresholdMZ

								alreadyInRange = True


						if alreadyInRange is False:
							if theIndex in theKset: print("WTFWTF!!!")
							theKset[theIndex] = [theMZval, maxThresholdMZ]
							previousK = theIndex

						break

			#print theKset

			AbTotal = 0
			adductAbTotal = dict([(adductPrecMass, 0) for adductPrecMass in sortedAdductKeys])

			for theK in theKset:
				mzVal = ms1peakData[theK]
				Ab = ms1peakData[theK+theSize]

				maxMZcutoff = theKset[theK][1]
				while (mzVal <= maxMZcutoff):


					if Ab > 0:
						adductPreCheck = round(mzVal*10)

						if adductPreCheck in adductRangeSet:
							for rangeSet in adductRangeSet[adductPreCheck]:
								minThresholdMZ = rangeSet[0][0]
								maxThresholdMZ = rangeSet[0][1]

								if (mzVal > minThresholdMZ) & (mzVal < maxThresholdMZ):
									isoPrecMass = rangeSet[1]

									if isoPrecMass == precursorMZofInterest:
										AbTotal += Ab
									else:
										adductAbTotal[isoPrecMass] += Ab

					theK += 1

					if theK >= theSize:
						break

					#try:
					mzVal = ms1peakData[theK]
					Ab = ms1peakData[theK+theSize]
					#except IndexError:
					#	print "theK: %s theSize: %s maxScan: %s mzVal: %s" % (theK, theSize, maxScan, mzVal)


			if AbTotal > 0:
				storedData.append([AbTotal, retentionTime, adductAbTotal.copy()])
				threeStrikes = 0
			else:
				threeStrikes += 1

		trueI = theI+direction*movingIndex
		if ((trueI > 1) & (trueI < maxScan) & (abs(retentionTime-precursorRIofInterest) < bridgeTooFarLimit)):
			movingIndex += 1
		else:
			threeStrikes = 10



def goBackCheck(theI, precursorMZofInterest, precursorRIofInterest, MS1data, showGraph=False):
	global isomean, isosd

	#if precursorMZofInterest == 820.8720627:
	#	showGraph = True
	#isotopeSet = dict([(precursorMZofInterest+isoMass, isotopeList[isoMass]) for isoMass in isotopeList] + [(precursorMZofInterest-isoMass, isotopeList[isoMass]) for isoMass in isotopeList])


	adductSet = {precursorMZofInterest: ["precursor"]}
	for pairName in massesToSearch:
		massOperations = massesToSearch[pairName]

		predictedM = ((precursorMZofInterest*massOperations[0]) - massOperations[1])/massOperations[2]
		potentialMZ = (predictedM*massOperations[3] + massOperations[4])/massOperations[5]

		if potentialMZ > 12:
			potentialMZ_5place = round(potentialMZ, 5)
			if potentialMZ_5place not in adductSet:
				adductSet[potentialMZ_5place] = [pairName]
			else:
				adductSet[potentialMZ_5place].append(pairName)


	#isotopeSet = dict([(precursorMZofInterest+isoMass, isotopeList[isoMass]) for isoMass in isotopeList])

	#isoFastLookup = dict([(math.floor(isoPrecMass*10), isoPrecMass) for isoPrecMass in isotopeSet])

	adductRangeSet = dict()
	for adductMass in adductSet:
		minThresholdMZ = (1-30.0/10**6)*adductMass
		maxThresholdMZ = (30.0/10**6+1)*adductMass

		#precheckKey = round(minThresholdMZ*10)

		loPrecheck = round(minThresholdMZ*10)
		hiPrecheck = round(maxThresholdMZ*10)

		if loPrecheck in adductRangeSet:
			adductRangeSet[loPrecheck].add(((minThresholdMZ, maxThresholdMZ), adductMass))
		else:
			adductRangeSet[loPrecheck] = set([((minThresholdMZ, maxThresholdMZ), adductMass)])

		if hiPrecheck in adductRangeSet:
			adductRangeSet[hiPrecheck].add(((minThresholdMZ, maxThresholdMZ), adductMass))
		else:
			adductRangeSet[hiPrecheck] = set([((minThresholdMZ, maxThresholdMZ), adductMass)])



	storedData = []

	sortedAdductKeys = sorted(list(adductSet))

	ms1tracer(MS1data, precursorMZofInterest, precursorRIofInterest, theI, sortedAdductKeys, adductRangeSet, storedData, 1)
	ms1tracer(MS1data, precursorMZofInterest, precursorRIofInterest, theI, sortedAdductKeys, adductRangeSet, storedData, -1)

	del adductSet[precursorMZofInterest]

	if showGraph: print("Checking for MZ: %s RI: %s Scan: %s %s" % (precursorMZofInterest, precursorRIofInterest, theI, len(storedData)))

	theAbs = []; outputAbs = []
	theRTs = []

	validMS1Status = [[], []]

	if len(storedData) >= 2:
		minAb = min(storedData, key=lambda psig: psig[0])[0]
		maxAbRT = max(storedData, key=lambda psig: psig[0])[1]

		for data in storedData:
			theAbs.append((data[0]/minAb)-1); outputAbs.append(data[0])
			theRTs.append(data[1])

		RTrange = max(theRTs) - min(theRTs)

		validMS1Status = [theRTs, outputAbs]

	if len(storedData) >= 8:  # a minimum trace length!!!

		#print "RT range: %s %s" % (min(theRTs), max(theRTs))
		#print theAbs
		#print theRTs


		# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
		p0 = [1., maxAbRT, 1.]
		#if showGraph: print " >Checking for MZ: %s RI: %s Scan: %s" % (precursorMZofInterest, precursorRIofInterest, theI)
		try:
			#if showGraph: print " >> Checking for MZ: %s RI: %s Scan: %s" % (precursorMZofInterest, precursorRIofInterest, theI)
			coeff, var_matrix = curve_fit(gauss, theRTs, theAbs, p0=p0)

			# Get the fitted curve
			hist_fit = gauss(theRTs, *coeff)
			outputFit = [(val+1)*minAb for val in hist_fit]

			peakCenter = coeff[1]
			peakSD = abs(coeff[2])

			peakStart = peakCenter - peakSD*3
			peakEnd = peakCenter + peakSD*3


			#originatingPeakData = dict()
			#origTotAb_forRatio = 0

			sumOfFit = 0; sumOfOriginal = 0; sumOfErrors = 0
			idealizedNumPoints = 0

			for k in range(len(theRTs)):
				currentRT = theRTs[k]
				if ((currentRT > peakStart) & (currentRT < peakEnd)):
					fittedPoint = hist_fit[k]
					originalPoint = theAbs[k]

					sumOfFit += fittedPoint
					sumOfOriginal += originalPoint

					#origTotAb_forRatio += outputAbs[k]
					#origTotAb_forRatio += originalPoint

					sumOfErrors += abs(fittedPoint - originalPoint)

					idealizedNumPoints += 1

					#originatingPeakData[currentRT] = originalPoint
					#if ((currentRT > peakCenter - peakSD*2) & (currentRT < peakCenter + peakSD*2)):


			idealizedPeakWidth = peakEnd - peakStart

			if sumOfFit > 0:

				if (sumOfErrors/sumOfFit < gaussianErr) & (idealizedPeakWidth < RTrange*1.5) & (idealizedNumPoints >= minPoints):

					adductOutputData = dict()

					isomean = peakCenter; isosd = peakSD

					if showGraph:
						fig = plt.figure()
						ax = fig.add_subplot(111, projection='3d')
						mzRepeat = [precursorMZofInterest for num in theRTs]

						#ax.scatter(theRTs, mzRepeat, outputAbs, c='blue')
						#ax.scatter(theRTs, mzRepeat, outputFit, c='r')

						ax.scatter(theRTs, mzRepeat, theAbs, c='blue')
						ax.scatter(theRTs, mzRepeat, hist_fit, c='r')


						# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
						print('Fitted mean = ', coeff[1])
						print('Fitted standard deviation = ', coeff[2])

					adductScoreCard = dict(); adductScoreCard_better = dict()
					status = ["UNASSIGNED", "NOTRACE", None]

					for adductMZ in adductSet:

						adductAbs = []; adductRTs = []
						adductPoints = 0
						adductPeakDataPaired = []; adductTotAb_forRatio = 0

						for i, data in enumerate(storedData):
							isoRT = data[1]
							isoAb = data[2][adductMZ]

							adductAbs.append(isoAb)

							if isoAb > 0:
								adductRTs.append(isoRT)

								if ((isoRT > peakStart) & (isoRT < peakEnd)):
									adductPoints += 1
									adductTotAb_forRatio += isoAb

									originalAb = theAbs[i]
									adductPeakDataPaired.append((isoAb, originalAb))


						if adductPoints >= idealizedNumPoints*0.75:  # must have 75% of the original number of idealized data points
							pairNameList = adductSet[adductMZ]

							for pairName in pairNameList:
								if pairName[0] in adductScoreCard:
									adductScoreCard[pairName[0]] += 1
								else:
									adductScoreCard[pairName[0]] = 1

							minadductAb = min([val for val in adductAbs if val > 0])
							modadductPeakData = [(abPair[0]/minadductAb)-1 for abPair in adductPeakDataPaired]

							spearR, thePval = spearmanr(modadductPeakData, [abPair[1] for abPair in adductPeakDataPaired])

							if ((thePval < 0.10) & (spearR > 0.5)):
								adductOutputData[adductMZ] = adductAbs

								for pairName in pairNameList:
									if pairName[0] in adductScoreCard_better:
										adductScoreCard_better[pairName[0]][0] += 1
										adductScoreCard_better[pairName[0]][1].add(adductMZ)
									else:
										adductScoreCard_better[pairName[0]] = [1, set([adductMZ])]

							if showGraph:
								print("the pair: %s m/z: %s Spearman: %s" % (str(pairNameList), adductMZ, spearR))

								modadductPeakDataSum = sum(modadductPeakData); origadductPeakDataSum = sum([abPair[1] for abPair in adductPeakDataPaired])
								adductAbRatio = origadductPeakDataSum/modadductPeakDataSum
								modadductAbs_toPlot = [adductAbRatio*((adductab/minadductAb)-1) for adductab in adductAbs if adductab > 0]
								adductmzRepeat = [adductMZ for num in adductRTs]

								if spearR > 0.5:
									ax.scatter(adductRTs, adductmzRepeat, modadductAbs_toPlot, c="g")
								else:
									ax.scatter(adductRTs, adductmzRepeat, modadductAbs_toPlot, c="yellow")


					if len(adductScoreCard_better) > 0:
						clusteredScoreCard = dict()
						for adductType in adductScoreCard_better:
							traceCount = adductScoreCard_better[adductType][0]

							if traceCount in clusteredScoreCard:
								clusteredScoreCard[traceCount].append(adductType)
							else:
								clusteredScoreCard[traceCount] = [adductType]

						if len(adductScoreCard_better) == 1:
							status[0] = "ASSIGNED"
							status[1] = "UNANIMOUS"

						elif len(adductScoreCard_better) > 1:
							maxValue = max(list(clusteredScoreCard))

							if len(clusteredScoreCard[maxValue]) == 1:
								status[0] = "ASSIGNED"
								status[1] = "MAJORITY"
							else:
								allAdductTraceSets = []
								for adductType in clusteredScoreCard[maxValue]:
									adductTraceMZset = adductScoreCard_better[adductType][1]
									allAdductTraceSets.append(adductTraceMZset)

								traceBool = all(x == allAdductTraceSets[0] for x in allAdductTraceSets)

								if traceBool:
									status[0] = "ASSIGNED"
									status[1] = "AMBIGUOUS"
								else:
									status[0] = "UNASSIGNED"
									status[1] = "CONFLICTING"

						status[2] = adductScoreCard_better.copy()



					validMS1Status = [peakStart, peakEnd, peakCenter, [theRTs, outputAbs, outputFit, adductOutputData], peakSD, status]

					if showGraph:
						print(adductScoreCard)
						print(adductScoreCard_better)

						print("%s %s %s %s %s" % (sumOfFit, sumOfOriginal, sumOfErrors, sumOfErrors/sumOfFit, sumOfErrors/sumOfOriginal))
						plt.show()
				#elif (sumOfErrors/sumOfFit < gaussianErr):
				#else:
				#	plt.scatter(theRTs, outputAbs, label='Test data')
				#	plt.scatter(theRTs, hist_fit, color='r', label='Fitted data',)
				#	plt.show()
			else:
				validMS1Status = [theRTs, outputAbs]


		except RuntimeError:
			validMS1Status = [theRTs, outputAbs]
			''

			#if showGraph:
			#	print "Unable to fit Gaussian!"
			#	plt.scatter(theRTs, outputAbs, label='Test data')
			#	plt.show()


	if len(theRTs) > 0:
		minRT = min(theRTs)
		maxRT = max(theRTs)

	else:
		minRT = 0
		maxRT = 0
	#if showGraph: print " <<<< Checking for MZ: %s RI: %s Scan: %s" % (precursorMZofInterest, precursorRIofInterest, theI)
	return validMS1Status, minRT, maxRT


def firstPassMSdata(localMZfile, totalScans):
	ms2rawCount = 0
	ms2precursorDict = dict()
	MS2dataKeys = []

	MS1cache = dict()

	for i in range(1, totalScans + 1):
		scanData = r.header(localMZfile, i)

		msLevel = scanData.rx("msLevel")[0][0]

		#if i == 40000: break
		if int(msLevel) == 1:
			ms1peakData = r.peaks(localMZfile, i)
			theSize = scanData.rx("peaksCount")[0][0]
			retentionTime = scanData.rx("retentionTime")[0][0]

			MS1cache[i] = [ms1peakData, theSize, retentionTime]


		if int(msLevel) == 2:
			ms2rawCount += 1

			precursorMZ = scanData.rx("precursorMZ")[0][0]
			precursorCharge = int(scanData.rx("precursorCharge")[0][0])
			retentionTime = scanData.rx("retentionTime")[0][0]
			precIntensity = scanData.rx("precursorIntensity")[0][0]

			if ((precursorCharge >= 0) & (precursorCharge <= 10) & (precursorMZ < maxPrecursorMZ)):
				peakData = r.peaks(localMZfile, i)
				theSize = r.dim(peakData)[0]

				spectrumData = []
				maxAb = 0
				totalAb = 0; precursorAb = 0

				for j in range(0, theSize):
					mzVal = peakData[j]
					Ab = peakData[theSize+j]

					if Ab > 0:
						spectrumData.append([mzVal, Ab])

						totalAb += Ab

						if (mzVal < (precursorMZ + 1)):
							findThePrecursor_ppm = (abs(mzVal - precursorMZ) / precursorMZ) * 1000000

							if findThePrecursor_ppm < 20:
								precursorAb += Ab

						if Ab > maxAb:
							maxAb = Ab


				if totalAb > 0:
					if precursorAb / totalAb <= 1-minFragmented:

						correctedSpectrumData = [[peak[0], (peak[1] / maxAb) * 999.0] for peak in spectrumData if peak[1] / maxAb > ms2minAb]

						MS2dataKeys.append([retentionTime, precursorMZ, precursorCharge, i, precIntensity])

						ms2precursorDict[(precursorMZ, retentionTime)] = correctedSpectrumData



		if i % 2000 == 0:
			print(" %s scans analyzed, %s MS2 spectra detected (%s sufficiently fragmented), %s MS1 scans recorded" % (i, ms2rawCount, len(ms2precursorDict), len(MS1cache)))

	print(">> %s scans analyzed, %s MS2 spectra detected (%s sufficiently fragmented), %s MS1 scans recorded <<" % (i, ms2rawCount, len(ms2precursorDict), len(MS1cache)))

	firstBin_maxThresholdMZ = (precursor_ppm/10**6+1)*MS2dataKeys[0][1]
	firstBin_minThresholdMZ = (1-precursor_ppm/10**6)*MS2dataKeys[0][1]

	theBins = [[[firstBin_minThresholdMZ, firstBin_maxThresholdMZ], [MS2dataKeys[0]]]]

	for key in MS2dataKeys[1:]:
		precursorMZ = key[1]
		noMatch = True

		#if precursorMZ == 820.8720627:
		#	print "Stuff %s" % (len(theBins))


		for aBin in theBins:
			binMinMZ = aBin[0][0]
			binMaxMZ = aBin[0][1]

			if ((precursorMZ > binMinMZ) & (precursorMZ < binMaxMZ)):
				aBin[1].append(key)
				noMatch = False
				#if precursorMZ == 820.8721823:
				#	print "Stuff 2 %s" % (str(theBins))
				newMaxThresholdMZ = (precursor_ppm/10**6+1)*precursorMZ

				if newMaxThresholdMZ > binMaxMZ:
					aBin[0][1] = newMaxThresholdMZ

				break

		if noMatch is True:
			minThresholdMZ = (1-precursor_ppm/10**6)*precursorMZ
			maxThresholdMZ = (precursor_ppm/10**6+1)*precursorMZ

			theBins.append([[minThresholdMZ, maxThresholdMZ], [key]])

	#print "Number of bins %s" % (len(theBins))


	#print theBins
	return ms2precursorDict, len(ms2precursorDict), MS2dataKeys, theBins, MS1cache


def smp_secondPassMS1verifier(MS1data, binnedBin, dataQu, progressQu):

	for theBin in binnedBin:
		alreadyScanned = []
		_totalSpectra = 0; _validSpectra = 0; _nonvalidSpectra = 0
		_nonvalidGroups = 0; _validGroups = 0

		# adduct dtuff
		_assignStatus = {"UNASSIGNED": 0, "ASSIGNED": 0}
		_reasonStatus = {"NOTRACE": 0, "UNANIMOUS": 0, "MAJORITY": 0, "CONFLICTING": 0, "AMBIGUOUS": 0}

		for MS2key in theBin[1]:
			_totalSpectra += 1

			precIntensity = MS2key[4]
			scanNumber = MS2key[3]
			precursorCharge = MS2key[2]
			precursorMZ = MS2key[1]
			retentionTime = MS2key[0]

			previouslyScannedBool = False

			for roster in alreadyScanned:
				ppmErr = abs(precursorMZ-roster[0])/(precursorMZ)*1000000

				validChrom = roster[4]
				if len(validChrom) > 2:
					mingroupRT = min(roster[1], validChrom[0])-1
					maxgroupRT = max(roster[2], validChrom[1])+1
				else:
					mingroupRT = roster[1] - 2
					maxgroupRT = roster[2] + 2

				if ppmErr < precursor_ppm:
					if ((retentionTime > mingroupRT) & (retentionTime < maxgroupRT)):
						previouslyScannedBool = True
						#print "Hey that worked! %s %s %s" % (precursorMZ, mingroupRT, maxgroupRT)

						if len(validChrom) > 2:

							peakSD = validChrom[4]
							peakCenter = validChrom[2]

							selectHi = peakCenter + peakSD * selectCount
							selectLo = peakCenter - peakSD * selectCount

							if ((retentionTime > selectLo) & (retentionTime < selectHi)):
								_validSpectra += 1
								roster[5].append((precursorMZ, retentionTime, precursorCharge, precIntensity, scanNumber))
						else:
							if len(roster) == 6:
								_nonvalidSpectra += 1
								roster[5].append((precursorMZ, retentionTime, precursorCharge, precIntensity, scanNumber))



			if previouslyScannedBool is False:


				validChrom, minRT, maxRT = goBackCheck(scanNumber, precursorMZ, retentionTime, MS1data)

				_groupData = [precursorMZ, minRT, maxRT, precursorCharge, validChrom]

				if len(validChrom) > 2:
					_validGroups += 1

					peakSD = validChrom[4]
					peakCenter = validChrom[2]

					selectHi = peakCenter + peakSD * selectCount
					selectLo = peakCenter - peakSD * selectCount

					if ((retentionTime > selectLo) & (retentionTime < selectHi)):
						_validSpectra += 1
						_groupData.append([(precursorMZ, retentionTime, precursorCharge, precIntensity, scanNumber)])
					else:
						_groupData.append([])


					# adduct prediction stuff
					adductStatusInfo = validChrom[5]

					_assignStatus[adductStatusInfo[0]] += 1
					_reasonStatus[adductStatusInfo[1]] += 1

				else:
					if len(validChrom[0]) > 0:
						_nonvalidGroups += 1
						_nonvalidSpectra += 1

						_groupData.append([(precursorMZ, retentionTime, precursorCharge, precIntensity, scanNumber)])




				alreadyScanned.append(_groupData)
				#print "%s %s %s" % (precursorMZ, minRT, maxRT)

		for roster in alreadyScanned:

			validChrom = roster[4]

			if ((len(validChrom) == 2) & (len(roster) == 6)):
				MS2Keys = roster[5]
				allTheRTs = [MS2scan[1] for MS2scan in MS2Keys]

				fakePeakCenter = np.mean(allTheRTs)

				roster[4].append(fakePeakCenter)


		dataQu.put(alreadyScanned)
		progressQu.put([(_validSpectra, _nonvalidSpectra), _totalSpectra, _validGroups, _nonvalidGroups, _assignStatus, _reasonStatus])


def analyzeTheFile(mzFile, totalScans, ms2precursorGroupDict, ms2precursorDict_batch, theFileName):
	ms2precursorDict, ms2DataCount, MS2dataKeys, binnedKeys, MS1cache = firstPassMSdata(mzFile, totalScans)

	groupsByMS1 = []
	validGroups = 0; nonvalidGroups = 0
	validSpectra = 0; nonvalidSpectra = 0
	totalGroups = 0; totalSpectra = 0


	adductAssignData = {"UNASSIGNED": 0, "ASSIGNED": 0}
	adductReasonData = {"NOTRACE": 0, "UNANIMOUS": 0, "MAJORITY": 0, "CONFLICTING": 0, "AMBIGUOUS": 0}

	print("\nScanning associated MS1 data for MS2 spectrum qualification and grouping...")


	binnedBins = [[] for i in range(SMPcount)]

	for i, aBin in enumerate(binnedKeys):
		binnedBins[i % SMPcount].append(aBin)


	dataQu = Queue()
	progressQu = Queue()

	smp_processes = [Process(target=smp_secondPassMS1verifier, args=(MS1cache, binnedBins[i], dataQu, progressQu)) for i in range(SMPcount)]


	[proc.start() for proc in smp_processes]

	itsAlive = any([proc.is_alive() for proc in smp_processes])

	while(itsAlive | (not dataQu.empty()) | (not progressQu.empty())):
		if(dataQu.empty() is not True):
			resultsData = dataQu.get()

			groupsByMS1.extend(resultsData)


		if(progressQu.empty() is not True):
			progressNum = progressQu.get()

			validSpectra += progressNum[0][0]
			nonvalidSpectra += progressNum[0][1]

			validGroups += progressNum[2]
			nonvalidGroups += progressNum[3]

			totalSpectra += progressNum[1]

			# adduct stuff
			_adductAssign = progressNum[4]
			_adductReason = progressNum[5]

			for _id in _adductAssign: adductAssignData[_id] += _adductAssign[_id]
			for _id in _adductReason: adductReasonData[_id] += _adductReason[_id]


			sys.stdout.write("\r%s spectra analyzed - %s gauss groups (%s spec) %s non-gauss (%s spec) >Unassig: %s (NoTrace: %s Confl: %s) | Assig: %s (Unanim: %s Majority: %s Ambig: %s)<" % (totalSpectra, validGroups, validSpectra, nonvalidGroups, nonvalidSpectra,
				adductAssignData["UNASSIGNED"], adductReasonData["NOTRACE"], adductReasonData["CONFLICTING"],
				adductAssignData["ASSIGNED"], adductReasonData["UNANIMOUS"], adductReasonData["MAJORITY"], adductReasonData["AMBIGUOUS"]))
			sys.stdout.flush()


		itsAlive = any([proc.is_alive() for proc in smp_processes])

	#end
	[proc.join() for proc in smp_processes]
	dataQu.close()
	progressQu.close()



	print("\nTotal raw MS2 spectra: %s" % len(ms2precursorDict))


	filteredms2precursorDict = dict()
	numTrueGroups = 0
	for group in groupsByMS1:
		chargeState = group[3]
		validChrom = group[4]

		if len(group) == 6:
			if len(group[5]) > 0:
				#chromatographicData = validChrom[3]
				#print group
				numTrueGroups += 1

				avgPrecMZ = np.mean([spectrumKey[0] for spectrumKey in group[5]])

				ms2precursorGroupDict.append([(avgPrecMZ, validChrom[2], theFileName), validChrom, group[5]])

				for specKey in group[5]:
					#_chg = specKey[2]
					_mz = specKey[0]
					_rt = specKey[1]
					_spectrum = list(ms2precursorDict[(_mz, _rt)])

					filteredms2precursorDict[(_mz, _rt)] = _spectrum



	print("Final chromatographic groups: %s" % (numTrueGroups))

	ms2precursorDict_batch[theFileName] = filteredms2precursorDict

	del ms2precursorDict

# --------- composite spectrum methods

def localizedKDE(data, x, div):
	pdfVal = 0
	for val in data:
		theH = ((precursor_ppm/2.0)/1e6)*val
		pdfVal += norm.pdf(x, val, theH)

	pdfVal /= float(div)

	return pdfVal


def compositeMaker(superSpectrum, relevantSpecCount, specDict):
	sortedMZs = []
	supSpecData = []

	subSpectrum = sorted(superSpectrum, key=lambda psig: psig[0][0], reverse=False)

	for val in subSpectrum:
		sortedMZs.append(val[0][0])
		supSpecData.append([val[0][0], val[0][1], val[1]])  # mz, abundance, specID

	#print sorted(theMZs)

	theMargin = 0.05
	#sortedMZs = sorted(theMZs)
	theRanges = [[sortedMZs[0]-theMargin, sortedMZs[0]+theMargin, [sortedMZs[0]]]]
	currentRangeIndex = 0
	for mz in sortedMZs[1:]:
		if (mz-theMargin) <= theRanges[currentRangeIndex][1]:  # if the lower range of the next mz value is within the range of the current set...
			theRanges[currentRangeIndex][1] = mz + theMargin   # ...then just extend it!
			theRanges[currentRangeIndex][2].append(mz)
		else:
			currentRangeIndex += 1
			theRanges.append([mz-theMargin, mz+theMargin, [mz]])

	trueCenters = []
	theStep = 0.00001
	theMZcount = len(sortedMZs)

	for currentRange in theRanges:
		toEvaluate = np.arange(currentRange[0], currentRange[1], theStep)

		allVals = localizedKDE(currentRange[2], toEvaluate, theMZcount)


		prevPoint = allVals[0]
		currentPoint = allVals[1]

		localizedHalfMin = ((1/float(theMZcount))*norm.pdf(currentRange[1], currentRange[1], ((precursor_ppm/2.0)/1e6)*currentRange[1]))/2

		for i in range(1, len(allVals)-1):
			nextPoint = allVals[i+1]

			if currentPoint > localizedHalfMin:

				if (currentPoint > prevPoint) & (currentPoint > nextPoint):
					#print "%s %s" % (toEvaluate[i], currentPoint)
					#qu.put([indexVal, ])
					#nestedCenters[indexVal].append([round(toEvaluate[i], 4), []])
					trueCenters.append([round(toEvaluate[i], 4), []])

			#print "%s %s" % (i, currentPoint)
			prevPoint = currentPoint
			currentPoint = nextPoint

	#finalSpectrum = dict()
	correlatedSpectrum = []
	F_key = []
	Fcount = 0

	for trueMZdata in trueCenters:
		trueMZval = trueMZdata[0]
		expMZval = []

		compoundToFragDict = dict()

		ppmTrigger = True; daTrigger = True

		if len(supSpecData) == 0: daTrigger = False

		toPop = []
		#for frag in supSpecData:
		index = 0
		while (ppmTrigger & daTrigger):
			try:
				frag = supSpecData[index]
			except IndexError:
				print("OVER INDEX!!!!!! %s" % index)

			ppm_err = ((trueMZval - frag[0])/trueMZval)*1e6
			#print "%s %s" % (trueMZval, frag[0])

			if abs(ppm_err) < precursor_ppm:
				toPop.append(frag)
				trueMZdata[1].append(frag[1])
				expMZval.append(frag[0])
				compoundToFragDict[frag[2]] = frag[1]  # specID --> abundance

			elif ppm_err < 0: ppmTrigger = False

			index += 1
			if ((frag[0] - trueMZval) > 1) | (index >= len(supSpecData)): daTrigger = False

		[supSpecData.remove(taken) for taken in toPop]


		#classNum = len(trueMZdata[1])
		if len(trueMZdata[1]) > 0:
			cpdMZval = round(np.mean(expMZval), 4)
			#cpdAbun = max(trueMZdata[1])
			theError = ((max(expMZval) - min(expMZval))/cpdMZval)*1000000


			fragSpectrum = [0 for i in range(0, relevantSpecCount)]

			for specID in compoundToFragDict:
				#print specID
				#print cpdDict
				fragSpectrum[specDict[specID]] = compoundToFragDict[specID]

			correlatedSpectrum.append(fragSpectrum)
			F_key.append([cpdMZval, Fcount, theError])  # the actual m/z value and then the index
			Fcount += 1
			'''
			if classNum not in finalSpectrum:
				finalSpectrum[classNum] = [[cpdMZval, cpdAbun]]
			else:
				finalSpectrum[classNum].append([cpdMZval, cpdAbun])
			'''

	#print supSpecData
	#if len(supSpecData) > 0:
	#	for aSpec in supSpecData:
	#		print "  %s" % aSpec

	return correlatedSpectrum, F_key


def compositeSpecWorker(toComposite):
	groupFolderName = toComposite[0]
	superSpectrum = toComposite[1]
	specCount = toComposite[2]
	specDict = toComposite[3]
	textInfo = toComposite[4]

	lowestMonoIso = textInfo[0]
	predictedChg = textInfo[1]
	specCount = textInfo[2]

	correlatedSpectrum, F_key = compositeMaker(superSpectrum, specCount, specDict)

	majspecWriter = open(os.path.join(groupFolderName, "majoritySpectrum.txt"), 'w')
	majspecWriter.write("Precursor m/z: %s\nCharge: %s\nConstituent spectra: %s\n" % (lowestMonoIso, predictedChg, specCount))

	majoritySpectrum = []
	for item in F_key:
		index = item[1]
		theF = item[0]
		theError = item[2]

		nonZeroAbs = [ab for ab in correlatedSpectrum[index] if ab != 0]
		occurence = len(nonZeroAbs)

		if occurence/float(specCount) >= 0.5:
			arithMeanAb = np.mean(nonZeroAbs)
			majoritySpectrum.append([theF, arithMeanAb, theError, occurence])

	if len(majoritySpectrum) > 0:
		basePeakAb = max(majoritySpectrum, key=lambda psig: psig[1])[1]

		finalMajoritySpec = []

		for peak in majoritySpectrum:
			reformattedAb = round((peak[1]/basePeakAb)*999.0, 2)
			majspecWriter.write("%s %s (%s/%s) %sppm\n" % (peak[0], reformattedAb, peak[3], specCount, round(peak[2], 2)))

			finalMajoritySpec.append([peak[0], reformattedAb, peak[3], specCount, round(peak[2], 2)])

		majspecWriter.close()

		return finalMajoritySpec
	else:
		return False




def smp_theCompositer(subrawData, compositeQu):
	counter = 0

	for toComposite in subrawData:

		groupFolderName = toComposite[0]
		textInfo = toComposite[4]

		lowestMonoIso = textInfo[0]
		predictedAdduct = textInfo[1]
		specCount = textInfo[2]
		takeItAll = textInfo[3]
		allavgTrets = textInfo[4]
		avgTretText = str(round(np.mean(allavgTrets),2))

		if takeItAll == False:
			XICoutput = "Gaussian"
		else:
			XICoutput = "Non-gaussian"

		finalMajoritySpec = compositeSpecWorker(toComposite)

		if finalMajoritySpec is not False:
			counter += 1

			PEPMASS = round(lowestMonoIso, 4)
			compositeQu.put([groupFolderName.split("/")[-1], predictedAdduct, PEPMASS, specCount, XICoutput, avgTretText, finalMajoritySpec])


def theCompositer(rawData, mgfWriter, groupProgress, rawGroupProgress, majProgress, adductString):
	totalToMake = sum([len(parcel) for parcel in rawData])
	compositeQu = Queue()

	smp_processes = [Process(target=smp_theCompositer, args=(rawData[i], compositeQu)) for i in range(0, SMPcount)]

	[proc.start() for proc in smp_processes]

	progress = 0
	itsAlive = any([proc.is_alive() for proc in smp_processes])

	while(itsAlive | (not compositeQu.empty())):
		if(compositeQu.empty() is not True):
			theData = compositeQu.get()
			progress += 1

			sys.stdout.write("\r Progress: %s supergroups from %s raw groups [%s] --> %s/%s majority spectra " % (groupProgress, rawGroupProgress, adductString[1:], progress+majProgress, totalToMake+majProgress))
			sys.stdout.flush()



			#mgfWriter.write("\nID: %s\nPrecursor m/z: %s\nCharge: %s\nConstituent spectra: %s\n" % (theData[0], theData[1], theData[2], theData[3]))
			mgfWriter.write("BEGIN IONS\nTITLE=%s\nCHARGE=%s\nPEPMASS=%s\nCOMMENT=Type:MajoritySpectrum; ConstituentSpectra:%s; XIC:%s; avgRT:%s\n" % (theData[0], theData[1], theData[2], theData[3], theData[4], theData[5]))
			for peak in theData[6]:
				mgfWriter.write("%s %s (%s/%s) %sppm\n" % (peak[0], peak[1], peak[2], peak[3], peak[4]))
			mgfWriter.write("END IONS\n\n")

			#if progress % 150 == 0:
			#	gc.collect()

		itsAlive = any([proc.is_alive() for proc in smp_processes])

	[proc.join() for proc in smp_processes]
	compositeQu.close()

	#print "\n"

	return totalToMake


#------------------------------------------------------------------------------------------------------------------------------


ms2precursorGroupDict = []
ms2precursorDict_batch = dict()

diagWriter = open(os.path.join(topfolderName, "diagnosticData.txt"), 'w')
diagWriter.write("Processed chromatograms:\n")

print("\n%s files found:" % (len(allFiles)))
validCount = 0
for fileTruePath, anMZXMLfile in allFiles:
	if "mzxml" in anMZXMLfile.lower():
		print(" << %s" % anMZXMLfile)
		diagWriter.write(">>%s<<\n" % anMZXMLfile)

		validCount += 1

diagWriter.write("\n")
print("%s are mzXML files" % validCount)

if validCount == 0:
	print("No mzXML files found! Exiting!")
	sys.exit()


#negPatt = re.compile(r'[\._]neg[\._]')


for fileTruePath, anMZXMLfile in allFiles:
	if "mzxml" in anMZXMLfile.lower():
		#negSearch = negPatt.search(fileTruePath.lower())
		#if negSearch is not None:
		#	polaritySet = "-"

		mzFile = r.openMSfile(fileTruePath, backend="Ramp")
		basicData = r.runInfo(mzFile)
		totalScans = basicData.rx("scanCount")[0][0]

		print("\nReading raw MS file %s..." % (anMZXMLfile))

		analyzeTheFile(mzFile, totalScans, ms2precursorGroupDict, ms2precursorDict_batch, anMZXMLfile)
		gc.collect()


#for chargeState in ms2precursorGroupDict:
finalBatches = []


NmassMatches = {0.0: set([0]), 0.1: set([0])}
for _chg in [1, 2, 3]:
	for numC13 in [1, 2]:
		trueOffset = numC13*1.008664/_chg
		toTenth = round(trueOffset, 1)
		NmassMatcherLo = math.floor(trueOffset*10)/10
		NmassMatcherHi = math.ceil(trueOffset*10)/10
		if toTenth in NmassMatches:
			NmassMatches[toTenth].add(trueOffset)  # ; print "%s %s" % (1.008664/_chg,NmassMatcher)
		else:
			NmassMatches[toTenth] = set([trueOffset])
		toTenthHi = round(toTenth + 0.1, 1)
		if toTenthHi in NmassMatches:
			NmassMatches[toTenthHi].add(trueOffset)  # ; print "%s %s" % (1.008664/_chg,NmassMatcher)
		else:
			NmassMatches[toTenthHi] = set([trueOffset])
		toTenthLo = round(toTenth - 0.1, 1)
		if toTenthLo >= 0:
			if toTenthLo in NmassMatches:
				NmassMatches[toTenthLo].add(trueOffset)  # ; print "%s %s" % (1.008664/_chg,NmassMatcher)
			else:
				NmassMatches[toTenthLo] = set([trueOffset])




mgfWriter = open(os.path.join(topfolderName, "allCompiledSpectra.MGF"), 'w')

for param in parameterDict:
	diagWriter.write("%s: %s\n" % (param, parameterDict[param]))

sortedSpectrumData = sorted(ms2precursorGroupDict)
#print [spectrum[0] for spectrum in sortedSpectrumData]

groupCount = 0; rawGroupCount = 0

compositeSpectralRawData = [[] for i in range(SMPcount)]
compositeCount = 0; alreadyCalculated = 0

print("\n\nCommence final grouping for %s raw groups..." % (len(sortedSpectrumData)))

finalAdductAssignData = {"UNASSIGNED": 0, "ASSIGNED": 0}
finalAdductReasonData = {"NOTRACE": 0, "UNANIMOUS": 0, "MAJORITY": 0, "CONFLICTING": 0, "AMBIGUOUS": 0}

while len(sortedSpectrumData) > 0:

	theBatch = [sortedSpectrumData[0]]
	#print theBatch

	precursorMZ_max = sortedSpectrumData[0][0][0]

	precursorPeakCenterRT = sortedSpectrumData[0][0][1]

	precursorRT_min = precursorPeakCenterRT-maxRTwindow
	precursorRT_max = precursorPeakCenterRT+maxRTwindow
	batchRTs = [precursorPeakCenterRT]

	# you need to find the spectrum closest to the peak center, that will be the representative spectrum!
	spectrumKeys = sortedSpectrumData[0][2]
	closestSpecKey = min(spectrumKeys, key=lambda psig: abs(psig[1]-precursorPeakCenterRT))
	batchRefSpec = ms2precursorDict_batch[sortedSpectrumData[0][0][2]][(closestSpecKey[0], closestSpecKey[1])]

	currentKey = sortedSpectrumData[0][0]
	sortedSpectrumData.pop(0)

	refDenom = 0
	fastBatchSpectrum = dict()
	for peak in batchRefSpec:
		theW = (peak[0]**0.25)*math.sqrt(peak[1])
		refDenom += theW**2

		roundDownKey = math.floor(peak[0])
		if roundDownKey not in fastBatchSpectrum:
			fastBatchSpectrum[roundDownKey] = [(peak, theW)]
		else:
			fastBatchSpectrum[roundDownKey].append((peak, theW))

	refDenom = math.sqrt(refDenom)

	toPop = []
	for i, spectrum in enumerate(sortedSpectrumData):

		precursorMZ = spectrum[0][0]
		precursorRT = spectrum[0][1]

		if precursorMZ > precursorMZ_max+5:
			break

		offsetKey = round((precursorMZ-precursorMZ_max), 1)

		if offsetKey in NmassMatches:
			potentialPPMs = []
			for trueOffset in NmassMatches[offsetKey]:
				isoPrecursorMZ_max = precursorMZ_max+trueOffset
				_ppmdiff = abs(precursorMZ-isoPrecursorMZ_max)/isoPrecursorMZ_max*1000000
				potentialPPMs.append([_ppmdiff, trueOffset])

			relevantBool = True
			ppm_difference = min(potentialPPMs)


		else:
			relevantBool = False

		if relevantBool:
			if ppm_difference[0] < 20:
				if ((precursorRT > precursorRT_min) & (precursorRT < precursorRT_max)):
					# use a cosine similarity to do a quick check!
					_specKeys = spectrum[2]
					_closestSpecKey = min(_specKeys, key=lambda psig: abs(psig[1]-precursorRT))
					_theSpec = ms2precursorDict_batch[spectrum[0][2]][(_closestSpecKey[0], _closestSpecKey[1])]

					_numerator = 0
					_theDenom = 0
					for peak in _theSpec:
						_theSpecW = (peak[0]**0.25)*math.sqrt(peak[1])
						_theDenom += _theSpecW**2

						roundDownKey = math.floor(peak[0])
						if roundDownKey in fastBatchSpectrum:

							for refPeak in fastBatchSpectrum[roundDownKey]:
								refMZ = refPeak[0][0]
								spectrumMatchPPM = (abs(peak[0]-refMZ)/refMZ)*1000000

								if spectrumMatchPPM < 20:

									_numerator += _theSpecW*refPeak[1]

					_theDenom = math.sqrt(_theDenom)
					finalSimilarity = _numerator/(refDenom*_theDenom)

					#print "final: %s _denom: %s refDenom: %s numerator: %s" % (finalSimilarity, _theDenom, refDenom, _numerator)
					if finalSimilarity > simThresh:

						theBatch.append(spectrum)
						batchRTs.append(precursorRT)
						toPop.insert(0, i)

						precursorMZ_max = spectrum[0][0]
						precursorRT_min = max(batchRTs)-maxRTwindow
						precursorRT_max = min(batchRTs)+maxRTwindow




	finalBatches.append(theBatch)

	groupFolderName = os.path.join(topfolderName, str(round(theBatch[0][0][0], 4))+"_"+str(round(theBatch[0][0][1], 2)))
	if os.path.exists(groupFolderName):
		#print "Folder already exists!"
		diagWriter.write("Folder already exists!\n")

		counter = 2
		stillExists = True
		while stillExists:

			tempgroupFolderName = groupFolderName + "_" + str(counter)
			if not os.path.exists(tempgroupFolderName):
				stillExists = False
				groupFolderName = tempgroupFolderName
				os.makedirs(groupFolderName)
				break
			else:
				counter += 1



	else:
		os.makedirs(groupFolderName)

	# ---------composite spectrum code
	specCount = 0
	superSpectrum = []
	specDict = dict()  # this is really an indexer for the correlated composite spectrum

	superIndexCount = 0
	# ---------------------------------

	finalAdductScoreCard = dict()
	for aGroup in theBatch:
		chromKey = aGroup[0][2]
		chromData = aGroup[1]

		if len(chromData) > 3:
			_adductScoreCard = chromData[5][2]

			if _adductScoreCard is not None:
				for adductType in _adductScoreCard:
					labeledAdductMZs = set([(adductTraceMZ, chromKey) for adductTraceMZ in _adductScoreCard[adductType][1]])

					if adductType in finalAdductScoreCard:
						finalAdductScoreCard[adductType][0] += _adductScoreCard[adductType][0]
						finalAdductScoreCard[adductType][1].update(labeledAdductMZs)
					else:
						finalAdductScoreCard[adductType] = [_adductScoreCard[adductType][0], labeledAdductMZs]


	finalStatus = ["UNASSIGNED", "NOTRACE"]
	predictedAdduct = "UNKNOWN"

	if len(finalAdductScoreCard) > 0:
		finalTallyAdduct = dict()
		for adductType in finalAdductScoreCard:
			traceCount = finalAdductScoreCard[adductType][0]

			if traceCount in finalTallyAdduct:
				finalTallyAdduct[traceCount].append(adductType)
			else:
				finalTallyAdduct[traceCount] = [adductType]


		if len(finalAdductScoreCard) == 1:
			finalStatus[0] = "ASSIGNED"
			finalStatus[1] = "UNANIMOUS"

			predictedAdduct = finalTallyAdduct[list(finalTallyAdduct)[0]][0]

		elif len(finalAdductScoreCard) > 1:
			maxValue = max(list(finalTallyAdduct))

			if len(finalTallyAdduct[maxValue]) == 1:
				finalStatus[0] = "ASSIGNED"
				finalStatus[1] = "MAJORITY"

				predictedAdduct = finalTallyAdduct[maxValue][0]

			else:
				allAdductTraceSets = []
				for adductType in finalTallyAdduct[maxValue]:
					adductTraceMZset = finalAdductScoreCard[adductType][1]
					allAdductTraceSets.append(adductTraceMZset)

				traceBool = all(x == allAdductTraceSets[0] for x in allAdductTraceSets)

				if traceBool:
					predictedAdduct = str(finalTallyAdduct[maxValue])

					finalStatus[0] = "ASSIGNED"
					finalStatus[1] = "AMBIGUOUS"
				else:
					finalStatus[0] = "UNASSIGNED"
					finalStatus[1] = "CONFLICTING"

	finalAdductAssignData[finalStatus[0]] += 1
	finalAdductReasonData[finalStatus[1]] += 1

	chromSet = set([aGroup[0][2] for aGroup in theBatch])

	gaussCount = 0
	nongaussCount = 0


	for aGroup in theBatch:
		chromData = aGroup[1]

		if len(chromData) > 3:
			gaussCount += 1
		else:
			nongaussCount += 1

	takeItAll = False  # generally you don't want to take all the spectra, only the stuff from gaussian XICs
	if gaussCount == 0:
		takeItAll = True  # ...that is unless there's NO gaussian XIC at all, so you need to take it all!

	theWriter = open(os.path.join(groupFolderName, "chromatographicData.csv"), 'w')
	theWriterCSV = csv.writer(theWriter, delimiter=',')
	theWriterCSV.writerow(["Predicted Adduct Type: %s\tPeak count: %s (%s gaussian, %s non-gaussian)\tChromatogram count: %s" % (predictedAdduct, len(theBatch), gaussCount, nongaussCount, len(chromSet))])
	theWriterCSV.writerow([])

	allavgMZs = []
	allavgTrets = []

	for aGroup in theBatch:
		#print "%s" % str(aGroup[0])
		diagWriter.write("%s\n" % str(aGroup[0]))

		chromKey = aGroup[0][2]
		chromData = aGroup[1]

		if len(chromData) > 3:
			gaussBool = True
			appendTxt = "g"
		else:
			gaussBool = False
			appendTxt = "n"

		spectrumKeys = aGroup[2]
		avgMZ = 0
		avgTret = 0


		# write each raw spectrum out
		for keys in spectrumKeys:
			_chg = keys[2]
			_mz = keys[0]
			_rt = keys[1]
			_precInt = keys[3]
			_precScanNo = keys[4]
			theSpectrum = ms2precursorDict_batch[chromKey][(_mz, _rt)]
			_mzText = str(round(_mz, 4))
			_rtText = str(_rt)

			specWriter = open(os.path.join(groupFolderName, _mzText + "_" + _rtText + "_" + appendTxt + ".txt"), 'w')


			adductType = "%s%s" % (str(_chg), polaritySet)
			specWriter.write("Precursor m/z: %s \nPrecursor Intensity: %s\nCharge: %s (%s predicted)\nRetention time: %s\nScan: %s\nNumPeaks: %s\nChromatogram: %s\n" % (_mzText, _precInt, adductType, predictedAdduct, _rtText, _precScanNo, str(len(theSpectrum)), chromKey))
			for peak in theSpectrum:
				_peakMZtext = str(round(peak[0], 4))
				_peakAbtext = str(peak[1])
				specWriter.write("%s %s\n" % (_peakMZtext, _peakAbtext))

				if (gaussBool | takeItAll): superSpectrum.append((peak, (_precScanNo, chromKey)))

			specWriter.close()

			# calculated for all XICs
			avgMZ += _mz
			avgTret += _rt

			if (gaussBool | takeItAll):
				# ---------composite spectrum code
				specDict[(_precScanNo, chromKey)] = superIndexCount
				superIndexCount += 1
				# ---------------------------------

				lastValidSpec = [chromKey, (_mz, _rt), _precScanNo, _precInt, _rtText]  # only used by singleton spectra in MGF file

		# calculated for all XICs
		avgMZ /= len(spectrumKeys)
		avgTret /= len(spectrumKeys)

		if (gaussBool | takeItAll):
			specCount += len(spectrumKeys)
			lastValidChromStatus = gaussBool  # really only used for singleton spectra to output gaussianness in MGF file

			# only recorded for 'valid' XICs!
			allavgMZs.append(avgMZ)
			allavgTrets.append(avgTret)

		if gaussBool:
			theWriterCSV.writerow(["Precursor m/z: %s PeakCenter: %s PeakStart: %s PeakEnd: %s Chromatogram: %s" % (round(avgMZ, 4), chromData[2], chromData[0], chromData[1], chromKey)])

			isoData = chromData[3][3]
			isotopeList = sorted(list(isoData))

			theWriterCSV.writerow(["Retention Time", "Raw MS1 Scan Abundance", "Fitted Gaussian"]+isotopeList)
			for i, val in enumerate(chromData[3][0]):
				isotopeReadout = [isoData[anIso][i] for anIso in isotopeList]
				theWriterCSV.writerow([val, chromData[3][1][i], chromData[3][2][i]]+isotopeReadout)
			theWriterCSV.writerow([])

		else:
			theWriterCSV.writerow(["Precursor m/z: %s XICcenter: %s Chromatogram: %s" % (round(avgMZ, 4), chromData[2], chromKey)])

			theWriterCSV.writerow(["Retention Time", "Raw MS1 Scan Abundance"])

			for _rt, _int in zip(chromData[0], chromData[1]):
				theWriterCSV.writerow([_rt, _int])

			theWriterCSV.writerow([])

	adductString = ">Unassig: %s (NoTrace: %s Conflict: %s) | Assig: %s (Unanim: %s Majority: %s Ambig: %s)" % (finalAdductAssignData["UNASSIGNED"], finalAdductReasonData["NOTRACE"], finalAdductReasonData["CONFLICTING"],
																														finalAdductAssignData["ASSIGNED"], finalAdductReasonData["UNANIMOUS"], finalAdductReasonData["MAJORITY"], finalAdductReasonData["AMBIGUOUS"])
	#print adductString
	diagWriter.write(adductString + "\n")


	# ---------composite spectrum code
	if specCount > 1:
		textInfo = [min(allavgMZs), predictedAdduct, specCount, takeItAll, allavgTrets]

		compositeSpectralRawData[compositeCount % SMPcount].append([groupFolderName, superSpectrum, specCount, specDict, textInfo])
		compositeCount += 1

		if compositeCount % 1000 == 0:
			calced = theCompositer(compositeSpectralRawData, mgfWriter, groupCount, rawGroupCount, alreadyCalculated, adductString)
			alreadyCalculated += calced

			compositeSpectralRawData = [[] for i in range(SMPcount)]

	else:
		#mgfWriter.write("\nID: %s\nPrecursor m/z: %s\nCharge: %s\nConstituent spectra: %s\n" % (groupFolderName.split("/")[-1], avgMZ, predictedAdduct, specCount))
		chromKey = lastValidSpec[0]
		theMZ = lastValidSpec[1][0]
		scanNo = lastValidSpec[2]
		precInt = lastValidSpec[3]
		rtText = lastValidSpec[4]
		validSpectrum = ms2precursorDict_batch[chromKey][lastValidSpec[1]]

		if lastValidChromStatus:
			XICoutput = "Gaussian"
		else:
			XICoutput = "Non-gaussian"

		PEPMASS = round(theMZ, 4)
		mgfWriter.write("BEGIN IONS\nTITLE=%s\nCHARGE=%s\nPEPMASS=%s\nCOMMENT=Type:SingletonSpectrum; SourceChromatogram:%s; Scan:%s; PrecIntensity:%s; XIC:%s; RT:%s\n" % (groupFolderName.split("/")[-1], predictedAdduct,
																																										PEPMASS, chromKey, scanNo, precInt, XICoutput, rtText))
		for peak in theSpectrum:
			_peakMZtext = str(round(peak[0], 4))
			_peakAbtext = str(round(peak[1], 2))
			mgfWriter.write("%s %s\n" % (_peakMZtext, _peakAbtext))
		mgfWriter.write("END IONS\n\n")


	groupCount += 1; rawGroupCount += len(theBatch)
	#sys.stdout.write("\r     Progress:  %s groups created                                                                 " % (groupCount))
	sys.stdout.write("\r Progress: %s supergroups from %s raw groups [%s] -> %s/%s majority spectra " % (groupCount, rawGroupCount, adductString[4:], alreadyCalculated, compositeCount))
	sys.stdout.flush()

	# ---------------------------------

	theWriter.close()

	for popit in toPop:
		sortedSpectrumData.pop(popit)

calced = theCompositer(compositeSpectralRawData, mgfWriter, groupCount, rawGroupCount, alreadyCalculated, adductString)
alreadyCalculated += calced

diagWriter.close()
mgfWriter.close()

print("\nComplete!")
