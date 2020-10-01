
import os, sys, re, msgpack, math, gc, subprocess, csv
from time import sleep
from shutil import copy
from multiprocessing import Process, Queue
#import numpy as np


precursorPatt = re.compile(r'Precursor\sm/z:\s(\d+\.\d+)')
spectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.\d+)')

#chromPatt = re.compile(r'(\d{1,7}\.\d{0,4}),(\d+\.\d+),(\d+\.\d+)')
#chromPatt = re.compile(r'(\d{1,7}\.\d{0,4}),\d+\.{0,1}\d*(?:,\d+\.{0,1}\d*)+')
AbPatt = re.compile(r',\d+\.{0,1}\d*')
isoMZPatt = re.compile(r'\d+\.\d+')

precursorMZPatt = re.compile(r'Precursor\sm/z:\s(\d+\.\d+)')

adductPatt = re.compile(r'Charge:\s(.+)')
adductPattSingleton = re.compile(r'Charge:\s.+\s\((.+)\spredicted\)')

constituentPatt = re.compile(r'Constituent\sspectra:\s(\d+)')
#oneChargePatt = re.compile(r'\d{1,2}')
peakcenterPatt = re.compile(r'PeakCenter:\s(\d+\.\d+)')

filenamePatt = re.compile(r'Chromatogram:\s(.+\..+)')
orderedFilePatt = re.compile(r'>>(.+\..+)<<')
polarityPatt = re.compile(r'polaritySet:\s(.)')


allMSdataPath = sys.argv[1]

justName = allMSdataPath
if allMSdataPath[-1] == os.sep:
	justName = allMSdataPath[:-1]

justName = justName.split(os.sep)[-1]


orderedFileList = []
with open(os.path.join(allMSdataPath, "diagnosticData.txt"), 'r') as diagReader:
	next(diagReader)

	for line in diagReader:
		if line.startswith(">>"):
			_chromName = re.search(orderedFilePatt, line).group(1)
			orderedFileList.append(_chromName)

		if line.startswith("polaritySet"):
			thePolarity = re.search(polarityPatt, line).group(1)
			thePolarity = "P" if thePolarity=="+" else "N"

		if line.startswith("simThresh"):
			break

print("Chromatograms analyzed in MetaboPique run:")

for _chromName in orderedFileList:
	print(" > %s" % (_chromName))



allMSdata = next(os.walk(allMSdataPath))
rootDirectory = allMSdata[0]

clusterCount = 0


specAnalysis = dict()

print("\nReading files...")
for clusterDir in allMSdata[1]:
	chromStatus = dict([(_chromName, "N\\A") for _chromName in orderedFileList])

	clusterPath = os.path.join(allMSdataPath, clusterDir)
	clustData = next(os.walk(os.path.join(allMSdataPath, clusterDir)))
	clusterCount += 1

	clustFiles = clustData[2]

	chromReader = open(os.path.join(clusterPath, "chromatographicData.csv"), "r")

	firstLine = next(chromReader).rstrip()

	for line in chromReader:
		if "enter: " in line:
			_initialInfo = line.rstrip()
			_precursorChrom = float(re.search(precursorPatt, _initialInfo).group(1))
			_identifier = re.search(filenamePatt, _initialInfo).group(1)

			if "XIC" in line:
				chromStatus[_identifier] = "Non-gaussian"
			else:
				chromStatus[_identifier] = "Gaussian"

	chromReader.close()

	exists = os.path.isfile(os.path.join(clusterPath, "majoritySpectrum.txt"))

	if exists:
		spectrumReader = open(os.path.join(clusterPath, "majoritySpectrum.txt"), "r")

		_lineCount = 0
		finalSpectrum = []
		for line in spectrumReader:
			if _lineCount == 0:
				precursorMZ = float(re.search(precursorMZPatt, line).group(1))
			if _lineCount == 1:
				predictedAdduct = re.search(adductPatt, line).group(1)
			if _lineCount == 2:
				constitCount = int(re.search(constituentPatt, line).group(1))
			if _lineCount > 2:
				spectrumData = re.search(spectrumPatt, line)

				if spectrumData is not None:
					finalSpectrum.append([float(spectrumData.group(1)), float(spectrumData.group(2))])

			_lineCount += 1

		srcFile = os.path.join(clusterPath, "majoritySpectrum.txt")

		specAnalysis[clusterDir] = [precursorMZ, predictedAdduct, finalSpectrum, srcFile, constitCount, chromStatus]

	else:
		#print("There are %s files in here!" % (len(clustFiles)))
		if len(clustFiles) == 2: lookforPatt = ".txt"
		else: lookforPatt = "_g.txt"

		for filename in clustFiles:


			if lookforPatt in filename:
				srcFile = os.path.join(clusterPath, filename)
				spectrumReader = open(srcFile, "r")

				precursorMZ = float(re.search(precursorPatt, next(spectrumReader)).group(1))	 # line 0

				_intensity = float(re.search(r'\d+\.\d+', next(spectrumReader)).group())	 # line 1

				line2 = next(spectrumReader)												 # line 2
				_charge = int(re.search(r'\d', line2).group())
				predictedAdduct = re.search(adductPattSingleton, line2).group(1)

				_rt = float(re.search(r'\d+\.\d+', next(spectrumReader)).group())			 # line 3
				_scanNo = float(re.search(r'\d+', next(spectrumReader)).group())			 # line 4
				next(spectrumReader)														 # line 5
				_id = re.search(filenamePatt, next(spectrumReader)).group(1)				 # line 6


				finalSpectrum = []
				for line in spectrumReader:

					spectrumData = re.match(spectrumPatt, line)

					if spectrumData is not None:
						finalSpectrum.append([float(spectrumData.group(1)), float(spectrumData.group(2))])

				specAnalysis[clusterDir] = [precursorMZ, predictedAdduct, finalSpectrum, srcFile, 1, chromStatus]
				#print(srcFile)

		if clusterDir not in specAnalysis:
			print(" >> Error detected in %s <<" % clusterDir)



	if clusterCount % 5 == 0:
		sys.stdout.write("\r   %s clusters" % (clusterCount))
		sys.stdout.flush()


print("\nComplete! %s clusters found" % (clusterCount))

searchType = input(" > Search type: (1) Direct (2) Hybrid (3) In-Source [1 default]? ")
searchType = 1 if searchType == "" else int(searchType)

if searchType not in [1,2,3]:
	sys.exit("Invalid selection made! Exiting!")

if searchType == 1:
	mspepOptions = "dzGaihv"
elif searchType == 2:
	mspepOptions = "dyGailv"
elif searchType == 3:
	mspepOptions = "duGaihv"

MFcutoff = input(" > Match Factor minimum cutoff [400 default]? ")
MFcutoff = 400 if MFcutoff == "" else float(MFcutoff)

SMPcount = input(" > Multiprocessing thread count [12 default]? ")
SMPcount = 12 if SMPcount == "" else int(SMPcount)


def smp_MetaboFinding(theBin, progressQu, aLetter, debug=False):
	#print len(theBin)

	file_path = 'inputFileBATCH_' + aLetter + '.msp'

	inputFileMSPEP = open(file_path, 'w')
	outputFilename = "outfile" + aLetter + ".tsv"

	if debug: print("Process %s: analyzing %s spectra..." % (aLetter, len(theBin)))


	for specKey in theBin:
		specData = specAnalysis[specKey]

		monoIsoMZ = specData[0]
		trueCharge = specData[1]
		theSpectrum = specData[2]
		srcFile = specData[3]
		#hitCount = GPfinder(monoIsoMZ, trueCharge, theSpectrum, specKey, srcFile)

		inputFileMSPEP.write("Name: %s\nPrecursorMZ: %s\nIon_mode: %s\nNum Peaks: %s\n" % (specKey, monoIsoMZ, thePolarity, len(theSpectrum)))

		for peak in theSpectrum:
			inputFileMSPEP.write("%s %s\n" % (peak[0], peak[1]))

		inputFileMSPEP.write("\n")



	inputFileMSPEP.close()


	zppm = 20

	if debug: print("Process %s: Launching MSPepSearch..." % (aLetter))

	mspepsearchreadout = subprocess.run(["wine", "MSPepSearch64.exe", mspepOptions,
		"/ZI", "1.6",
		"/ZPPM", str(zppm),
		"/MPPM", "40",
		"/MzLimits", "50", "-1",
		#"/DiffIK", "1", # excludes library match (1K means first segment of InChI)
		"/HITS", "10",
		"/OutBestHitsOnly", "/MatchPolarity",
		"/LIB", "nist_msms", "/HiPri", "/INP", #"/LibInMem",
		file_path,  # this is the input file
		"/OUTTAB", outputFilename,
		"/OutNumComp", "/OutNumMP", "/OutMaxScore",
		"/OutNumPk", "/OutSrchNumPk", "/OutInstrType",
		"/OutCE", "/OutPrecursorMZ",
		"/OutPrecursorType", "/OutChemForm", "/OutIK",
		"/OutNISTrn", "/MinMF", "200"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	theResults = open(outputFilename, "r")
	if debug: print("Process %s: MSpepsearch analysis complete! Analyzing results..." % aLetter)
	processedResults = dict()
	lineCount = 0
	for line in theResults:
		parsedstuff = line.split("\t")
		#print(parsedstuff)

		if lineCount == 3:
			keyIndex = dict([(colName, i) for i, colName in enumerate(parsedstuff)])


		if (lineCount >= 4) & (len(parsedstuff) > 5):

			inputName = parsedstuff[keyIndex["Unknown"]]
			lib_precursorMZ = parsedstuff[keyIndex["Lib Precursor m/z"]]
			compoundName = parsedstuff[keyIndex["Peptide"]]
			theNISTNO = parsedstuff[keyIndex["NIST r.n."]]
			InChIKey = parsedstuff[keyIndex["InChIKey"]]
			theMF = float(parsedstuff[keyIndex["Score"]])
			adductType = parsedstuff[keyIndex["Prec.Type"]]
			InstType = parsedstuff[keyIndex["Instr.Type"]]

			specCount = specAnalysis[inputName][4]
			predictedAdduct = specAnalysis[inputName][1]
			mp_precursorMZ = round(float(specAnalysis[inputName][0]), 4)


			if searchType == 2:
				deltaMass = parsedstuff[keyIndex["DeltaMass"]]
			elif searchType == 3:
				deltaMass = round(mp_precursorMZ - float(lib_precursorMZ), 4)

			if specCount == 1:
				specType = "Single"
			else:
				specType = "Majority (%s constituent spectra)" % specCount


			if theMF > MFcutoff:

				finalOutput = [theMF, compoundName, theNISTNO, InChIKey, lib_precursorMZ, adductType, InstType, specType, predictedAdduct, mp_precursorMZ]
				if searchType != 1:
					finalOutput.append(deltaMass)

				if inputName in processedResults:
					processedResults[inputName].append(finalOutput)
				else:
					processedResults[inputName] = [finalOutput]


		lineCount += 1

	if debug: print("Process %s: Analysis complete! %s out of %s experimental spectra with search results" % (aLetter, len(processedResults), len(theBin)))

	filteredResults = []
	for inputName in processedResults:
		topResult = max(processedResults[inputName])
		theMF, compoundName, theNISTNO, InChIKey, lib_precursorMZ, adductType, InstType, specType, predictedAdduct, mp_precursorMZ = topResult[:10]
		outputLine = [inputName, specType, predictedAdduct, mp_precursorMZ, compoundName, lib_precursorMZ, InChIKey, theMF, theNISTNO, adductType, InstType]

		if searchType != 1:
			deltaMass = topResult[-1]
			outputLine.append(deltaMass)

		filteredResults.append(outputLine)
	progressQu.put([len(theBin), len(processedResults), filteredResults])

	os.remove(file_path)
	os.remove(outputFilename)



def MetaboFinding(SMPcount, allParcels):
	totalProgress = 0
	successProgress = 0

	outputHeader = ["Cluster"] + orderedFileList + ["Query Spectrum Type", "Predicted Adduct Type", "Input precursor m/z", "Compound Name", "Lib precursor m/z", "InChIKey", "placeHolderMF", "NISTNO", "Lib adduct", "Instrument"]

	if searchType == 1:
		outputHeader[-4] = "Direct MF"
		appendText = "Direct"
	elif searchType == 2:
		outputHeader[-4] = "Hybrid MF"
		outputHeader.append("DeltaMass")
		appendText = "Hybrid"
	elif searchType == 3:
		outputHeader[-4] = "In-Source MF"
		outputHeader.append("Precursor m/z difference")
		appendText = "InSource"

	finalSpreadsheet = open(justName + "_NISTsearchResults_" + appendText + ".csv", 'w')
	finalSheetCSV = csv.writer(finalSpreadsheet, delimiter=",")

	letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R"]

	parcelMax = len(allParcels)


	progressQu = Queue()

	smp_processes = [[Process(target=smp_MetaboFinding, args=(allParcels[i], progressQu, letters[i])), letters[i]] for i in range(SMPcount)]
	nextParcelIndex = SMPcount - 1

	[proc[0].start() for proc in smp_processes]

	#sleep(5)
	itsAlive = True
	#for proc in smp_processes:
	#	print proc.is_alive()



	finalSheetCSV.writerow(outputHeader)

	while(itsAlive | (not progressQu.empty())):

		if(progressQu.empty() is not True):
			progressNum = progressQu.get()

			totalProgress += progressNum[0]
			successProgress += progressNum[1]

			sys.stdout.write("\r >>> %s spectra analyzed, %s successfully IDed! <<<" % (totalProgress, successProgress))
			sys.stdout.flush()

			collatedResults = progressNum[2]

			for line in collatedResults:
				clustName = line[0]
				chromStatus = specAnalysis[clustName][5]

				chromReadout = []
				for _chromName in orderedFileList:
					chromReadout.append(chromStatus[_chromName])

				finalSheetCSV.writerow([clustName] + chromReadout + line[1:])


		itsAlive = False
		for i in range(SMPcount):
			proc = smp_processes[i][0]

			if proc.is_alive() is False:

				nextParcelIndex += 1
				if nextParcelIndex < parcelMax:
					proc.join()
					newProc = Process(target=smp_MetaboFinding, args=(allParcels[nextParcelIndex], progressQu, smp_processes[i][1]))

					newProc.start()

					smp_processes[i][0] = newProc

			else:
				itsAlive = True


	#end
	[proc[0].join() for proc in smp_processes]
	progressQu.close()

	finalSpreadsheet.close()

# --------------------------------------- main ------------------------------------------
allParcels = []
_binned = [[] for i in range(SMPcount)]

for i, specKey in enumerate(specAnalysis):

	_binned[i % SMPcount].append(specKey)

	if (i + 1) % 2000 == 0:
		allParcels.extend(_binned)
		_binned = [[] for i in range(SMPcount)]

allParcels.extend(_binned)

print("\nAnalyzing %s potential metabolites..." % clusterCount)

MetaboFinding(SMPcount, allParcels)

print("\nComplete!")



