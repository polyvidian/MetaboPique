#!/usr/bin/env python2
import sys, os, re, math
from collections import OrderedDict
import tkinter as Tk
#from Tkinter import *
from tkinter import ttk
from tkinter import font
import matplotlib
matplotlib.use('TkAgg')
#from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

try: from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as PlotNav
except: from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as PlotNav
# implement the default mpl key bindings
#from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

precursorPatt = re.compile(r'Precursor\sm/z:\s(\d+.\d+)')
spectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.\d+)')
#majSpectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.\d+)\s\((\d+)/(\d+)\)\s(\d+\.\d+)ppm')
majSpectrumPatt = re.compile(r'(\d{2,4}\.\d{1,4})\s(\d+\.\d+)\s(.+)')

#chromPatt = re.compile(r'(\d{1,7}\.\d{0,4}),(\d+\.\d+),(\d+\.\d+)')
#chromPatt = re.compile(r'(\d{1,7}\.\d{0,4}),\d+\.{0,1}\d*(?:,\d+\.{0,1}\d*)+')
chromPatt = re.compile(r'(\d{1,7}\.\d{0,4}),\d+\.{0,1}\d*(?:,\d+\.{0,1}\d*)*')
AbPatt = re.compile(r',\d+\.{0,1}\d*')
isoMZPatt = re.compile(r'\d+\.\d+')

filenamePatt = re.compile(r'Chromatogram:\s(.+\..+)')
peakCountPatt = re.compile(r'Peak\scount:\s(\d+)\s')
chromCountPatt = re.compile(r'Chromatogram\scount:\s(\d+)')

allMSdataPath = sys.argv[1]

allMSdata = next(iter(os.walk(allMSdataPath)))
rootDirectory = allMSdata[0]

storedMSdata = dict()
clusterCount = 0; spectrumCount = 0
storeKeys = []

onlySpectra = dict()
allSpectra = dict()
onlyChrom = dict()
chromLoc = dict()  # location of the MS2 spectra on the MS1 chromatographic data

colors = ['r', 'g', 'c', 'deeppink', 'purple', 'forestgreen', 'teal', 'cyan', 'brown', 'orangered', 'magenta', 'cornflowerblue', 'aqua', 'azure', 'orchid', 'navy', 'gold', 'coral', 'lime']

print("Reading files...")
for clusterDir in allMSdata[1]:
	clusterPath = os.path.join(allMSdataPath, clusterDir)
	#clustData = os.walk(os.path.join(allMSdataPath, clusterDir)).next()
	clusterCount += 1


	chromLoc[clusterDir] = [[], [], [], [], []]

	#clustFiles = clustData[2]
	storedMSdata[clusterDir] = [[], []]
	storeKeys.append(clusterDir)

	#for filename in clustFiles:
	#	if ((".txt" in filename) & ("_" in filename)):
	#		spectrumCount += 1


	chromReader = open(os.path.join(clusterPath, "chromatographicData.csv"), "r")
	#storedMSdata[clusterDir]["chromData"] = []

	firstLine = chromReader.readline().rstrip()
	_peakNo = int(re.search(peakCountPatt, firstLine).group(1))
	_chromNo = int(re.search(chromCountPatt, firstLine).group(1))

	onlyChrom[clusterDir] = [set(), [], firstLine]

	chromReader.close()

	storedMSdata[clusterDir][0] = [_peakNo, _chromNo]  # number of peaks, and then number of chromatograms

	if clusterCount % 10 == 0:
		#sys.stdout.write("\r   %s clusters, %s spectra" % (clusterCount, spectrumCount))
		sys.stdout.write("\r   %s clusters" % (clusterCount))
		sys.stdout.flush()

#print "\nComplete! %s clusters containing %s spectra found" % (clusterCount, spectrumCount)
print("\nComplete! %s clusters" % (clusterCount))


#f = plt.figure(figsize=(10,6), dpi=100, facecolor='white')
f = Figure(figsize=(10, 6), dpi=100, facecolor='white')

f2 = Figure(figsize=(6, 6), dpi=100, facecolor='white')


threeDthing = []; lineCollection = dict(); redlineCollection = []



def onpick(event):
	print("Hey this works!")
	global canvas, allLines
	# on the pick event, find the orig line corresponding to the
	# legend proxy line, and toggle the visibility
	legline = event.artist
	_id = legline.get_label()

	for origline in allLines[_id]:
		vis = not origline.get_visible()
		origline.set_visible(vis)
		# Change the alpha on the line in the legend so we can see what lines
		# have been toggled
		if vis:
			legline.set_alpha(1.0)
		else:
			legline.set_alpha(0.2)

	canvas.draw()


def on_key_press(event):
	global canvas
	#print "keypress!!!"
	key_press_handler(event, canvas, theTool)


def treeview_sort_column(tv, col, reverse):
	l = [(tv.set(k, col), k) for k in tv.get_children('')]
	#l.sort(reverse=reverse)
	l.sort(key=lambda t: int(t[0]), reverse=reverse)

	# rearrange items in sorted positions
	for index, (val, k) in enumerate(l):
		tv.move(k, '', index)

	# reverse sort next time
	tv.heading(col, command=lambda: \
			treeview_sort_column(tv, col, not reverse))


ax = ''
allSpecPreviousExistence = ""
lined = dict()


def itemClicked(event):
	specKey = str(left_tree.focus())
	global f, f2, canvas, specCanvas, onlySpectra, threeDthing, spectrumWindow, lineCollection, redlineCollection, ax, allSpecPreviousExistence, allSpectra, allLines
	print(specKey)

	chromKey = specKey.split("|")[0]

	if allSpecPreviousExistence is not chromKey:
		if len(lineCollection) > 0:
			for specNum in list(lineCollection):
				#for aLine in lineCollection[specNum]:
				#	aLine.remove()
				del lineCollection[specNum]


		del redlineCollection[:]

		f2.clf()

		# all spectra 3D plot------------------------------------------------
		if spectrumWindow.winfo_exists() == 0:
			spectrumWindow = Tk.Toplevel()  # Makes the window
			spectrumWindow.wm_title("All Spectra")
			spectrumWindow.geometry("700x650")

			topFrame = Tk.Frame(spectrumWindow, width=700, height=700)
			topFrame.pack(side="top", fill="both", expand=True)
			specCanvas = FigureCanvasTkAgg(f2, master=topFrame)

			specCanvas.get_tk_widget().pack(side="top", fill="both", expand=True)
			theTool2 = PlotNav(specCanvas, topFrame)
			theTool2.update()

		# on demand loading of all spectral data---------------------------------------
		if chromKey not in allSpectra:
			allSpectra[chromKey] = []

			clustData = next(iter(os.walk(os.path.join(allMSdataPath, chromKey))))

			#print os.path.join(allMSdataPath, chromKey)
			clusterPath = os.path.join(allMSdataPath, chromKey)

			majSpecBool = False

			internalSpecCount = 1
			for filename in clustData[2]:
				#print filename
				if ((".txt" in filename) & ("_" in filename)):
					spectrumReader = open(os.path.join(clusterPath, filename), "r")
					_textOut = ""

					_spectrum = []
					_count = 0
					for line in spectrumReader:
						#_textOut += line

						if _count == 0:
							_precursor = float(re.search(r'\d+\.\d+', line).group())
							_textOut += line
						elif _count == 1:
							_intensity = float(re.search(r'\d+\.\d+', line).group())
							_textOut += line
						elif _count == 2:
							_charge = int(re.search(r'\d', line).group())
							_textOut += line
						elif _count == 3:
							_rt = float(re.search(r'\d+\.\d+', line).group())
							_textOut += line
						elif _count == 4:
							_scanNo = float(re.search(r'\d+', line).group())
							_textOut += line
						elif _count == 6:
							_id = re.search(filenamePatt, line).group(1)
							_textOut += line
						else:
							spectrumData = re.match(spectrumPatt, line)
							if spectrumData != None:
								_spectrum.append([float(spectrumData.group(1)), float(spectrumData.group(2))])
								#_textOut += theMZ + " " + str(round(theAb, 2))
							else:
								_textOut += line

						_count += 1

					spectrumReader.close()

					storedMSdata[chromKey][1].append([_precursor, _charge, _rt, internalSpecCount])

					chromLoc[chromKey][0].append(_rt)
					chromLoc[chromKey][1].append(_precursor)
					chromLoc[chromKey][2].append(_intensity)
					chromLoc[chromKey][3].append("S%s" % (internalSpecCount))
					chromLoc[chromKey][4].append(_id)

					numRows = int(math.ceil(len(_spectrum)/6.0))
					_specOut = ["" for i in range(numRows)]


					for i, peak in enumerate(_spectrum):
						strMZ = str(peak[0])
						strAb = str(round(peak[1], 2))

						_specOut[i % numRows] += strMZ+" "+strAb+"\t\t\t"

					_specOut.insert(0, "\n"+"m/z\tAbundance\t\t"*6+"\n"+"-------------------------------\t\t\t"*6)

					for aline in _specOut:
						_textOut += aline+"\n"

					onlySpectra[chromKey+"|"+str(internalSpecCount)] = [_precursor, _charge, _rt, _spectrum, _textOut, internalSpecCount]
					allSpectra[chromKey].append([_spectrum, internalSpecCount])

					internalSpecCount += 1

				if "majority" in filename:
					majSpecBool = True

					majSpectrumReader = open(os.path.join(clusterPath, filename), "r")
					_textOut = ""

					_spectrum = []
					_count = 0
					for line in majSpectrumReader:
						#_textOut += line

						if _count == 0:
							_monoIso = float(re.search(r'\d+\.\d+', line).group())
							_textOut += line
						elif _count == 1:
							if "UNKNOWN" not in line:
								_charge = int(re.search(r'\d', line).group())
							else:
								_charge = 0

							_textOut += line

						elif _count == 2:
							_constituentCount = int(re.search(r'\d+', line).group())
							_textOut += line
						else:
							majspectrumData = re.match(majSpectrumPatt, line)
							if majspectrumData is not None:
								_spectrum.append([float(majspectrumData.group(1)), float(majspectrumData.group(2)), majspectrumData.group(3)])

								#_textOut += theMZ + " " + str(round(theAb, 2))
							else:
								_textOut += line

						_count += 1

					majSpectrumReader.close()

					numRows = int(math.ceil(len(_spectrum)/5.0))
					_specOut = ["" for i in range(numRows)]

					_maxMZmaj = max(_spectrum)[0]
					_minMZmaj = min(_spectrum)[0]
					for i, peak in enumerate(_spectrum):
						strMZ = str(peak[0])
						strAb = str(round(peak[1], 2))
						if len(strAb) > 5:
							strAb = str(round(peak[1], 1))
						anno = peak[2]

						_specOut[i % numRows] += strMZ+" "+strAb+"\t"+anno+"\t\t\t"

					_specOut.insert(0, "\n"+"m/z\tAb   Freq  Error\t\t\t"*5+"\n"+"---------------------------------------\t\t\t\t"*5)

					for aline in _specOut:
						_textOut += aline+"\n"

					onlySpectra[chromKey+"|Majority"] = [_monoIso, _charge, "N/A", _spectrum, _textOut, "Majority", _maxMZmaj, _minMZmaj]


			# add to the tree!
			if majSpecBool:
				left_tree.insert(chromKey, "end", chromKey+"|Majority", text="MajoritySpectrum", values=("", "", _monoIso, "N/A"))

			for spectrumData in storedMSdata[chromKey][1]:
				spectrumCount = spectrumData[3]
				specTreeName = chromKey + "|" + str(spectrumCount)
				_prec = spectrumData[0]
				_rt = spectrumData[2]

				left_tree.insert(chromKey, "end", specTreeName, text="Spectrum"+str(spectrumCount), values=("", "", _prec, _rt))

			# ------------------ on demand chromatographic data loading
			chromReader = open(os.path.join(clusterPath, "chromatographicData.csv"), "r")

			firstLine = chromReader.readline().rstrip()
			chromReader.readline()

			_chromStats = []
			for line in chromReader:
				if "enter: " in line:
					_initialInfo = line.rstrip()
					_precursorChrom = float(re.search(precursorPatt, _initialInfo).group(1))
					_identifier = re.search(filenamePatt, _initialInfo).group(1)
				elif "Retention" in line:
					_isoMZs = [float(isoMZ) for isoMZ in re.findall(isoMZPatt, line)]

				else:
					chromData = re.match(chromPatt, line)
					if chromData != None:
						theAbs = re.findall(AbPatt, line)
						_chromStats.append([float(chromData.group(1))] + [float(abVal[1:]) for abVal in theAbs])

				if line.rstrip()=="":
					#storedMSdata[chromKey]["chromData"].append([_initialInfo, _precursorChrom, _chromStats])
					_sortedChromStats = sorted(_chromStats)
					onlyChrom[chromKey][0].add(_identifier)
					onlyChrom[chromKey][1].append([_initialInfo, _precursorChrom, _sortedChromStats, _identifier, _isoMZs])
					_chromStats = []
					_initialInfo = ""
					_precursorChrom = "WRONG"

			#if len(storedMSdata[chromKey]["chromData"]) > 1:
			#	print storedMSdata[chromKey]["chromData"]
			chromReader.close()

		#----------------------------------------------------------- 3D drawing of all spectra in the cluster

		allTheSpectra = allSpectra[chromKey]

		#if len(allTheSpectra) > 1:
		ax = f2.add_subplot(111, projection='3d')
		f2.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02, wspace=0.02)

		ax.set_xlabel('m/z')
		ax.set_ylabel('Spectrum Index')
		ax.set_zlabel('Intensity')


		allMZ = []; allIndex = []

		for threedspec in allTheSpectra:
			theSpectrum = threedspec[0]

			_theIndex = threedspec[1]
			lineCollection[_theIndex] = []

			allIndex.append(_theIndex)

			for peak in theSpectrum:
				_theMZ = peak[0]
				_theAbs = peak[1]

				if _theAbs > 50:

					line = art3d.Line3D(*zip((_theMZ, _theIndex, 0), (_theMZ, _theIndex, _theAbs)), marker=' ', markevery=(1, 1))
					ax.add_line(line)

					lineCollection[_theIndex].append(line)

					allMZ.append(_theMZ)

		ax.set_xlim(min(allMZ), max(allMZ))
		ax.set_ylim(0, len(allTheSpectra)+1)
		ax.set_yticks(allIndex)
		ax.set_zlim(0, 1100)

		toolbar2 = f2.canvas.toolbar
		toolbar2.update()
		specCanvas.draw()

		allSpecPreviousExistence = chromKey


	# ------------------------------- 2D spectrum viewer
	if specKey in onlySpectra:

		#plt.clf()
		if len(threeDthing) > 0:
			a3dthing = threeDthing.pop()
			f.delaxes(a3dthing)
			#a3dthing.remove()
			del a3dthing



		#a = f.add_subplot(111)
		a = f.add_subplot(111)

		f.subplots_adjust(left=0.06, right=0.94, top=0.9, bottom=0.1)
		a.stem([0], [0], markerfmt=' ')
		a.spines['top'].set_visible(False)
		a.spines['right'].set_visible(False)
		#print specKey
		precursorMZ = onlySpectra[specKey][0]
		filteredSpectrum = onlySpectra[specKey][3]
		textOut = onlySpectra[specKey][4]
		specNum = onlySpectra[specKey][5]

		localMajorityBool = False
		if chromKey+"|Majority" in onlySpectra:
			localMajorityBool = True
			majMZmax = onlySpectra[chromKey+"|Majority"][6]
			majMZmin = onlySpectra[chromKey+"|Majority"][7]


		if len(redlineCollection) > 0:
			previousSpecNum = redlineCollection[0][1][1]
			for aLine in lineCollection[previousSpecNum]:
				ax.add_line(aLine)

			[rLine[0].remove() for rLine in redlineCollection]

			redlineCollection = []

		if specNum != "Majority":
			for aLine in lineCollection[specNum]:
				aLine.remove()

		sortedSpectrum = sorted(filteredSpectrum)
		xData = []; yData = []

		for anAttr in sortedSpectrum:
			xData.append(anAttr[0])
			yData.append(anAttr[1]/10)

		if localMajorityBool is False:
			theRange = (max(xData)+10)-(min(xData)-10)
			a.set_xlim(min(xData)-10, max(xData)+10)
		else:
			theRange = (majMZmax+10)-(majMZmin-10)
			a.set_xlim(majMZmin-10, majMZmax+10)
			print("worked!")

		a.arrow(float(precursorMZ), 3, 0.0, -1, fc="k", ec="k", head_width=theRange*0.01, head_length=2)
		a.stem(xData, yData, 'silver', markerfmt=' ')

		a.set_ylim(0, 110)
		a.axvline(0, color='black', linewidth=2)

		a.grid(color='lightgray', alpha=0.7)

		a.set_xlabel("m/z", fontsize=12)
		a.set_ylabel("Relative Intensity", fontsize=12)

		whatToPlot = []
		theRecords = []  # you use this to check if there is any peak that is in the range of your current peak that would get in the way of the text (i.e. it's a higher intensity)
		for mz_intensity in sortedSpectrum:
			#this code will only allow text to be plotted either if the preceding plotted text is far enough away ( > theRange*0.05), or if it has a higher intensity than the preceding plotted text (which it subsequently removes)
			originalIntensity = mz_intensity[1]

			_xcoor = mz_intensity[0]
			_ycoor = originalIntensity/10.0

			if ((originalIntensity > 50) & (specNum != "Majority")):
				line = art3d.Line3D(*zip((_xcoor, specNum, 0), (_xcoor, specNum, originalIntensity)), marker=' ', c='r', markevery=(1, 1))
				ax.add_line(line)
				redlineCollection.append([line, (_xcoor, specNum, originalIntensity)])

			if len(whatToPlot) > 0:
				last_xcoor = whatToPlot[-1][0]
				last_ycoor = whatToPlot[-1][1]

				#print "%s --> %s" % (_xcoor, _xcoor - last_xcoor)

				if (_xcoor - last_xcoor < theRange*0.05) & (_ycoor > last_ycoor) & (_ycoor > 1):  # if it is NOT far enough away BUT it is taller than the last text then go ahead
					whatToPlot.pop()
					whatToPlot.append([_xcoor, _ycoor])
					#print "popped! %s" % (last_xcoor)
				elif (_xcoor - last_xcoor > theRange*0.05) & (_ycoor > 1):  # if it is far enough away from the last text then go ahead
					goAhead = True
					for checkPrevious in theRecords:
						if (_xcoor - checkPrevious[0] < theRange*0.05) & (_ycoor < checkPrevious[1]):  # check to see if anything in range that is also greater than the current intensity
							goAhead = False
					if goAhead == True: whatToPlot.append([_xcoor, _ycoor])
			else:
				if (_ycoor > 1): whatToPlot.append([_xcoor, _ycoor])

			theRecords.append([_xcoor, _ycoor])

		#a.arrow(float(precursorMZ), 3, 0.0, -1, fc="k", ec="k", head_width=theRange*0.05, head_length=5)


		for dataToPlot in whatToPlot:
			a.annotate(str(dataToPlot[0]), xy=(dataToPlot[0], dataToPlot[1]), ha='center', size='x-small')


		toolbar1 = f.canvas.toolbar
		toolbar1.update()
		canvas.draw()

		toolbar2 = f2.canvas.toolbar
		toolbar2.update()
		specCanvas.draw()

		threeDthing.append(a)

		specData.delete(1.0, Tk.END)
		specData.insert(0.0, textOut)

	# ------------------------------- MS1 chromatogram plot
	elif specKey in onlyChrom:

		if len(threeDthing) > 0:
			a3dthing = threeDthing.pop()
			f.delaxes(a3dthing)
			#a3dthing.remove()
			del a3dthing

		f.clf()

		ax2 = f.add_subplot(111, projection='3d')

		ax2.set_xlabel('retention time (seconds)')
		ax2.set_ylabel('m/z')
		ax2.set_zlabel('Raw Abundance')

		IDset = onlyChrom[specKey][0]
		fileNames = dict([(_id, (i*0.01, i)) for i, _id in enumerate(sorted(IDset))])
		allLines = dict([(_id, []) for _id in IDset])

		#print fileNames
		theData = onlyChrom[specKey][1]
		specLocations = chromLoc[specKey]

		allMZs = []
		textOut = onlyChrom[specKey][2]+"\n\n"

		#print theData
		#allTrueTraceMZs = set([math.floor(aSet[1]*10) for aSet in theData])

		allTrueTraceMZs = dict([(_id, set()) for _id in IDset])
		for aSet in theData:
			_identifier = aSet[3]
			quickSearchVal = round(aSet[1]*10)
			allTrueTraceMZs[_identifier].add(quickSearchVal)
		#print allTrueTraceMZs

		for aSet in theData:
			setLines = []
			_readOutData = aSet[0]
			textOut += _readOutData + "\n"
			_identifier = aSet[3]
			_isoMZs = aSet[4]

			theOffset = fileNames[_identifier][0]
			#theOffset = 0
			colorIndex = fileNames[_identifier][1]

			theLength = len(aSet[2])
			_prec = [aSet[1] for i in range(theLength)]
			_rt = [aSet[2][i][0] for i in range(theLength)]
			_real = [aSet[2][i][1] for i in range(theLength)]

			setLines = []
			#print(aSet[2])
			if len(aSet[2][0]) > 2:
				_fitted = [aSet[2][i][2] for i in range(theLength)]
				theLineFitted, = ax2.plot(_rt, _prec, _fitted, c='mediumblue', linestyle='--', label=_identifier)
				setLines.append(theLineFitted)

			# remember this HAS to come AFTER the fitted gaussian line because the legend color will be all messed up
			theLineReal, = ax2.plot(_rt, _prec, _real, c=colors[colorIndex % len(colors)], label=_identifier)
			setLines.append(theLineReal)

			allMZs.extend(_prec)

			for i, isoMZ in enumerate(_isoMZs):
				if round(isoMZ*10) not in allTrueTraceMZs[_identifier]:
					_isoTrace = [aSet[2][j][i+3] for j in range(theLength)]
					_isoPrec = [isoMZ for i in range(theLength)]
					ghostLine, = ax2.plot(_rt, _isoPrec, _isoTrace, c=colors[colorIndex % len(colors)], label=_identifier, linestyle=':')
					allMZs.extend(_isoPrec)

					setLines.append(ghostLine)


			allLines[_identifier].extend(setLines)


		handles, labels = ax2.get_legend_handles_labels()
		by_label = OrderedDict(zip(labels, handles))
		#theLeg = ax2.legend(list(by_label.values()), list(by_label), loc=2, prop={'size': 7}, bbox_to_anchor=(1, 1))  # remember this moves the legend to not overlap the graph
		theLeg = ax2.legend(by_label.values(), by_label.keys(), loc=2, prop={'size': 9})

		#modYcoors = []
		for x, y, z, theLabel, _id in zip(specLocations[0], specLocations[1], specLocations[2], specLocations[3], specLocations[4]):
			theOffset = fileNames[_id][0]
			colorIndex = fileNames[_id][1]
			theText = ax2.text(x, y, z, theLabel, size=9, color='k', label=_id)

			#modYcoors.append(y+theOffset)

			theScat = ax2.scatter(x, y, z, c=colors[colorIndex % len(colors)], marker='.', s=40, label=_id)
			allLines[_id].extend([theText, theScat])

		[legline.set_picker(5) for legline in theLeg.get_lines()]  # 5 pts tolerance


		minMZ = min(allMZs + specLocations[1]) - 0.05; maxMZ = max(allMZs + specLocations[1]) + 0.05
		ax2.set_ylim(minMZ, maxMZ)

		threeDthing.append(ax2)

		specData.delete(1.0, Tk.END)
		specData.insert(0.0, textOut)
		#ax2.legend(loc=2, prop={'size': 9})

		#ax.set_axis_off
		#ax.scatter(pc1[ctrlN:], pc2[ctrlN:], pc3[ctrlN:], c='b', marker='^',s=60)

		canvas.mpl_connect('pick_event', onpick)
		canvas.mpl_connect('button_press_event', lambda event:canvas._tkcanvas.focus_set())
		canvas.mpl_connect('key_press_event', on_key_press)
		canvas.draw()



root = Tk.Tk()  # Makes the window
root.wm_title("IsoPique Viewer v3.1.0")  # Makes the title that will appear in the top left
root.geometry("1900x1000") # this is probably the secret to preventing the jumpy stuff for the navigation bar!!!!
#root.wm_geometry("")
root.config(background = "#FFFFFF") #sets background color to white

default_font = font.nametofont("TkDefaultFont")
default_font.configure(size=8)
customFont = font.Font(family="Helvetica", size=10)

spectrumWindow = Tk.Toplevel()  # Makes the window
spectrumWindow.wm_title("All Spectra")
spectrumWindow.geometry("700x650")


#put widgets here
leftFrame = Tk.Frame(root, width=160, height = 1200) # this is ABSOLUTELY the secret to preventing a jumpy window!!! expand MUST BE FALSE and the fill MUST BE "y"!!!
leftFrame.pack(side="left", fill="y", expand=False)
#leftFrame.grid(row=0, column=0, padx=5, pady=2)

left_tree = ttk.Treeview(leftFrame, height="50")

ysb = ttk.Scrollbar(orient=Tk.VERTICAL, command=left_tree.yview)
left_tree['yscroll'] = ysb.set
ysb.pack(side="right", in_=leftFrame,fill="y")
#xsb = ttk.Scrollbar(orient=HORIZONTAL, command=left_tree.xview)
#left_tree['yscroll'] = ysb.set
#left_tree['xscroll'] = xsb.set
#left_tree.configure(yscroll=ysb.set, xscroll=xsb.set)


#left_tree.grid(row=0, column=0)

#xsb.grid(row=1, column=0, sticky='ew', in_=leftFrame)

#left_tree.pack(expand=YES, fill=BOTH)
left_tree.heading("#0", text="Group")#, command=lambda: \
			#treeview_sort_column(left_tree, "#0", False))
left_tree.column('#0', stretch="yes", minwidth=100, width=150)
left_tree.pack(side="left", fill="both", expand=True)

#left_tree.configure(yscroll=ysb.set, xscroll=xsb.set)

#left_tree.tag_configure('ttk', background='yellow')
#left_tree.tag_bind('ttk', '<1>', itemClicked); # the item clicked can be found via tree.focus()

left_tree.bind('<<TreeviewSelect>>', itemClicked)

left_tree["columns"]=("pNum","cNum","mz","rt")
left_tree.column("pNum", width=35)
left_tree.column("cNum", width=35)
left_tree.column("mz", width=110)
left_tree.column("rt", width=85)
left_tree.heading("pNum", text="Peaks", command=lambda: treeview_sort_column(left_tree, "pNum", True))
left_tree.heading("cNum", text="Chrom", command=lambda: treeview_sort_column(left_tree, "cNum", True))
left_tree.heading("mz", text="Prec m/z")
left_tree.heading("rt", text="Ret Time")

#left_tree.insert("" , 0,    text="Line 1", values=("1A","1b"))

superOverlayer = dict()

storeKeys = sorted(storeKeys)

for cluster in storeKeys:
	peakCount = storedMSdata[cluster][0][0]
	chromCount = storedMSdata[cluster][0][1]
	left_tree.insert("", "end", cluster, text = cluster, values=(peakCount,chromCount, "", ""))





rightFrame = Tk.Frame(root, width=1110, height = 900)
rightFrame.pack(side="top", fill="both", expand=True)
#alsoRightFrame = Frame(root, width=900, height=900)
#alsoRightFrame.pack(side="bottom",fill="x", expand=False)

#rightFrame.grid(row=0, column=1, padx=10, pady=2)

specData = Tk.Text(rightFrame, width = 160, height = 20, takefocus=0, font=customFont)
#specData.grid(row=4, column=0, padx=10, pady=2)
specData.pack(side="bottom", in_=rightFrame,fill="both", expand=True)

canvas = FigureCanvasTkAgg(f, master=rightFrame)

canvas.get_tk_widget().pack(side="top", fill="both")
theTool = PlotNav(canvas, rightFrame)
theTool.update()

canvas.draw()

topFrame = Tk.Frame(spectrumWindow, width=700, height=700)
topFrame.pack(side="top", fill="both", expand=True)
specCanvas = FigureCanvasTkAgg(f2, master=topFrame)

specCanvas.get_tk_widget().pack(side="top", fill="both", expand=True)
theTool2 = PlotNav(specCanvas, topFrame)
theTool2.update()

specCanvas.draw()

#canvas.get_tk_widget().grid(column=0,row=1)
#canvas.get_tk_widget().pack_forget()



root.mainloop() #start monitoring and updating the GUI. Nothing below here runs.

print("Complete!")
