// Short description:
// 1. image is openend (.czi format!), microscopy images stained for DAPI (blue color channel #2) and 
// 		gH2AX (green color channel #1). 
// 2. DAPI color channel has its background removed and is gauss-filtered. Then Default thresholding is applied.
// 		Small islands of size smaller than what is set are removed.
// 3. An Overlay of all images (segmentation, DAPI & gH2AX) is displayed, user selects suitable cells for analysis.
// 4. All selected cells are processed in a loop. Measured quantities:
//		- Nucleus size
// 		- number of foci
// 		- mean +/- SD size of the foci as full area at half maximum intensity
//		- Index of ROI
// 5.	All ROIs and the segmented DAPI channel are stored in a separate directory.

// Copyright 2019, Johannes Müller, OncoRay Dresden, johannes.mueller@uniklinikum-dresden.de
// 
// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software without 
//    specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
//////////////////////////////////////// DEFAULT SETTINGS//////////////////////////////////////////////////////
/*
saveImgs 		= true;
measure_fSize 	= true;
prominence 	= 2000.0; // determines how sensitive the peak finder is to smaller maxima
CutOff  = 0.5; // value for relative cutoff (Cutoff * (Max - Min)) to determine foci area: 1 = full area, 0 = Zero Area, 0.5 = full area at half maximum

ChFoci = 1; // which channel contains the foci staining?
ChDAPI =2; // which channel contains the DAPI staining?
*/

Validation_mode = false; // set this to "true" if you only want to validate the foci sizes.

////////////////////////////////////////START//////////////////////////////////////////////////////

// Clean up
run("Close All");
roiManager("reset");
run("Clear Results");


// Input GUI
#@ File (label="Input image") filename

#@ String (visibility=MESSAGE, value="Processing parameters", required=false) a
#@ Integer (label="Foci prominence", min=0, max=1000, value=10) prominence
#@ Float (label="Intensity cutoff", min=0, max=1.0, value=CutOff) CutOff

#@ String (visibility=MESSAGE, value="Image parameters", required=false) aaa
#@ Float (label="Pixel size (microns)", value=0.16), pixSize
#@ Integer (label="gH2AX channel", min=1, max=3, value=ChFoci) ChFoci
#@ Integer (label="DAPI channel", min=1, max=3, value=ChDAPI) ChDAPI

#@ String (visibility=MESSAGE, value="Cell inclusion parameters", required=false) aa
#@ Float (label="Minimal nucleus size (microns)", style="slider", min=0, max=100, stepSize=0.1, value = 30) min_size
#@ Float (label="Brightness: lower percentile", style="slider", min=0, max=1, stepSize=0.01, value=0) lower_perc
#@ Float (label="Brightness: upper percentile", style="slider", min=0, max=1, stepSize=0.01, value=1) upper_perc

#@ Boolean (label="Save cell-wise images",  saveImgs=saveImgs) saveImgs


dir = File.getParent(filename);


if (isOpen("Measurements")) {
	close("Measurements");	
}

//Allocate arrays and set measurements
run("Set Measurements...", "area mean min median center display redirect=None decimal=2");
CellMask   = "CellMask";
SegMap     = "SegMap";

// Create Results table
run("Table...", "name=[Measurements] width=800 height=600");
print("[Measurements]", "\\Headings:Label \tArea \tN_Foci \tType \tFSize_Mean \tFSize_Std \tMutualdistance");

// Load Data
run("Bio-Formats (Windowless)", "open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
Image = File.nameWithoutExtension;
rename(Image);

savepath = File.getParent(filename) + "\\" + Image + "_results\\";
File.makeDirectory(savepath);

// get datatype of Image
selectWindow(Image);
bD = bitDepth();

// is it an RGB image?
if (bD == 24) {
	run("Make Composite");
	ChFoci = 2;
	ChDAPI = 3;
}

// were pixelsizes set?
getPixelSize(unit, pixelWidth, pixelHeight);
if (unit == "inch") {
	run("Set Scale...", "distance=1 known=" + pixSize + " unit=microns global");
}

CellSeg(Image, ChDAPI);  // segment DAPI image with Stardist 2D
cleanROIs();  // remove cells that don't pass brightness/area criterion

reply = getBoolean("This is your label map of segmented nuclei. Happy with it?");
if (reply == 0) {
	exit();	
}


NCells = roiManager("count");
Cellname = newArray(roiManager("count"));  // this array stores the names of all rois for later ID

// Process remaining nuclei
selectWindow(Image);

for (i = 0; i < NCells; i++) {
	selectWindow(Image);
	setSlice(ChDAPI);
	
	roiManager("select", i);
	Cellname[i] = call("ij.plugin.frame.RoiManager.getName", i);  // add the name of this cell to the collection
	
	// Measure the size
	run("Measure");
	SizeNucleus = getResult("Area", nResults()-1);

	// A duplicate of the cell is created here. Can be saved, or discarded (see option above)
	// add first channel to duplicate. This is done channelwise to avoid channel/slice/timepoint confusion
	run("Duplicate...", " ");
	rename(Cellname[i]);
	run("Add Slice", "add=slice");
	run("Add Slice", "add=slice");


	// add 2nd channel (foci) to duplicate
	selectWindow(Image);
	setSlice(ChFoci);
	run("Copy");
	run("Blue");
	selectWindow(Cellname[i]);
	setSlice(2);
	run("Paste");
	run("Green");
	
	// add 3rd channel (DAPI mask) to stack
	setSlice(3);
	run("Set...", "value=1 slice");
	run("Grays");

	// Actually measure number of foci
	// The following part of the script is basically the re-evaluation (to be written)
	setSlice(2);
	run("Select None"); // remove selection so that selection type is recognized properly
	run("Find Maxima...", "prominence="+prominence+" strict exclude output=[Point Selection]");
	
	// Store number of foci
	if (selectionType() == -1) {
		nFoci = 0;
		mean = NaN;
		stdDev = NaN;
	} else {
		getSelectionCoordinates(x, y);
		nFoci = x.length;
	}


	// This part measures the size of a foci
	if (nFoci > 0) {
		selectWindow(Cellname[i]);

		// get background
		setSlice(3);
		selectForeground(Cellname[i]);
		run("Create Selection");
		resetThreshold();
		setSlice(2);

//		BG = getPercentile(Cellname[i], 25);

		// Divide image in ROIs that contain one foci each.
		// If there's only one foci, particle segmentation fails and the
		// entire nucleus is considered as a single area of interest.
		if (nFoci > 1) {
			run("Find Maxima...", "prominence=" + prominence + " strict exclude output=[Segmented Particles]");
			setThreshold(128, 255);  // make sure correct area is selected
			run("Create Selection");
			close();
		} else {
			run("Restore Selection");
		}
		
		// Make a copy of DAPI mask and imprint segmented particles (aka foci ROIs) in this mask
		selectWindow(Cellname[i]);
		setSlice(3);
		run("Copy");
		run("Add Slice", "add=slice");
		run("Paste");
		run("Restore Selection");

		/*
		*	divide the binary nucleus area into separated regions with zero-value pixels;
		*	This way, a part of the nucleus area can clearly be assigned to each foci.
		*	Within this separated "puzzle piece", the area of a foci can then be measured
		*/
		if (nFoci > 1) {
			run("Clear Outside", "slice");
		}
		
		/*
		 * Select i-th cell that's currently analyzed (or cell subimage, respectively)
		 * Select Foreground (union of now clearly separated foci regions) and split into separate ROIs.
		 * Remove the initial ROI (union of all pieces) so that only the puzzle pieces remain
		 */

		selectWindow(Cellname[i]);
		selectForeground(Cellname[i]);
		roiManager("Add");
		roiManager("Select", roiManager("count")-1);
		if (nFoci > 1) {
			roiManager("Split");
			roiManager("Select", NCells); // delete selection that wasn't split
			roiManager("Delete");
		}
		
		Foci_sizes = newArray(nFoci);
		Foci_x = newArray(nFoci);
		Foci_y = newArray(nFoci);

		/*
		 * Iterate over all foci regions (a.k.a. puzzle pieces)
		 * ROIs Nr. 0 - NCells are attributes of the previous cell segmentation
		 * ROIs Nr. Ncells + 1 - NROIs are foci regions
		 * Measure Area and other stuff of each foci
		*/
		for (j = NCells; j < roiManager("count"); j++) {
			processPiece(Cellname[i], j, 2, CutOff, Validation_mode);			
			run("Copy");
			close();
			selectWindow(Cellname[i]);
			run("Paste");
			resetThreshold();
			// store area of each foci
			Foci_sizes[j-NCells] = getResult("Area", nResults() - 1);
			Foci_x[j-NCells] = getResult("XM", nResults() - 1);
			Foci_y[j-NCells] = getResult("YM", nResults() - 1);
		}

		// measure mean size/std of all foci
		Array.getStatistics(Foci_sizes, min, max, mean, stdDev);

		// Median distance of foci to its neighbours
		MutDist = MutualDistance(Foci_x, Foci_y);

		// clean up ROI Manager: Keep Cell ROIs (0-NCells), remove puzzle piece ROIs (NCells - N)
		N = roiManager("count");
		for (j = NCells; j < N; j++) {
			roiManager("Select", NCells); // delete selection that wasn't split
			roiManager("Delete");
		}

		// do the point selection again so that images can be stored with the point selection.
		// May be easier for later introspection.
		resetThreshold();
		setSlice(2);
		run("Find Maxima...", "prominence="+prominence+" strict exclude output=[Point Selection]");
	}

	// If images of each cell should be stored (no matter how many foci).
	if (saveImgs){
		if (Validation_mode) {
			// If data comes from validation mode, add this to the filename.
			saveAs(".tif", savepath + Cellname[i] + "_valid");
			rename(Cellname[i]);
			close();
		} else {
			// If data comes from automated analysis, specify this in the file name
			saveAs(".tif", savepath + Cellname[i] + "_auto");
			rename(Cellname[i]);
			close();
		}

	} else {
		close();
	}

	// print cell-averaged results to table
	print("[Measurements]", Cellname[i] +"\t" + 
							SizeNucleus +"\t"+ 
							nFoci+"\t" + 
							"Average" + "\t" + 
							mean +"\t"+ 
							stdDev + "\t" + 
							MutDist);

	// print foci-wise results in results to table
	for (f = 0; f < Foci_sizes.length; f++) {
		print("[Measurements]", Cellname[i] +"\t" +
								" " +"\t" +
								" " +"\t" +
								"Single" + "\t" +
								Foci_sizes[f] +"\t"+
								" ");
	}
}

// When results are stored, add different ending to table.
// Otherwise, original measurements would be overwritten in valid mode.
if (Validation_mode) {
	roiManager("Save", savepath + "RoiSet_" + Image + ".zip");
	selectWindow("Measurements");
	saveAs("Measurements", savepath + "results_" + Image + "_valid.csv");
} else {
	roiManager("Save", savepath + "RoiSet_" + Image + ".zip");
	selectWindow("Measurements");
	saveAs("Measurements", savepath + "results_" + Image + "_auto.csv");
}

run("Close All");
selectWindow("Results");
run("Close");

function MutualDistance(X, Y){
	/*
	 * This function takes Vectors <x> and <y> with the center of mass
	 * coordinates of measured foci. The euclidian distance of each foci
	 * to all other foci is then measured. This gives and estimate
	 * on how densely foci are present in the analyzed cell.
	 */

	N = X.length;
	output = newArray();

	// if there's only one foci, this measurement makes no sense
	if (N == 1) {
		return NaN;
	}

	// otherwise, iterate over coordinates
	for (i = 0; i < N; i++) {
		x = X[i];
		y = Y[i];
		for (j = 0; j < N; j++) {
			_x = X[j];
			_y = Y[j];
			d = Math.sqrt( Math.pow(x - _x, 2) + Math.pow(y - _y, 2));	//euclidian distance between x and _x
			print(d);
			if(d > 0.001){
				output = Array.concat(output, d); 	// append only non-zero result to array. In other words: ignore distance to self
			}
		}
	}
	
	// get median distance and return
	output = Array.sort(output);
	return output[floor(output.length/2)];
}

function processPiece(Image, ROI, Slice, p, valid_mode){
	// this function looks at a part of a cell that contains one (and only one) Foci
	// and measures its area.
	// Input:
	// Image - name of the image to process.
	// ROI - ROI of this piece
	// Slice - If the Image has more than one slice, specify it
	// p - Threshold specification: percentage of max. intensity above which a pixel counts as foci

	selectWindow(Image);
	roiManager("Select", ROI);

	// duplicate piece
	run("Duplicate...", "duplicate");
	rename("PuzzlePiece");
	setSlice(Slice);

	// if automated analysis is desired: first bracket.
	// else: automated analysis is replaced by manual foci delineation.
	run("Measure");
	Max = getResult("Max", nResults() -1);
	BG = getPercentile("PuzzlePiece", 5);

	if (valid_mode == false) {
		setThreshold(BG + p * (Max - BG), Max);
		run("Convert to Mask", "background=Dark black");
		selectForeground("PuzzlePiece");
		run("Measure");
		run("16-bit");
		run("Multiply...", "value=1000 slice");
		run("Select None");
	} else {
		run("Clear Outside", "slice");
		run("Select None");
		setMinAndMax(BG, Max);
		run("In [+]"); // enlarge window a bit.
		run("In [+]");
		run("In [+]");
		run("In [+]");
		setTool("freehand"); // set free ROI tool
		waitForUser("Input needed", "draw outline of Foci in image. click 'ok' when done.");
		run("Clear Outside", "slice");
		run("Set...", "value=100000");
		selectForeground("PuzzlePiece");
		run("Measure");
		run("Select None");
	}
}

function PercOfArray(my_array, perc){
	// returns a percentile of a numeric vector <my_array>
	my_array = Array.sort(my_array);
	index = floor(perc * (my_array.length -1));

	return my_array[index];
}

function getPercentile(Image, perc){

	selectWindow(Image);

	// get histogram
	nBins = 256;
	getHistogram(values, counts, nBins);

	// get total number of counts
	total = 0;
	for (i = 0; i < counts.length; i++) {
		total = total + counts[i];
	}

	// go through histogram until threshold is reached
	c = counts[0];
	i = 0;
	do {
		i = i+1;
		c = c + counts[i];
	} while (c < perc/100 * total);

	print(perc +"% threshold of " + Image + ": " + 0.5*(values[i+1] + values[i]));

	return 0.5*(values[i+1] + values[i]);
}

function enhanceVisibility(Image, DAPIchannel, fociChannel){
	selectWindow(Image);
	run("Make Composite");
	setSlice(DAPIchannel);
	run("Grays");
	setSlice(fociChannel);
	run("Green");
	run("RGB Color");
	run("Enhance Contrast", "saturated=0.03");
}

function CellSeg(Input, channel) {
	// Segment Nuclei
	run("Duplicate...", "duplicate channels="+channel);
	rename("CellMask");
	run("Subtract Background...", "rolling=120");
	run("Gaussian Blur...", "sigma=3");
	setAutoThreshold("Default dark");
	run("Convert to Mask");
	run("Watershed");

	return "CellMask";
}

function selectForeground (Input){
	selectWindow(Input);
	setThreshold(1, 100000);
	run("Create Selection");
}

function CellSeg(Image, ChDAPI){
	// run StarDist

	print(ChDAPI);

	selectWindow(Image);
	setSlice(ChDAPI);
	run("Duplicate...", "title=DAPI");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], " +
			"args=['input':'DAPI', "+
			"'modelChoice':'Versatile (fluorescent nuclei)', "+
			"'normalizeInput':'true', "+
			"'percentileBottom':'5.0', "+
			"'percentileTop':'97.30000000000001', "+
			"'probThresh':'0.3500000000000002', "+
			"'nmsThresh':'0.4', "+
			"'outputType':'Both', "+
			"'nTiles':'1', "+
			"'excludeBoundary':'15', "+
			"'roiPosition':'Automatic', "+
			"'verbose':'false', "+
			"'showCsbdeepProgress':'false', "+
			"'showProbAndDist':'false'], "+
			"process=[false]");
}

function cleanROIs(){
	NCells = roiManager("count");

	
	// convert Measurements to arrays
	roiManager("Measure");
	NucleiSizes = newArray(NCells);
	NucleiBrightness = newArray(NCells);
	for (i = 0; i < NCells; i++) {
		NucleiSizes[i] = getResult("Area", i);
		NucleiBrightness[i] = getResult("Mean", i);
	}

	// filter ROI list according to cell inclusion parameters
	lower_brightness = PercOfArray(NucleiBrightness, lower_perc);
	upper_brightness = PercOfArray(NucleiBrightness, upper_perc);
	run("Clear Results");

	selectWindow("Label Image");

	to_be_removed = newArray();
	// remove bad ROIs from list
	for (i = 0; i < NCells; i++) {
		roiManager("select", i);
		roiManager("measure");
		area = getResult("Area", nResults - 1);
		mean_brightness = getResult("Mean", nResults - 1);

		// Is the brightness of this nucleus outside specified range?
		if ((mean_brightness < lower_brightness) || (mean_brightness > upper_brightness)) {
			to_be_removed = Array.concat(to_be_removed, i);
			print("brightness wrong!");
			run("Clear");
		}

		// Is the area of this nucleus too small?
		if (area < min_size) {
			to_be_removed = Array.concat(to_be_removed, i);
			print("too small!");
			run("Clear");
		}
	}

	// remove from ROI Manager
	roiManager("select", to_be_removed);
	roiManager("delete");
	run("Clear Results");
	
}
