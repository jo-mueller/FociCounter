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

saveImgs 		= true;
measure_fSize 	= true;
prominence 	= 2000.0; // determines how sensitive the peak finder is to smaller maxima
CutOff  = 0.5; // value for relative cutoff (Cutoff * (Max - Min)) to determine foci area: 1 = full area, 0 = Zero Area, 0.5 = full area at half maximum

ChFoci = 1; // which channel contains the foci staining?
ChDAPI =2; // which channel contains the DAPI staining?

Validation_mode = false; // set this to "true" if you only want to validate the foci sizes.

////////////////////////////////////////START//////////////////////////////////////////////////////

// Input GUI
#@ File (label="Input image") filename
dir = File.getParent(filename);

#@ Boolean (label="Save cell-wise images",  saveImgs=saveImgs) saveImgs
#@ Boolean (label="measure foci size",  measure_fSize=measure_fSize) measure_fSize
#@ Integer (label="Prominence", min=0, max=10000, value=prominence) prominence
#@ Float (label="Cutoff", min=0, max=1.0, value=CutOff) CutOff

#@ Integer (label="gH2AX channel", min=1, max=2, value=ChFoci) ChFoci
#@ Integer (label="DAPI channel", min=1, max=2, value=ChDAPI) ChDAPI


// Clean up
run("Close All");
roiManager("reset");
run("Clear Results");

//Allocate arrays and set measurements
run("Set Measurements...", "area min median display redirect=None decimal=2");
CellMask   = "CellMask";
SegMap     = "SegMap";

// Create Results table
run("Table...", "name=[Measurements] width=800 height=600");
print("[Measurements]", "\\Headings:Label \tArea \tN_Foci	\tType \tFSize_Mean \tFSize_Std");

// Load Dataun("Bio-Formats (Windowless)", "open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
Image = File.nameWithoutExtension
rename(Image);

savepath = File.getParent(filename) + "\\" + Image + "_results\\";
File.makeDirectory(savepath);

// Segment DAPI Nuclei
setSlice(2);
CellMask = CellSeg(Image, ChDAPI);

// Get selection of segmented cells
selectForeground(CellMask);
roiManager("Add");
run("RGB Color");

// Display as overlay to pick suitable cells/exclude dividing cells
// Change LUT for visibility
enhanceVisibility(Image, ChDAPI, ChFoci);

// Time for the user to do the work
selectWindow(CellMask);
run("Add Image...", "image=["+Image+" (RGB)] x=0 y=0 opacity=70");
waitForUser("Use FloodFill Tool to mark well-segmented cells for analysis.");

// Post-process result
run("Remove Overlay");
run("Select None");
run("8-bit");
close(Image+" (RGB)");

// Get selected cells
setThreshold(1, 244);
run("Convert to Mask");

// Manage ROIs
roiManager("reset");
run("Analyze Particles...", "size=0-Infinity pixel show=Nothing add");

nFoci = newArray(roiManager("count"));
Cellname     = newArray(roiManager("count")); // Array for all the cell-ROIs

// Process selected cells
selectWindow(Image);
run("Gaussian Blur...", "sigma=1");

NCells = roiManager("count");
for (i = 0; i < NCells; i++) {
	selectWindow(Image);
	setSlice(ChDAPI);
	
	roiManager("select", i);
	Cellname[i] = call("ij.plugin.frame.RoiManager.getName", i);
	
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
	selectWindow(Cellname[i]);
	setSlice(2);
	run("Paste");

	// add 3rd channel (DAPI mask) to stack
	selectWindow(CellMask);
	roiManager("select", i);
	run("Copy");
	selectWindow(Cellname[i]);
	setSlice(3);
	run("Paste");
	run("Set...", "value=100000 slice"); // this is a pseudo mask for the DAPI channel

	// Actually measure number of foci
	setSlice(2);
	run("Select None"); // remove selection so that selection type is recognized properly
	run("Find Maxima...", "prominence="+prominence+" strict exclude output=[Point Selection]");
	
	if (selectionType() == -1) {
		nFoci = 0;
		mean = NaN;
		stdDev = NaN;
	} else {
		getSelectionCoordinates(x, y);
		nFoci = x.length;
	}

	// This part measures the size of a foci
	if (measure_fSize && nFoci > 0) {
		selectWindow(Cellname[i]);

		// get background
		setSlice(3);
		selectForeground(Cellname[i]);
		run("Create Selection");
		resetThreshold();
		setSlice(2);
//		BG = getPercentile(Cellname[i], 25);


		// Check the number of foci; If there's only one, particle segmentation fails.
		if (nFoci > 1) {
			run("Find Maxima...", "prominence="+prominence+" strict exclude output=[Segmented Particles]");
			run("Create Selection");
			close();
		} else {
			run("Restore Selection");
		}

		// Add another channel to saveImg
		selectWindow(Cellname[i]);
		setSlice(3);
		run("Copy");
		run("Add Slice", "add=slice");
		run("Paste");
		run("Restore Selection");

		// divide the foci with zero-value pixels;
		// This way, a part of the nucleus area can clearly be assigned to each foci.
		// Within this separated "puzzle piece", the area of a foci can then be measured
		if (nFoci > 1) {
			run("Clear", "slice");
		}
		
		// get area pieces: the ROIs of the pieces are added to the ROI manager (that already contains the cell ROIs)
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

		// loop over pieces;
		// The ROIs 0 - NCells belong to the cells.
		// The ROIs NCells + 1 - NROIs are puzzle pieces.
		for (j = NCells; j < roiManager("count"); j++) {
			processPiece(Cellname[i], j, 2, CutOff, Validation_mode);			
			run("Copy");
			close();
			selectWindow(Cellname[i]);
			run("Paste");
			resetThreshold();
			// store area of each foci
			Foci_sizes[j-NCells] = getResult("Area", nResults() - 1);
		}

		// measure mean size/std of all foci
		Array.getStatistics(Foci_sizes, min, max, mean, stdDev);

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

	// print averaged results to table
	print("[Measurements]", Cellname[i] +"\t" + 
							SizeNucleus +"\t"+ 
							nFoci+"\t" + 
							"Average" + "\t" + 
							mean +"\t"+ 
							stdDev);

	// print results for each foci in results file
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
	setThreshold(128, 100000);
	run("Create Selection");
}
