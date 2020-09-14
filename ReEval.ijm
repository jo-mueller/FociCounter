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


////////////////////////////////////////START//////////////////////////////////////////////////////

// Input GUI
#@ File (label="Input Directory", style="directory") DataDir

#@ Boolean (label="Save (new) images?",  saveImgs=saveImgs) saveImgs

#@ Boolean (label="measure foci size",  measure_fSize=measure_fSize) measure_fSize
#@ Integer (label="Prominence", min=0, max=10000, value=prominence) prominence
#@ Float (label="Cutoff", min=0, max=1.0, value=CutOff) CutOff


// Clean up
close("*");
roiManager("reset");
run("Clear Results");

// Create Results table
if (isOpen("Measurements")) {
	close("Measurements");	
}

run("Table...", "name=[Measurements] width=800 height=600");
print("[Measurements]", "\\Headings:Label \tArea \tN_Foci \tType \tFSize_Mean \tFSize_Std \tMutualdistance");

// Allocate variables and set measurements
run("Set Measurements...", "area min median center display redirect=None decimal=2");

// Create directory for new results
Sample = split(DataDir, "\\");
Sample = Sample[Sample.length - 1];
NewOutDir = DataDir + "\\" + Sample + 
			"_ReEval_CutOff" + d2s(CutOff*10, 0) + "_prominence" + d2s(prominence, 0) + "\\";
if (!File.exists(NewOutDir)) {
	File.makeDirectory(NewOutDir);
}

// get analyzed cell images
Filelist = getFileList(DataDir);

for (i = 0; i < Filelist.length; i++) {
	// skip non-images
	if (!endsWith(Filelist[i], "tif")) {
		continue;
	}

	open(DataDir + "\\" + Filelist[i]);
	CellName = File.nameWithoutExtension;
	rename(CellName);

	// Measure the Nucleus size
	setSlice(3);
	run("Measure");
	SizeNucleus = getResult("Area", nResults()-1);


	// Delete last slice
	setSlice(4);
	run("Delete Slice");

	// Actually measure number of foci
	// The following part of the script is basically the re-evaluation (to be written)
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
	if (nFoci > 0) {
		/* 
		 *  Add another channel to output image.
		 *	Copy-paste segmented nucleus there
		 *	Then project the foci-wise regions into this new slice */
		selectWindow(CellName);
		setSlice(3);
		run("Copy");
		run("Add Slice", "add=slice");
		run("Paste");
		run("Select None");
		
		// Divide image in ROIs that contain one foci each.
		// If there's only one foci, particle segmentation fails and the
		// entire nucleus is considered as a single area of interest.
		if (nFoci > 1) {
			setSlice(2);
			run("Find Maxima...", "prominence=" + prominence + " strict exclude output=[Segmented Particles]");
			setThreshold(128, 255);  // make sure correct area is selected
			run("Create Selection");
			SegMap = getTitle();
		} else {
			run("Restore Selection");
		}

		selectWindow(CellName);
		run("Restore Selection");
		close(SegMap);
		

		/*
		*	divide the binary nucleus area into separated regions with zero-value pixels;
		*	This way, a part of the nucleus area can clearly be assigned to each foci.
		*	Within this separated "puzzle piece", the area of a foci can then be measured
		*/
		if (nFoci > 1) {
			setSlice(4);
			run("Clear Outside", "slice");
		}

		selectWindow(CellName);
		selectForeground(CellName);
		roiManager("Add");
		roiManager("Select", roiManager("count")-1);
		if (nFoci > 1) {
			roiManager("Split");
			roiManager("Select", 0); // delete union of all regions (index = 0)
			roiManager("Delete");
		}
		
		Foci_sizes = newArray(nFoci);
		Foci_x = newArray(nFoci);
		Foci_y = newArray(nFoci);

		/*
		 * Iterate over all foci regions (a.k.a. puzzle pieces)
		 * Measure Area and other stuff of each foci
		*/
		for (j = 0; j < roiManager("count"); j++) {
			processPiece(CellName, j, 2, CutOff);	

			// Copy Foci segmentation into Cell-Image
			run("Copy");
			close();
			selectWindow(CellName);
			run("Paste");
			resetThreshold();
			
			// store area of each foci
			Foci_sizes[j] = getResult("Area", nResults() - 1);
			Foci_x[j] = getResult("XM", nResults() - 1);
			Foci_y[j] = getResult("YM", nResults() - 1);
		}

		// measure mean size/std of all foci
		Array.getStatistics(Foci_sizes, min, max, mean, stdDev);

		// Median distance of foci to its neighbours
		MutDist = MutualDistance(Foci_x, Foci_y);

		// Clean up RoiManager
		roiManager("reset");

		// do the point selection again so that images can be stored with the point selection.
		// May be easier for later introspection.
		resetThreshold();
		run("Select None");
		setSlice(2);
		run("Find Maxima...", "prominence="+prominence+" strict exclude output=[Point Selection]");
	}
	
	// If images of each cell should be stored (no matter how many foci).
	if (saveImgs){
		saveAs(".tif", NewOutDir + CellName);
	}


	// print cell-averaged results to table
	print("[Measurements]", CellName +"\t" + 
							SizeNucleus +"\t"+ 
							nFoci+"\t" + 
							"Average" + "\t" + 
							mean +"\t"+ 
							stdDev + "\t" + 
							MutDist);

	// print foci-wise results in results to table
	for (f = 0; f < Foci_sizes.length; f++) {
		print("[Measurements]", CellName +"\t" +
								" " +"\t" +
								" " +"\t" +
								"Single" + "\t" +
								Foci_sizes[f] +"\t"+
								" ");
	}
}

selectWindow("Measurements");
saveAs("Measurements", NewOutDir + Sample + "_ReEval.csv");


function processPiece(Image, ROI, Slice, p){
	// this function looks at a part of a cell that contains one (and only one) Foci
	// and measures its area.
	// Input:
	// Image - name of the image to process.
	// ROI - ROI of this piece
	// Slice - If the Image has more than one slice, specify it. Should be the slice with the intensity information.
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

	setThreshold(BG + p * (Max - BG), Max);
	run("Convert to Mask", "background=Dark black");
	selectForeground("PuzzlePiece");
	run("Measure");
	run("16-bit");
	run("Multiply...", "value=1000 slice");
	run("Select None");
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


function MutualDistance(X, Y){
	/*
	 * This function takes Vectors <x> and <y> with the center of mass
	 * coordinates of measured foci. The euclidian distance of each foci
	 * to all other foci is then measured. This gives and estimate
	 * on how densely foci are present in the analyzed cell.
	 */

	N = X.length;
	output = newArray();
	
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
	Array.show(output);
	return output[floor(output.length/2)];
}

function selectForeground (Input){
	selectWindow(Input);
	setThreshold(128, 100000);
	run("Create Selection");
}

