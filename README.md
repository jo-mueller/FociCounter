# FociCounter
Performs auto-segmentation of cell nuclei and measures the number and size of radiation induced DNA lessions within the nucleus.

## Usage
In order to use this script, download [Fiji](https://imagej.net/software/fiji/downloads) and open the script by drag & dropping the file `FociEval_V2020_2_JM.ijm` onto the Fiji toolbar. Click `Run` to start the script. 

A user interface will then pop up:
<img src="./docs/imgs/GUI.PNG" width="50%" height="50%">]

Parameters to be set and what they do:

* `Input image`: file path to a n-channel image, out of which one should be a DAPI staining and the other one showing the DNA lessions ("foci").
* `Foci prominence`: Determines how sensitive the Foci detection should be. High numbers lead to fewer detected foci, low numbers yield more foci detections. It is highly encouraged to optimize this parameter before using this script. To do so, use the `Process > Find Maxima > ..` function from the Fiji toolbar, set the output type to `Point Selection` and check the `Preview point selection` box. ou can then try different prominence value. 
* `Intensity cutoff`: Determines which part around a detected foci will be counted as part of the foci and which is not. If set to `0.5` (which is suggested), this parameter corresponds to a [full-width at half maximum](https://en.wikipedia.org/wiki/Full_width_at_half_maximum) area criterion.
* `Pixel size`: Check if correct - to validate, open your image and press `Ctrl + I` to inspect the metadata of your image.
* `gH2AX` channel: Number of slice in the stack of channels that contains the image with the foci.
* `DAPI` channel: Number of slice in the stack of channels that contains the DAPI staining. (*Note*: Start counting at 1)
* `Minimal nucleus size`: Exclude nuclei smaller than this size from the analysis
* `Foci brightness (lower percentile)`: The brightness of all foci in the image is measured and nuclei with a brightness lower than this value will be excluded from the analysis
* `Foci brightness (upper percentile)`: The brightness of all foci in the image is measured and nuclei with a brightness higher than this value will be excluded from the analysis. It makes sense to use this parameter to exclude, i.e., dividing nuclei from the analysis.
* `Boundary exclusion radius`: Cells this close to the edge of the image will be excluded from the analysis.
* `Batch mode`: Hides the processing in the background if checked - also makes the analysis faster.
* `Save cell wise images`: Saves a copy of the DAPI staining, the foci staining and the derived segmentations of the nucleus and the foci in a separate image for introspection.
