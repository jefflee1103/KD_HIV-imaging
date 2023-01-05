// Get directory of the image
outDir = getDirectory("image");

// Get image name string and create output path
name = getTitle;
basename = substring(name,0,lengthOf(name)-5)
outpath = outDir+basename

// Convert ROIs to a Label image
run("ROIs to Label image");

// Save the Label image 
saveAs("Tiff", outpath+"_mask.tiff");

// Remove ROIs
roiManager("Delete");

// Close all windows
run("Close All");


