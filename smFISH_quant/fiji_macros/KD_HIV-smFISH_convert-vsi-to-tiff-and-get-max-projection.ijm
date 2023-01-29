// Define input and output directories
input = getDirectory("Choose Source Directory ");
output = getDirectory("Choose Destination Directory ");

// Open files
list = getFileList(input);

setBatchMode(true);

for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	path = input + list[i];
	if(endsWith(list[i], ".vsi")) {
    	run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
    	fileName = substring(list[i],0,lengthOf(list[i])-4);
    	saveAs("Tiff", output+fileName);
		run("Z Project...", "projection=[Max Intensity]");
		saveAs("Tiff", output+fileName+"_MAX.tif");
		run("Close All");
	}
}

setBatchMode(false)