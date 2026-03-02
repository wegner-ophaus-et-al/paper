root = ""
list = getFileList(root);
list = Array.sort(list);
for (i = 0; i < list.length; i++) {
    sample_name = list[i];
    // Print index in case, to continue in case of error
    print(i);
    input = root + File.separator + sample_name;
    open(input + File.separator + "raw" + File.separator + "granules.tif");
    alfaID = getImageID(); // Store the image ID of alfa.tif
    getDimensions(width, height, _, _, _);
    run("Enhance Contrast...", "saturated=30");
    run("Duplicate...", "title=[Segmentation Helper]");
    run("Gaussian Blur...", "sigma=5");
	run("Find Edges");
    open(input + File.separator + "raw" + File.separator + "nanos-rna.tif");
    open(input + File.separator + "segmentation" + File.separator + "granules.tif");
    granule_segID = getImageID();
    //open(input + File.separator + "segmentation" + File.separator + "granules.tif");
    
    // Tile windows
    run("Tile");
    
    setTool("wand");
    waitForUser("Select unwanted segmentations");
	// Check if there is a selection before clearing
    if (selectionType() != -1) {
        run("Clear", "slice");
    }
    run("Select None");
    
//    // Set multipoint tool and wait for user to select best granule
//    setTool("multipoint");
//    waitForUser("Select best situated granule");
//    
//    // Get the coordinates from the multipoint selection
//    getSelectionCoordinates(xpoints, ypoints);
//    
//    // Save coordinates to text file
//    if (xpoints.length > 0) {
//        coordFile = input + File.separator + "granule_coords.txt";
//        f = File.open(coordFile);
//        print(f, "X\tY");
//        for (j = 0; j < xpoints.length; j++) {
//            print(f, xpoints[j] + "\t" + ypoints[j]);
//        }
//        File.close(f);
//    }
//    
    selectImage(granule_segID);
    saveAs("tiff", input + File.separator + "segmentation" + File.separator + "granules.tif");
    close();
    
//    //Load and refine the nucleus segmentation
//    open(input + File.separator + "segmentation" + File.separator + "nucleus.tif");
//    nuc_segID = getImageID();
//    run("Tile");
//    waitForUser("Refine nucleus segmentation");
//    
//    if (selectionType() != -1) {
//        run("Set...", "value=1");
//    }
//    selectImage(nuc_segID);
//    saveAs("tiff", input + File.separator + "segmentation" + File.separator + "nucleus.tif");
//    close();
//    
    //Load and refine the cell segmentation
    
//    open(input + File.separator + "segmentation" + File.separator + "cell.tif");
    newImage("cell", "8-bit black", width, height, 1);
    setMinAndMax(0, 1);
    cell_segID = getImageID();
    run("Tile");
    waitForUser("Refine cell segmentation");
    if (selectionType() != -1) {
		run("Set...", "value=1");
    }
    selectImage(cell_segID);
    saveAs("tiff", input + File.separator + "segmentation" + File.separator + "cell.tif");
    close();
    
//    // Select the alfa.tif window using its image ID
//    selectImage(alfaID);
//    makeRectangle(0, 0, 100, 150);
//    setTool("rectangle");
//    run("Enhance Contrast...", "saturated=30");
//    waitForUser("Select a noise region for SNR calculations ");
//    run("Measure");
//	saveAs("Results", input + File.separator + "noise.csv");
//    
//    
//    
//    run("Clear Results");
    close("*");
}
setBatchMode(false);
