root = "/Volumes/HELHEIM/analyzed_data/heat_shock/raw_data_KT/20260114_RNAscope_bactin-561_Nos-633_Vasa-488_Tg24_12hpf_heat-shock_vs_ctrl"
list = getFileList(root);
list = Array.sort(list);
run("Close All");
for (i = 53; i < list.length; i++) {
    sample_name = list[i];
    // Print index in case, to continue in case of error
    print(i);
    input = root + File.separator + sample_name;
    open(input + File.separator + "raw" + File.separator + "nucleus.tiff");
    open(input + File.separator + "raw" + File.separator + "granule.tiff");
    
    getDimensions(width, height, _, _, _);
    run("Enhance Contrast...", "saturated=1");
    run("Duplicate...", "title=[Segmentation Helper]");
    run("Enhance Contrast...", "saturated=5");
    run("Gaussian Blur...", "sigma=5");
	run("Find Edges");
	resetMinAndMax;
//    open(input + File.separator + "raw" + File.separator + "nanos-rna.tif");
    
    open(input + File.separator + "segmentations" + File.separator + "granules.tiff");
    granule_segID = getImageID();
    setMinAndMax(0, 1);
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
    setTool("multipoint");
    waitForUser("Select best situated granule");
    
    // Get the coordinates from the multipoint selection
    getSelectionCoordinates(xpoints, ypoints);
    
    // Save coordinates to text file
    if (xpoints.length > 0) {
        coordFile = input + File.separator + "granule_coords.txt";
//        File.close(f);
        f = File.open(coordFile);
        
        print(f, "X\tY");
        for (j = 0; j < xpoints.length; j++) {
            print(f, xpoints[j] + "\t" + ypoints[j]);
        }
        File.close(f);
        print(f);
    }
    run("Clear Results");
   
     
    
    
    
    selectImage(granule_segID);
    saveAs("tiff", input + File.separator + "segmentations" + File.separator + "granule.tiff");
    close();
    
    //Load and refine the nucleus segmentation
    open(input + File.separator + "segmentations" + File.separator + "nucleus.tiff");
    nuc_segID = getImageID();
    run("Tile");
    waitForUser("Refine nucleus segmentation");
    
    if (selectionType() != -1) {
        run("Set...", "value=1");
    }
    selectImage(nuc_segID);
    saveAs("tiff", input + File.separator + "segmentations" + File.separator + "nucleus.tiff");
    close();
    
    //Load and refine the cell segmentation
    
//    open(input + File.separator + "segmentation" + File.separator + "cell.tif");
	newImage("cell", "16-bit black", width, height, 1);
    setMinAndMax(0, 1);
    cell_segID = getImageID();
    run("Tile");
    waitForUser("Refine cell segmentation");
    if (selectionType() != -1) {
		run("Set...", "value=1");
    }
    selectImage(cell_segID);
    saveAs("tiff", input + File.separator + "segmentations" + File.separator + "cell.tiff");
    close();
    

    close("*");
}
setBatchMode(false);