root = "/Users/julian/local_files/20260428_stress-granules_HS/F617"
list = getFileList(root);
list = Array.sort(list);
run("Close All");
for (i = 0; i < list.length; i++) {
    sample_name = list[i];
    // Print index in case, to continue in case of error
    print("Index: " + i + "/" + list.length);
    input = root + File.separator + sample_name;
    open(input + File.separator + "images" + File.separator + "stress.tif");
    open(input + File.separator + "images" + File.separator + "nls.tif");
    
    getDimensions(width, height, _, _, _);
    run("Enhance Contrast...", "saturated=1");
    run("Duplicate...", "title=[Segmentation Helper]");
    run("Enhance Contrast...", "saturated=5");
    run("Gaussian Blur...", "sigma=5");
	run("Find Edges");
	resetMinAndMax;
      

	newImage("cell", "8-bit black", width, height, 1);
    setMinAndMax(0, 1);
    cell_segID = getImageID();
    run("Tile");
    setTool("brush");
    waitForUser("Refine cell segmentation");
    if (selectionType() != -1) {
		
		roiManager("Add");
		count_roimanager = roiManager("count");
		selectImage(cell_segID);
		roiManager("Select", count_roimanager -1);
		run("Set...", "value=1");
    }
    selectImage(cell_segID);
    
    maskDir = input + File.separator + "masks" + File.separator;
    if (!File.isDirectory(maskDir)) {
    File.makeDirectory(maskDir);
    
	}

    saveAs("tiff", maskDir + "cell.tif");
    
    close();
    

    close("*");
}
setBatchMode(false);