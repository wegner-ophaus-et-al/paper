root = "/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_gra/10hpf";

list = getFileList(root);
list = Array.sort(list);

setBatchMode(true);

	for (i = 0; i < list.length; i++) { //
		g = "processed_" + list[i];
		roiManager("Deselect");
		roiManager("Delete");
		processFile();
}

function generateMaskfromROI(roi_index, mask_name, output_dir, mask_width, mask_height) {
	newImage(mask_name, "8-bit black", mask_width, mask_height, 1);
	roiManager("Select", roi_index);
	run("Set...", "value=1");
	saveAs("Tiff", output_dir + File.separator + mask_name + ".tif");
	
}

function processFile() {
	working_dir = root + File.separator + list[i];
	list_tiff_files = getFileList(working_dir);
	tiff_file = "";
	tiff_file_suffix = ".tif";
	for (q = 0; q < list_tiff_files.length; q++) {
		if (endsWith(list_tiff_files[q], tiff_file_suffix)) {
			tiff_file = list_tiff_files[q];
		}
	}
	if (tiff_file == "") {
		continue
	}
	roiManager("Open", working_dir + File.separator + "csv" + File.separator + "RoiSet.zip");
//	print(working_dir + tiff_file);
	open(working_dir + File.separator + tiff_file);
	tif_width = getWidth();
	tif_height = getHeight();
	close();
	
	
	output = working_dir + "masks";
	File.makeDirectory(output);
	print(output);
	
	roi_count = roiManager("count");
	
	if (roi_count == 4) {
		
		// Irradiation area
		generateMaskfromROI(0, "irradiated", output, tif_width, tif_height);
		close();
		
		// background area
		generateMaskfromROI(2, "background", output, tif_width, tif_height);
		close();
		
		// correction area
		generateMaskfromROI(3, "correction", output, tif_width, tif_height);
		close();
	} else if (roi_count == 3) {
		// Irradiation area
		generateMaskfromROI(0, "irradiated", output, tif_width, tif_height);
		close();
		
		// background area
		generateMaskfromROI(1, "background", output, tif_width, tif_height);
		close();
		
		// correction area
		generateMaskfromROI(2, "correction", output, tif_width, tif_height);
		close();
	} else if (roi_count == 5) {
		// Irradiation area
		generateMaskfromROI(2, "irradiated", output, tif_width, tif_height);
		close();
		
		// background area
		generateMaskfromROI(3, "background", output, tif_width, tif_height);
		close();
		
		// correction area
		generateMaskfromROI(4, "correction", output, tif_width, tif_height);
		close();
	} else {
		exit("Unexpected number of ROIs: " + roi_count)
	}
}
setBatchMode(false);




//// function description
//
//roiManager("Delete");
//open("/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_dnd/10hpf/20220902_p2_dnd/csv/RoiSet.zip");
//roiManager("Open", "/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_dnd/10hpf/20220902_p2_dnd/csv/RoiSet.zip");
//open("/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_dnd/10hpf/20220902_p2_dnd/p2_dnd.tif");
//selectImage("p2_dnd.tif");
//// Get image dimensions
//close;
//
//newImage("bleach_mask", "8-bit black", 400, 400, 1);
//roiManager("Select", 0);
//run("Set...", "value=1");
//
//saveAs("Tiff", "/Volumes/HELHEIM/analyzed_data/diffusivity/FRAP_germ-granule_components/FRAP_dnd/10hpf/20220902_p2_dnd/masks/bleach_mask.tif");
//close;
//newImage("background_mask", "8-bit black", 400, 400, 1);
//roiManager("Select", 2);
////Same for the total_signal_mask
//roiManager("Select", 3);

