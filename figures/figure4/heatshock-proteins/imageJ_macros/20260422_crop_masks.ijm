input = "/Users/julian/local_files/20260505_heat-shock/F618"


list = getFileList(input);
list = Array.sort(list);

	for (i = 0; i < list.length; i++) {
		g = list[i];
		is_directory = File.isDirectory(input + File.separator + g);
		if (is_directory == true) {	
		
			raw_file_dir = input + File.separator + g + File.separator + "original_raw";
			raw_file_dir_filelist = getFileList(raw_file_dir);
//			print(raw_file_dir_filelist);
			for (j = 0; j < raw_file_dir_filelist.length ; j++) {
				if (File.exists(raw_file_dir + File.separator + raw_file_dir_filelist[j]) == true && startsWith(raw_file_dir_filelist[j], ".") == false) {
//					print("Processing: " + raw_file_dir_filelist[j]);
					tif_file_name = raw_file_dir_filelist[j];
					processFile();	
				}

			}
//		File.exists(path);
//		if (endsWith(list[i], ".tif")){
//		processFile();
		}
		
}

function processFile() {
	open(raw_file_dir + File.separator + tif_file_name);
	og_img_id = getImageID();
	og_name = getInfo("image.filename");
	output_path = input + File.separator + g + "masks";
	File.makeDirectory(output_path);
	
	// MIP that sh*t
	run("Z Project...", "projection=[Max Intensity]");
	mip_image_id = getImageID();
	mip_name = getInfo("image.filename");
	run("Tile");
	print("index: ", i + "/"+ list.length);
	window_name_array = newArray(og_img_id, mip_image_id);
	for (p = 0; p < window_name_array.length; p++) {
		// Set exposure
		selectImage(window_name_array[p]);
		Stack.setDisplayMode("composite");
		run("Grays");
		run("Enhance Contrast", "saturated=0.01");
		run("Next Slice [>]");
		run("Cyan");
		run("Enhance Contrast", "saturated=0.00001");
		getMinAndMax(_, max);
		setMinAndMax(250, max);
		run("Next Slice [>]");
		run("Magenta");
		run("Enhance Contrast", "saturated=0.00001");
		getMinAndMax(_, max);
		setMinAndMax(250, max);
		run("Previous Slice [<]");
	}
	
	selectImage(mip_image_id);
	breaker = false;
	q = 0;
	while (breaker == false) {
		selectImage(mip_image_id);
		makeRectangle(100, 100, 256, 256);
		
		// Get screen width
		sw = screenWidth;
		
		// Ask for mermaid
		choice_array = newArray("Add", "Next Image");
		Dialog.createNonBlocking("Add this mask?");
		Dialog.enableYesNoCancel("Yes", "Next image");
		Dialog.setLocation(sw/4,0);
		Dialog.show();
		choice = Dialog.getYesNoCancel();
		

		if (choice == "yes") {
		run("Create Mask");
		new_string = "%" + q + ".tif";
		print(output_path + File.separator + replace(og_name, ".tif", new_string));
		saveAs("Tiff", output_path + File.separator + replace(og_name, ".tif", new_string));
		close;
		q += 1;
		}
		else {
			breaker = true;
			close("*");
			}
	}
	

	
	
}
