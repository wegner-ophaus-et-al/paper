function subtractBackground(radius) {
    id = getImageID();
    Stack.getDimensions(width, height, channels, slices, frames);

    for (c = 1; c <= channels; c++) {
        Stack.setPosition(c, 1, 1);
        run("Duplicate...", "title=bg_ch" + c);
        run("Median...", "radius=" + radius);
        selectImage(id);
        Stack.setPosition(c, 1, 1);
        imageCalculator("Subtract", id, "bg_ch" + c);
        close("bg_ch" + c);
    }
}



function adjustChannel(channel_number, name_var) {
	Stack.setChannel(channel_number);
	
	getStatistics(area, mean, min, max, std);
	
	nBins = 65536;
	getHistogram(values, counts, nBins);
	
	total = 0;
	for (i = 0; i < counts.length; i++)
	    total += counts[i];
	
	saturated_pixels = 0.05; //4;
	target = total * (1-(saturated_pixels/100));
	cumulative = 0;
	percentile_val = max; 
	
	for (i = 0; i < counts.length; i++) {
	    cumulative += counts[i];
	    if (cumulative >= target) {
	        percentile_val = i;
	        break;
	    }
	}
	
	setMinAndMax(0, percentile_val);
	saveAs("PNG", export_dir + File.separator + drug + File.separator + name_var +".png");
	
}

run("Set Scale...", "distance=1 known=0.1032 unit=µm global");
run("Scale Bar...", "width=5 height=2 thickness=6 font=0 bold overlay");
Stack.setDisplayMode("color");

export_dir = "/Volumes/Kur/paper_data_sorted/heatshock_stress-granule_markers/figures/representatives/images/F617";
drug = "heatshock_denoise";


File.makeDirectory(export_dir + File.separator + drug);

subtractBackground(120); 

//Set all channel colors
Stack.setChannel(1);
run("Yellow");
Stack.setChannel(3);
run("Cyan");
Stack.setChannel(2);
run("Magenta");

adjustChannel(1, "nls");
adjustChannel(3, "granulito");
adjustChannel(2, "stress_protein");

Stack.setDisplayMode("composite");
Stack.setActiveChannels("011");
saveAs("PNG", export_dir + File.separator + drug + File.separator + "merge" +".png");



