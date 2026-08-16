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



function adjustChannel(channel_number, name_var, saturated_pixels) {
	Stack.setChannel(channel_number);
	
	getStatistics(area, mean, min, max, std);
	
	nBins = 65536;
	getHistogram(values, counts, nBins);
	
	total = 0;
	for (i = 0; i < counts.length; i++)
	    total += counts[i];
	
//	saturated_pixels = 4;//4;
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
	print(name_var + ": " +"0, " + percentile_val )
	saveAs("PNG", export_dir + File.separator + drug + File.separator + name_var +".png");
	
}

run("Scale Bar...", "width=5 height=2 thickness=6 font=0 bold overlay");
print(File.name);
Stack.setDisplayMode("color");

Stack.setChannel(1);
run("Cyan");
Stack.setChannel(4);
run("Magenta");


export_dir = "/Volumes/Kur/paper_data_sorted/ALFAtag/19xALFAtag_confo_KT/representatives";
drug = "a1dd8afd_patA";

File.makeDirectory(export_dir + File.separator + drug);

//subtractBackground(50); 

adjustChannel(1, "granules", 0.008);
adjustChannel(4, "alfa", 9);
//adjustChannel(3, "vasa_rna");

Stack.setDisplayMode("composite");
Stack.setActiveChannels("1001");
saveAs("PNG", export_dir + File.separator + drug + File.separator + "merge" +".png");



