// ============================================================
// Segmentation Correction Macro
// Iterates over all sample folders in the data/ directory.
// For each folder:
//   1. Opens nls.tif + cell.tif mask  → user corrects → saves cell.tif
//   2. Opens nls.tif + nucleus.tif mask → user corrects → saves nucleus.tif
// ============================================================

// ---------- configuration ----------
//#@ File (label="Data directory", style="directory") dataDir
dataDir = "/Volumes/icb_remote/Documents/JW/py/projects/20260512_endoPGCs/data"

// -----------------------------------

setBatchMode(false);

dataPath = dataDir + File.separator;
folders  = getFileList(dataPath);

for (f = 0; f < folders.length; f++) {

    folderName = folders[f];

    // skip anything that is not a directory
    if (!File.isDirectory(dataPath + folderName)) continue;

    // strip trailing separator that getFileList may add
    sampleDir = dataPath + folderName;
    if (endsWith(sampleDir, File.separator))
        sampleDir = substring(sampleDir, 0, lengthOf(sampleDir) - 1);

    imagesDir = sampleDir + File.separator + "images" + File.separator;
    masksDir  = sampleDir + File.separator + "masks"  + File.separator;

    nlsPath    = imagesDir + "nls.tif";
    cellPath   = masksDir  + "cell.tif";
    nucleusPath= masksDir  + "nucleus.tif";

    // sanity checks – skip silently if files are missing
    if (!File.exists(nlsPath)) {
        print("Skipping " + folderName + " – nls.tif not found");
        continue;
    }

    // --------------------------------------------------------
    // ROUND 1 : nls.tif  +  cell.tif
    // --------------------------------------------------------
    correctOrDrawMask(nlsPath, cellPath, folderName, "cell", "yellow");

    // --------------------------------------------------------
    // ROUND 2 : nls.tif  +  nucleus.tif
    // --------------------------------------------------------
    correctOrDrawMask(nlsPath, nucleusPath, folderName, "nucleus", "cyan");

    print("Done with folder: " + folderName);
}

print("\n=== All folders processed ===");


// ============================================================
// Helper function – correct an existing mask OR let the user
// draw a new one when the mask file is absent or empty.
//
//   nlsPath    – full path to nls.tif (reference image)
//   maskPath   – full path to the mask tif to correct / create
//   folderName – human-readable label for dialogs
//   label      – "cell" or "nucleus"
//   colour     – overlay stroke colour string, e.g. "yellow"
// ============================================================
function correctOrDrawMask(nlsPath, maskPath, folderName, label, colour) {

    // ---- open reference image ----
    open(nlsPath);
    nlsID = getImageID();
    rename("nls [" + label + "] – " + folderName);
    run("Enhance Contrast", "saturated=0.35");
    getDimensions(w, h, ch, sl, fr);

    hasExistingSelection = false;

    // ---- try to load existing mask ----
    maskExists = File.exists(maskPath);
    maskEmpty  = true;   // assume empty until proven otherwise

    if (maskExists) {
        open(maskPath);
        maskID = getImageID();

        // check whether the mask actually contains any non-zero pixels
        getStatistics(area, mean);
        if (mean > 0) {
            maskEmpty = false;
            // convert mask pixels to a ROI selection
            setThreshold(1, 65535);
            run("Create Selection");
            roiManager("Add");
            hasExistingSelection = true;
        }
        close();   // close mask image – we work on nls from here
    }

    // ---- overlay existing ROI (if any) ----
    selectImage(nlsID);
    if (hasExistingSelection) {
        roiManager("Select", 0);
        Overlay.addSelection(colour, 2);   // add as coloured overlay
        roiManager("Select", 0);           // keep active so user can edit
    }

    // ---- prompt the user ----
    if (!maskExists || maskEmpty) {
        // no usable mask – ask user to draw from scratch
        waitForUser("Draw " + label + " segmentation",
            "Folder: " + folderName + "\n \n" +
            "No " + label + " segmentation was found (file missing or empty).\n" +
            "Please draw the " + label + " boundary using any selection tool\n" +
            "(Freehand, Polygon, Brush, Wand …), then click OK.");
    } else {
        // existing mask loaded – let user correct it
        waitForUser("Correct " + label + " segmentation",
            "Folder: " + folderName + "\n \n" +
            "The " + colour + " overlay shows the current " + label + " segmentation.\n" +
            "Edit the selection as needed (Brush, Wand, Polygon, etc.),\n" +
            "then click OK to save and continue.");
    }

    // ---- validate that a selection now exists ----
    selectImage(nlsID);
    selType = selectionType();   // -1 = no selection

    if (selType == -1) {
        // user dismissed without drawing – warn and skip saving
        showMessage("Warning",
            "No selection present for " + label + " in:\n" + folderName +
            "\n\nMask file will NOT be updated.");
        print("WARNING: no selection drawn for " + label + " in " + folderName + " – skipping save");
    } else {
        // ---- rasterise the (possibly new/edited) selection into a mask image ----
        newImage(label + "_corrected", "16-bit black", w, h, 1);
        corrID = getImageID();
        run("Restore Selection");      // bring the selection from nls image
        run("Fill", "slice");          // fill with 65535 (white)
        run("Select None");
        resetThreshold();

        saveAs("Tiff", maskPath);
        print("Saved " + label + " mask: " + maskPath);
        close();   // close corrected mask image
    }

    // ---- clean up ----
    if (roiManager("count") > 0) {
        roiManager("Deselect");
        roiManager("Delete");
    }
    selectImage(nlsID);
    run("Remove Overlay");
    run("Select None");
    close();
}
