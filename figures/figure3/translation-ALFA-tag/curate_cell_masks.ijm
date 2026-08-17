// =====================================================================
//  curate_cell_masks.ijm
//
//  Walks through a folder of samples
//
//      <root>/<8hexdigits>__<realname>/original_raw/<realname>.lsm
//      <root>/<8hexdigits>__<realname>/masks/cell.tif
//
//  For every sample it
//    1. opens the .lsm and shows channel CHANNEL,
//    2. loads masks/cell.tif and puts it on the image as a selection,
//    3. activates the selection brush tool and waits,
//    4. on OK: turns the (edited) selection back into a binary mask and
//       overwrites masks/cell.tif  (a one-time .bak copy is kept).
//
//  Only the first ID_CHARS characters of the folder name are ever shown
//  (window title, dialog title, log) so the blinding is preserved.
// =====================================================================

// ------------------------- settings ---------------------------------
CHANNEL        = 3;              // channel to display
RAW_SUBDIR     = "original_raw";
RAW_EXT        = ".lsm";
MASK_SUBDIR    = "masks";
MASK_NAME      = "cell.tif";
ID_CHARS       = 8;              // how many characters of the name to show
START_AT       = 1;              // resume: 1 = first sample, 5 = fifth, ...
USE_BIOFORMATS = true;           // false -> plain open()
PROJECT_Z      = false;          // true  -> max-project z-stacks first
AUTO_CONTRAST  = true;
MAKE_BACKUP    = true;           // keep cell.tif.bak of the original mask
BRUSH_SIZE     = 15;             // 0 = leave the current brush size alone
// --------------------------------------------------------------------

root = getDirectory("Choose the folder that contains the sample folders");

setOption("BlackBackground", true);
run("Colors...", "foreground=white background=black selection=yellow");
run("Close All");
print("\\Clear");

// ---- collect sample folders ----
entries = getFileList(root);
Array.sort(entries);
dirs = newArray(0);
for (i = 0; i < entries.length; i++) {
    if (File.isDirectory(root + entries[i])) {
        nm = replace(entries[i], "\\\\", "/");
        if (endsWith(nm, "/")) nm = substring(nm, 0, lengthOf(nm) - 1);
        if (nm != "." && nm != "..") dirs = Array.concat(dirs, nm);
    }
}
print("Found " + dirs.length + " sample folders in " + root);

quit = false;
for (k = START_AT - 1; k < dirs.length && !quit; k++) {

    id      = dirs[k];
    shortID = id;
    if (lengthOf(shortID) > ID_CHARS) shortID = substring(shortID, 0, ID_CHARS);

    rawDir   = root + id + File.separator + RAW_SUBDIR + File.separator;
    maskDir  = root + id + File.separator + MASK_SUBDIR + File.separator;
    maskPath = maskDir + MASK_NAME;

    // ---- find the raw file ----
    raw = "";
    if (File.isDirectory(rawDir)) {
        fl = getFileList(rawDir);
        for (j = 0; j < fl.length; j++)
            if (raw == "" && endsWith(toLowerCase(fl[j]), RAW_EXT)) raw = rawDir + fl[j];
    }

    if (raw == "") {
        print(shortID + "   ->  no " + RAW_EXT + " found, skipped");
    } else {
        showProgress(k, dirs.length);

        // ---- open raw, immediately hide the real file name ----
        if (USE_BIOFORMATS)
            run("Bio-Formats Importer",
                "open=[" + raw + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
        else
            open(raw);
        imgID = getImageID();
        rename(shortID);                       // <-- blinding

        if (PROJECT_Z) {
            getDimensions(w, h, ch, sl, fr);
            if (sl > 1) {
                run("Z Project...", "projection=[Max Intensity]");
                projID = getImageID();
                selectImage(imgID); close();
                imgID = projID;
                selectImage(imgID);
                rename(shortID);
            }
        }

        getDimensions(w, h, ch, sl, fr);
        if (ch > 1) {
            Stack.setDisplayMode("color");
            Stack.setChannel(minOf(CHANNEL, ch));
        }
        if (AUTO_CONTRAST) run("Enhance Contrast", "saturated=0.35");   // before the ROI!

        // ---- existing mask -> selection ----
        hasRoi = false;
        if (File.exists(maskPath)) {
            open(maskPath);
            mID = getImageID();
            getDimensions(mw, mh, mc, ms, mf);
            if (mw != w || mh != h)
                print(shortID + "   ->  WARNING: mask is " + mw + "x" + mh +
                      ", image is " + w + "x" + h);
            if (bitDepth() == 8) setThreshold(1, 255); else setThreshold(1, 65535);
            run("Create Selection");
            resetThreshold();
            if (selectionType() != -1) {
                roiManager("reset");
                roiManager("add");
                hasRoi = true;
            }
            selectImage(mID); close();
        } else {
            print(shortID + "   ->  no mask yet, starting from scratch");
        }

        selectImage(imgID);
        if (hasRoi) {
            roiManager("select", 0);
            roiManager("reset");
            if (ch > 1) Stack.setChannel(minOf(CHANNEL, ch));   // in case the ROI moved us
        }
        if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
        selectImage(imgID);

        // ---- hand over the selection brush ----
        setTool("brush");
        // If this line throws an error in your installation, set BRUSH_SIZE = 0
        // and set the size by double-clicking the tool icon instead.
        if (BRUSH_SIZE > 0) eval("script", "ij.gui.Toolbar.setBrushSize(" + BRUSH_SIZE + ");");

        Dialog.createNonBlocking("Sample " + shortID + "   (" + (k + 1) + "/" + dirs.length + ")");
        Dialog.addMessage("Edit the cell mask with the selection brush:\n" +
                          "    drag           = add to selection\n" +
                          "    Alt + drag     = remove from selection\n" +
                          "    double-click tool icon = brush size\n \n" +
                          "OK = write mask,   Cancel = abort the macro.");
        Dialog.addCheckbox("Save mask for this sample", true);
        Dialog.addCheckbox("Stop after this sample", false);
        Dialog.show();
        doSave = Dialog.getCheckbox();
        quit   = Dialog.getCheckbox();

        // ---- selection -> binary mask -> overwrite file ----
        if (doSave) {
            selectImage(imgID);
            getDimensions(w, h, ch, sl, fr);

            if (!File.isDirectory(maskDir)) File.makeDirectory(maskDir);
            if (MAKE_BACKUP && File.exists(maskPath) && !File.exists(maskPath + ".bak"))
                File.copy(maskPath, maskPath + ".bak");

            if (selectionType() == -1) {
                newImage("newmask", "8-bit black", w, h, 1);        // empty mask
                print(shortID + "   ->  no selection, wrote EMPTY mask");
            } else {
                run("Create Mask");        // 8-bit, 255 inside the selection
                run("Grays");              // drop the inverting LUT
                run("Select None");
            }
            saveAs("Tiff", maskPath);
            close();
            print(shortID + "   ->  mask saved");
        } else {
            print(shortID + "   ->  skipped, file unchanged");
        }

        run("Close All");
    }
}

run("Close All");
showProgress(1);
print("Done.");
