# load libraries
library("EBImage")

# read in the image
f <- system.file("images", "sample.png", package = "EBImage")
img <- readImage(f)

display(img, method = "raster")

nuc <- readImage(system.file("images", "nuclei.tif", package = "EBImage"))
display(nuc, method = "browser")

imgcol <- readImage(system.file("images", "sample-color.png", package = "EBImage"))

# print the image objects
print(img, short = TRUE)
print(imgcol, short = TRUE)

# examine the structure of the image
str(img)

# examine the pixel values
imageData(img)[1:3, 1:6]

imageData(imgcol)[1:3, 1:6, 1:3]
hist(imgcol)

# load image from folder
strawb <- readImage("/Users/mspriggs/Library/CloudStorage/OneDrive-Illumina,Inc./Documents/Applied_Bioinformatics/modules/machine_learning_for_metabolomics/wed/prac/Image Analysis Practical/strawb1.jpg")
display(strawb, method = "raster")

# convert to grayscale
strawb_gray <- channel(strawb, "gray")
display(strawb_gray, method = "raster")

# image adjustment
strawb_bright <- strawb_gray + 0.4
display(strawb_bright, method = "raster")

# make sure pixel values are between 0 and 1
strawb_bright[strawb_bright > 1] <- 1
strawb_bright[strawb_bright < 0] <- 0

display(strawb_bright, method = "raster")

strawb_invert <- normalize(1 - strawb_gray)
display(strawb_invert, method = "raster")

strawb_neg <- max(strawb_gray) - strawb_gray
display(strawb_neg, method = "raster")

strawb_contrast <- strawb_gray * 2
display(strawb_contrast)

strawb_exp <- strawb_gray^0.5
display(strawb_exp, method = "raster")

strab_comb <- combine(strawb_gray, strawb_gray + 0.4, strawb_gray * 2, strawb_gray^0.5)
display(strab_comb, method = "raster", all = TRUE)

writeImage(strawb_gray, "strawb_gray.jpg", quality = 90)

# applying linear filters
w <- makeBrush(size = 51, shape = "gaussian", sigma = 7)
strawb_smooth <- filter2(getFrame(strawb_gray, 1), w)
display(strawb_smooth, method = "raster")

nuc_gblur <- gblur(nuc, sigma = 5)
display(nuc_gblur, method = "raster", all = TRUE)

# edge detection
lap_filter <- matrix(1, nrow = 3, ncol = 3)
lap_filter[2, 2] <- -8
strawb_filt <- filter2(strab_comb, lap_filter)
display(strawb_filt, method = "raster", all = TRUE)

# image thresholding
strawb_thresh <- strawb_gray > 0.5
display(strawb_thresh, method = "raster")

# adaptive thresholding
disc <- makeBrush(7, shape = "disc")
disc <- disc / sum(disc)
offset <- 0.05
nuc_bg <- filter2(nuc, disc)
nuc_thresh <- nuc > (nuc_bg + offset)
display(nuc_thresh, method = "raster", all = TRUE)

display(thresh(nuc, w = 15, h = 15, offset = 0.05), all = TRUE)

# image segmentation
nmask <- watershed(distmap(nuc), 2)
display(colorLabels(nmask), method = "raster", all = TRUE)
