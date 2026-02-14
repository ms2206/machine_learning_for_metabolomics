### Machine Learning for Metabolomics ####

### Image Analysis Practical - Answers

### Maria Anastasiadi



# Install EBImage from Bioconductor
if (!require("BiocManager", quietly = TRUE, version=3.20))
  install.packages("BiocManager")

BiocManager::install("EBImage", force = TRUE)

browseVignettes("EBImage")


library("EBImage")


f = system.file("images", "sample.png", package="EBImage")
img = readImage(f)
print(img, short=TRUE)
dim(img)
str(img)
imageData(img)[1:5, 1:6]
hist(img)

# Show coloured image
imgcol = readImage(system.file("images", "sample-color.png", package="EBImage"))
display(imgcol)
print(imgcol, short=TRUE)
dim(imgcol)
str(imgcol)
imageData(imgcol)[1:5, 1:6, 2]
hist(imgcol)

#Display the loaded image
EBImage::display(img)
display(img, method="raster")
text(x = 20, y = 20, label = "Parrots", adj = c(0,1), col = "orange", cex = 2)

## Load image from file
strawb <- readImage("strawb1.jpg")
display(strawb, method="raster")
display(strawb, method="browser")

strawb.gray <-channel(strawb, mode="gray")
display(strawb.gray)

# Display all images 
nuc = readImage(system.file("images", "nuclei.tif", package = "EBImage"))
EBImage::display(nuc, method = "raster", all = TRUE)
EBImage::display(1 - nuc, all = TRUE)

# Display a single image
EBImage::display(1 - nuc, frame = 2)


# Image manipulation

# Produce negative image
strawb_neg = max(strawb.gray) - strawb.gray
display( strawb_neg )


# Flip dark areas to light and vice versa
strawb.inv <- normalize(-strawb.gray)
display(strawb.inv)

# Increase the brightness
strawb_bright = strawb.gray +0.4
# strawb_bright <- normalize(strawb_bright) If you do this it will result in the original picture
display( strawb_bright )

strawb_bright[strawb_bright > 1] <- 1 # Replace values more than 1 with 1
display( strawb_bright )

strawb_contr = strawb.gray * 1.5
strawb_contr <- normalize(strawb_contr)
display( strawb_contr)
strawb_contr[strawb_contr > 1] <- 1 # Replace values more than 1 with 1
display( strawb_contr )


# adjust the contrast and the gamma factor through multiplication and exponentiation
strawb_exp = strawb.gray ^ 0.5
display(strawb_exp)

# combine
strawb_comb = combine(
  strawb.gray, # original gray image
  strawb.gray + 0.3, # brighter image
  strawb.gray * 2, # increased contrast
  strawb.gray ^ 0.5 # gamma correction
)
display(strawb_comb, all=TRUE)

## Save strawberry image in grayscale
writeImage(strawb.gray, "strawb1_gray.jpeg", quality = 85)

## Gaussian (linear) filter
w = makeBrush(size = 51, shape = "gaussian", sigma = 7)
strawb.smooth = filter2(getFrame(strawb.gray, 1), w)
display(strawb.smooth)


## In one go
nuc_gblur = gblur(nuc, sigma = 5)
display(nuc_gblur, all=TRUE )

## Edge detection with Laplacian filter
lap.filt = matrix(1, nrow = 3, ncol = 3)
lap.filt[2, 2] = -8
strawb_filt = filter2(strawb_comb, lap.filt)
display(strawb_filt)

display(strawb_comb, all=TRUE)

## The increased contrast transformation results in best edge detection

writeImage(strawb_filt, "strawb1_Lap_Filt.jpeg", quality = 100)

## Global thresholding
nuc_thresh = nuc >0.05
display(nuc_thresh, all=TRUE)

## NOTE: If you want to add labels to images with multiple panels you can do it like so:

nx = 2
ny = 3
width = dim(nuc_thresh)[1]
height = dim(nuc_thresh)[2]
x_offset = y_offset = 20
n = numberOfFrames(nuc_thresh, 'render')

## actual plotting
display(nuc_thresh, method = "raster", all = TRUE, nx = 2)
text(x = rep(seq(from = 0, by = width, length.out = nx), ny) + x_offset,
     y = rep(seq(from = 0, by = height, length.out = ny), each = nx) + y_offset,
     label = LETTERS[1:n], 
     adj = c(0,1), col = "orange", cex = 2)

## Apply thresholding in the strawberry image
strawb_thresh = strawb.gray > 0.5
display(strawb_thresh)

## Adaptive Thresholding
disc = makeBrush(31, "disc")
disc = disc / sum(disc)
offset = 0.05
nuc_bg = filter2( nuc, disc )
nuc_th = nuc > nuc_bg + offset
display(nuc_th, all=TRUE)


strawb_bg = filter2( strawb.gray, disc )
strawb_th = strawb.gray > strawb_bg + offset
display(strawb_th)

## Image segmentation
nmask = watershed( distmap(nuc_th), 2 )
display(colorLabels(nmask), all=TRUE)



