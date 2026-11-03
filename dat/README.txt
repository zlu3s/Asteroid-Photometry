This directory contains data for three targets and bias/flats in four folders.

GJ3069.raw is BD time series data. test.fits is junk, and the raw directory contains images straight off the detector.
Note: header transference was misbehaving tonight. If you are looking for a descriptive header, look in the raw folder.
All BD data was taken in the R band. A 15 second exposure was used.

niobe.raw is asteroid color and time series data. See the raw directory for good headers, including timestamps.
Filenames include the filter letter code used for that image. Data is niobe*.fits. A 15 second exposure was used.
My apologies for the troubles, it seems like the codebase malfunctions in some fun new way every time I'm up here.
0036 and 0060 are messed up by telescope motion.

SN.raw is supernova color images. Same deal as the first two. Data is of format SN*.fits, the file names have the
filter letter code used for that image. Exposures by band were RVBU 30/30/60/90, respectively. SNV0010 was messed up

rflat3 is the flatfield and bias images. There are one hundred biases, labeled r.ZERO*.fits. There are also a number
of flatfield images for each filter - 20 B (as inspection will show, the first 6 apparent b flats are saturated. 
Please double check all flat fields for saturation, as there was some trouble in obtaining them), 6 R (very bright),
26 U (rather faint), and 6 V (very bright).
For reasons unknown to me, the filter wheel had the B and V filters confused. I *think* the flats are labelled correctly,
but I am not absolutely positive. The headers in the raw database should have the correct filters, as they are read
directly from the telescope.

all the data is 3x3 binned, which means that badpixels.3x3.dat (a list of bad pixels) will work for all images.
the flats and bias frames will also work for all images, making sure to match with the filters of the flatfields.
data is less pre-processed than customary for this class, we just don't have the time on sky to set everything up neatly.

All data was taken on the Mont4k instrument. One image pixel is approximately 0.43''.