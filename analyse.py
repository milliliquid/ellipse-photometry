import numpy as np
import matplotlib.pyplot as plt
from photutils.aperture import EllipticalAperture
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from photutils.isophote import EllipseSample, EllipseFitter
from photutils.aperture import EllipticalAperture
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from photutils.segmentation import (detect_sources)
from tabulate import tabulate

class datum:
    def __init__(self, value, err):
        self.value = value
        self.err = err

inputs = np.genfromtxt("inputs.txt", str,delimiter="=")
if int(inputs[0][1]) == 1:
    # Retrieves from inputs file
    importFName = str(inputs[1][1])
    interval = int(inputs[2][1])
    zp = float(inputs[3][1])
    zpErr = float(inputs[4][1])
    x0 = float(inputs[5][1])
    y0 = float(inputs[6][1])
    sma = float(inputs[7][1])
    eps = float(inputs[8][1])
    pa = float(inputs[9][1])
else:
    # Gets FITS file name for import
    importFName = str(input("FITS File Name: "))
    
    # Retrieves through manual input
    interval = int(input("Measurement Interval (px): "))
    zpInput = str(input("Zero point with error (separate with comma): "))
    zp, zpErr = zpInput.split(",")
    zp = float(zp)
    zpErr = float(zpErr)
    ellipseParamsInput = str(input("Initial ellipse guess params\nFormat: x,y,semi-major axis,ellipticity,position angle\n"))
    x0, y0, sma, eps, pa = [float(i) for i in ellipseParamsInput.split(",")]
    
# Adds .fits extension if not included
importFNameBare = importFName.split(".")[0]
importFName =  importFNameBare + ".fits"
# Names for the output files
exportDataFName = importFNameBare + "_data.dat"
exportParamsFName = importFNameBare + "_params.dat"
exportBackgroundFName = importFNameBare + "_bg.dat"

# Open FITS file
fn = get_pkg_data_filename(importFName)
f = fits.open(fn)
data = fits.getdata(fn, ext=0)

# Extract pscale
header = fits.getheader(fn, ext=0)
pscale = header[16]
print("\n" + tabulate([["Plate Scale", pscale]], tablefmt="plain") + "\n")

# Get background data
sigmaClip = SigmaClip(sigma=3, maxiters=100)
threshold = detect_threshold(data, n_sigma=3, sigma_clip=sigmaClip)
segmentImg = detect_sources(data, threshold, n_pixels=4)
footprint = circular_footprint(radius=50)
mask = segmentImg.make_source_mask(footprint=footprint)
bkgMean, bkgMedian, bkgStd = sigma_clipped_stats(data, sigma=3.0, mask=mask)
print(tabulate([["Mean", bkgMean], ["Median", bkgMedian], ["Std Dev", bkgStd]], headers=["Background:",""], tablefmt="plain") + "\n")

# Removes background
for n in data:
    n[n < 0] = bkgMean
cleanData = data - bkgMean

# Fit isophote
sample = EllipseSample(image=cleanData, sma=sma, x0=x0, y0=y0, astep=0.1, eps=eps, position_angle=pa * np.pi / 180)
fitter = EllipseFitter(sample)
isophote = fitter.fit()
err = bkgStd

# Photometry function
def do_photometry(storage, a, data, err, isophote, eps, pa):
    for n,i in enumerate(a):
        aper = EllipticalAperture(positions=(isophote.x0, isophote.y0), a=i, b=i - eps * i, theta=pa)
        photometry = aper.do_photometry(data=data, error=np.full_like(data, err))
        if aper.bbox.ixmin < 0 or aper.bbox.iymin < 0 or aper.bbox.ixmax >= len(data[0]) - 1 or aper.bbox.iymax >= len(data) - 1:
            a = np.delete(a,np.arange(n,len(a),1))
            break
        storage.append(datum(photometry[0][0], photometry[1][0]))
    return storage

a = np.arange(5,5001,interval)
count = []
countErr = []
countEPos = []
countENeg = []
countPaPos = []
countPaNeg = []

for n,i in enumerate(a):
    aper = EllipticalAperture(positions=(isophote.x0, isophote.y0), a=i, b=i - isophote.eps * i, theta=isophote.pa)
    photometry = aper.do_photometry(data=cleanData, error=np.full_like(data, err))
    if aper.bbox.ixmin < 0 or aper.bbox.iymin < 0 or aper.bbox.ixmax >= len(data[0]) - 1 or aper.bbox.iymax >= len(data) - 1:
        a = np.delete(a,np.arange(n,len(a),1))
        break
    count.append(photometry[0][0])
    countErr.append(photometry[1][0])
countEPos = do_photometry(countEPos, a, cleanData, err, isophote, isophote.eps + isophote.ellip_err, isophote.pa)
countENeg = do_photometry(countENeg, a, cleanData, err, isophote, isophote.eps - isophote.ellip_err, isophote.pa)
countPaPos = do_photometry(countPaPos, a, cleanData, err, isophote, isophote.eps, isophote.pa + isophote.pa_err)
countPaNeg = do_photometry(countPaNeg, a, cleanData, err, isophote, isophote.eps, isophote.pa - isophote.pa_err)

ellipse_err = [np.sqrt(((countEPos[index].value - countENeg[index].value) / 2) ** 2 + ((countPaPos[index].value - countPaNeg[index].value) / 2) ** 2 + countErr[index] ** 2) for index in np.arange(np.min([len(countErr), len(countEPos), len(countENeg), len(countPaPos), len(countPaNeg)]))]

# Convert to relative magnitude
relMag = np.array([zp - 2.5 * np.log10(n) for n in count[:len(ellipse_err)]])
relMagErr = np.array([np.sqrt((2.5 * ellipse_err[n] / (count[n] * np.log(10))) ** 2 + zpErr ** 2) for n in np.arange(len(ellipse_err))])

# Print params
print(tabulate([["Centre X-coordinate", isophote.x0], ["Centre X-coordinate Error", isophote.x0_err], ["Centre Y-coordinate", isophote.y0], ["Centre Y-coordinate Error", isophote.y0_err], ["Ellipticity", isophote.eps], ["Ellipticity Error", isophote.ellip_err], ["Position Angle (Deg)", isophote.pa * 180 / np.pi], ["Position Angle Error (Deg)", isophote.pa_err * 180 / np.pi]], headers=["Ellipse Parameters:",""], tablefmt="plain") + "\n")

# Exporting the data
exportData = np.column_stack((a[:len(ellipse_err)], relMag, relMagErr))
np.savetxt(exportDataFName, exportData, delimiter=',')
np.savetxt(exportParamsFName,[[pscale, isophote.x0, isophote.x0_err, isophote.y0, isophote.y0_err, isophote.eps, isophote.ellip_err, isophote.pa * 180 / np.pi, isophote.pa_err * 180 / np.pi]], delimiter=",", header="Plate Scale, Centre X-coordinate, Centre X-coordinate Error, Centre Y-coordinate, Centre Y-coordinate Error, Ellipticity, Ellipticty Error, Position Angle, Position Angle Error")
np.savetxt(exportBackgroundFName,[[bkgMean, bkgMedian, bkgStd]], delimiter=",", header="Background Mean, Background Median, Background Standard Deviation")