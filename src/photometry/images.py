import numpy as np
from astropy.io import fits
import sys, os

home = os.environ.get("PIPELINE")

class MultiCalibrateImage:
    def __init__(self,paths):
        self.imgs = [fits.getdata(p) for p in paths]
        self.master = np.median(self.imgs, axis=0) 


class RawFitsImage:
    def __init__(self, img_path, bias, flat):
        self.img_path = img_path
        self.data = fits.getdata(self.img_path)
        self.header = fits.getheader(self.img_path)
        self.bias = bias
        self.flat = flat
        self.flat_field_norm = (self.flat - self.bias) / np.median(self.flat - self.bias)
        self.reduced_image = self.data / self.flat_field_norm

    def save_reduced_image(self, path):
        reduced_image = self.reduced_image
        fits.writeto(path, reduced_image, header=self.header, overwrite=True)
