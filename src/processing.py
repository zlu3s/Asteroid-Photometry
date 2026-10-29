import numpy as np
from astropy.io import fits

def reduce_image(img_path, bias_path=None, flat_path=None, dark_path=None):
    image = fits.getdata(img_path)
    bias = flat = dark = 0
    if bias_path: bias = fits.getdata(bias_path) 
    if flat_path: flat = fits.getdata(flat_path) 
    if dark_path: dark = fits.getdata(dark_path)

    flat_field_norm = (flat - bias) / np.median(flat - bias)
    reduced_image = image / flat_field_norm
    return reduced_image


def find_ref_star(image):

