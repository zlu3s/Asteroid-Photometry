from find_stars import find_stars, identify_stars
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy import units as u
import os, sys

home = os.getenv("PIPELINE")
file = f"{home}/dat/reduced_images/niobeB.0110.fits"
table = find_stars(file)
ra = float(star1['ra'])
dec = float(star1['dec'])
target = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

job = Gaia.cone_search_async(target, radius=.5*u.arcmin)
results = job.get_results()
print(results)
