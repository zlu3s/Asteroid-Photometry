#!/usr/bin/env python3
"""
find_stars.py - Detect stars in a FITS image, identify them (Gaia + SIMBAD),
and compute their absolute magnitudes.

Usage:
    python find_stars.py image.fits
    python find_stars.py image.fits --csv stars.csv
    python find_stars.py image.fits --fwhm 4 --sigma 5 --maglim 17 --plot
    python find_stars.py image.fits --no-crossmatch      # detection only, works offline

Requirements:
    pip install numpy astropy photutils astroquery pyvo matplotlib
    (matplotlib only for --plot; astroquery/pyvo only for the cross-match.
     The cross-match needs an internet connection and a FITS file with a WCS.)

As a module:
    from find_stars import find_stars, identify_stars
    table = find_stars("image.fits")
    table = identify_stars("image.fits", table)
    print(table)

Absolute magnitude (Gaia G band):
    M_G = G + 5*log10(parallax_mas) - 10
    This does NOT correct for interstellar extinction, and it is only
    trustworthy when the parallax is measured well (see --min-plx-snr).
"""

import argparse
import sys
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder

warnings.simplefilter("ignore", category=FITSFixedWarning)

SIMBAD_TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap"


# ----------------------------------------------------------------------------
# Image loading and star detection
# ----------------------------------------------------------------------------
def load_image(path):
    """Return (2D data array, header) from the first HDU holding a 2D image."""
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = np.asarray(hdu.data, dtype=float)
                while data.ndim > 2:  # collapse cubes / RGB to first plane
                    data = data[0]
                return data, hdu.header
    raise ValueError(f"No 2D image data found in {path}")


def find_stars(path, fwhm=3.0, sigma=5.0, box_size=64, max_stars=None):
    """
    Detect stars in a FITS image.

    Returns an astropy Table with: id, x, y, ra, dec (if WCS), peak, flux,
    sharpness, roundness. Sorted brightest (by flux) first.
    """
    data, header = load_image(path)
    data = np.nan_to_num(data, nan=np.nanmedian(data), posinf=0.0, neginf=0.0)

    # Subtract an (possibly uneven) sky background
    try:
        bkg = Background2D(
            data,
            box_size=box_size,
            filter_size=(3, 3),
            sigma_clip=SigmaClip(sigma=3.0),
            bkg_estimator=MedianBackground(),
        )
        data_sub = data - bkg.background
        noise = bkg.background_rms_median
    except Exception:
        _, median, std = sigma_clipped_stats(data, sigma=3.0)
        data_sub = data - median
        noise = std

    finder = DAOStarFinder(fwhm=fwhm, threshold=sigma * noise)
    sources = finder(data_sub)

    if sources is None or len(sources) == 0:
        return Table(names=["id", "x", "y", "peak", "flux"], dtype=[int, float, float, float, float])

    sources.sort("flux", reverse=True)
    if max_stars:
        sources = sources[:max_stars]

    out = Table()
    out["id"] = np.arange(1, len(sources) + 1)
    out["x"] = np.round(sources["xcentroid"], 2)
    out["y"] = np.round(sources["ycentroid"], 2)

    wcs = WCS(header)
    if wcs.has_celestial:
        sky = wcs.pixel_to_world(sources["xcentroid"], sources["ycentroid"])
        out["ra"] = np.round(sky.ra.deg, 6)
        out["dec"] = np.round(sky.dec.deg, 6)

    out["peak"] = np.round(sources["peak"], 2)
    out["flux"] = np.round(sources["flux"], 2)
    out["sharpness"] = np.round(sources["sharpness"], 3)
    out["roundness"] = np.round(sources["roundness1"], 3)
    return out


# ----------------------------------------------------------------------------
# Cross-match: Gaia (magnitude, parallax) + SIMBAD (names)
# ----------------------------------------------------------------------------
def _field_cone(path):
    """Return (center SkyCoord, radius in degrees) covering the whole image."""
    data, header = load_image(path)
    ny, nx = data.shape
    wcs = WCS(header)
    center = wcs.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    corners = wcs.pixel_to_world([0, nx - 1, 0, nx - 1], [0, 0, ny - 1, ny - 1])
    radius = center.separation(corners).deg.max() * 1.05
    return center, float(1.5)


def _query_gaia(center, radius_deg, maglim):
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = -1
    adql = f"""
        SELECT source_id, ra, dec, parallax, parallax_error,
               phot_g_mean_mag, bp_rp
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),
                           CIRCLE('ICRS', {center.ra.deg:.6f}, {center.dec.deg:.6f}, {radius_deg:.5f}))
          AND phot_g_mean_mag < {maglim}
    """
    return Gaia.launch_job_async(adql).get_results()


def _simbad_names(gaia_ids):
    """Map Gaia DR3 source_id (as str) -> SIMBAD main identifier, where one exists."""
    import pyvo

    names = {}
    if not gaia_ids:
        return names
    tap = pyvo.dal.TAPService(SIMBAD_TAP)
    ids = list(gaia_ids)
    for i in range(0, len(ids), 200):  # chunk to keep queries small
        chunk = ids[i:i + 200]
        in_list = ",".join(f"'Gaia DR3 {g}'" for g in chunk)
        query = f"""
            SELECT basic.main_id AS main_id, ident.id AS gaia_ident
            FROM ident JOIN basic ON ident.oidref = basic.oid
            WHERE ident.id IN ({in_list})
        """
        try:
            res = tap.run_sync(query).to_table()
        except Exception as e:
            print(f"Warning: SIMBAD query failed for a chunk ({e})", file=sys.stderr)
            continue
        for row in res:
            gid = str(row["gaia_ident"]).replace("Gaia DR3", "").strip()
            names[gid] = " ".join(str(row["main_id"]).split())
    return names


def identify_stars(path, table, match_radius=2.0, maglim=19.0, min_plx_snr=5.0):
    """
    Add Gaia/SIMBAD identification and absolute magnitude to a detection table.

    Adds: gaia_id, name, G, bp_rp, parallax, plx_err, dist_pc, abs_mag_G,
          match_sep_arcsec, reliable
    """
    if "ra" not in table.colnames:
        raise ValueError(
            "This FITS file has no usable WCS (sky coordinates), so stars can't be "
            "cross-matched. Plate-solve it first (e.g. astrometry.net or ASTAP)."
        )
    if len(table) == 0:
        return table

    center, radius = _field_cone(path)
    if radius > 2.0:
        raise ValueError(
            f"Field radius is {radius:.1f} deg - too large for a single Gaia query. "
            "Crop the image or lower the field of view."
        )

    print(f"Querying Gaia DR3 (G < {maglim}) in a {radius * 60:.1f} arcmin radius cone...")
    gaia = _query_gaia(center, radius, maglim)
    print(f"  {len(gaia)} Gaia sources in the field")

    n = len(table)
    gaia_id = np.full(n, "", dtype=object)
    G = np.full(n, np.nan)
    bp_rp = np.full(n, np.nan)
    plx = np.full(n, np.nan)
    plx_err = np.full(n, np.nan)
    sep_as = np.full(n, np.nan)

    if len(gaia) > 0:
        det = SkyCoord(table["ra"] * u.deg, table["dec"] * u.deg)
        cat = SkyCoord(np.asarray(gaia["ra"]) * u.deg, np.asarray(gaia["dec"]) * u.deg)
        idx, sep, _ = det.match_to_catalog_sky(cat)
        ok = sep.arcsec < match_radius
        for i in np.where(ok)[0]:
            g = gaia[int(idx[i])]
            gaia_id[i] = str(int(g["source_id"]))
            G[i] = float(g["phot_g_mean_mag"])
            bp_rp[i] = float(g["bp_rp"]) if np.isfinite(g["bp_rp"]) else np.nan
            plx[i] = float(g["parallax"]) if np.isfinite(g["parallax"]) else np.nan
            plx_err[i] = float(g["parallax_error"]) if np.isfinite(g["parallax_error"]) else np.nan
            sep_as[i] = sep.arcsec[i]

    # Names from SIMBAD where the star is catalogued; otherwise fall back to Gaia ID
    matched = sorted({g for g in gaia_id if g})
    print(f"Looking up {len(matched)} matched stars in SIMBAD...")
    simbad = _simbad_names(matched)
    name = np.array(
        [simbad.get(g, f"Gaia DR3 {g}") if g else "" for g in gaia_id], dtype=object
    )

    # Distance and absolute magnitude (only meaningful for positive, well-measured parallax)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = plx / plx_err
        reliable = np.isfinite(snr) & (plx > 0) & (snr >= min_plx_snr)
        dist = np.where(plx > 0, 1000.0 / plx, np.nan)
        abs_mag = np.where(plx > 0, G + 5 * np.log10(plx) - 10, np.nan)

    table["gaia_id"] = gaia_id
    table["name"] = name
    table["G"] = np.round(G, 3)
    table["bp_rp"] = np.round(bp_rp, 3)
    table["parallax"] = np.round(plx, 4)
    table["plx_err"] = np.round(plx_err, 4)
    table["dist_pc"] = np.round(dist, 1)
    table["abs_mag_G"] = np.round(abs_mag, 3)
    table["match_sep_arcsec"] = np.round(sep_as, 3)
    table["reliable"] = reliable
    return table


# ----------------------------------------------------------------------------
# Plotting and CLI
# ----------------------------------------------------------------------------
def plot_result(path, table):
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval

    data, _ = load_image(path)
    vmin, vmax = ZScaleInterval().get_limits(data)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
    ax.scatter(table["x"], table["y"], s=60, facecolors="none", edgecolors="lime", linewidths=0.8)
    ax.set_title(f"{len(table)} stars detected")
    plt.tight_layout()
    plt.show()


def main():
    p = argparse.ArgumentParser(
        description="Find stars in a FITS image, get their names (Gaia/SIMBAD) and absolute magnitudes."
    )
    p.add_argument("fits_file", help="Path to the FITS image")
    p.add_argument("--fwhm", type=float, default=3.0, help="Star FWHM in pixels (default 3.0)")
    p.add_argument("--sigma", type=float, default=5.0, help="Detection threshold in sigma (default 5.0)")
    p.add_argument("--box", type=int, default=64, help="Background mesh size in pixels (default 64)")
    p.add_argument("--max", type=int, default=None, help="Keep only the N brightest detections")
    p.add_argument("--no-crossmatch", action="store_true", help="Skip Gaia/SIMBAD lookup (offline mode)")
    p.add_argument("--match-radius", type=float, default=2.0, help="Cross-match radius in arcsec (default 2.0)")
    p.add_argument("--maglim", type=float, default=19.0, help="Gaia G magnitude limit for the query (default 19)")
    p.add_argument("--min-plx-snr", type=float, default=5.0,
                   help="Min parallax/error to flag a distance as reliable (default 5)")
    p.add_argument("--csv", help="Also write results to this CSV file")
    p.add_argument("--plot", action="store_true", help="Show the image with detected stars circled")
    args = p.parse_args()

    try:
        table = find_stars(args.fits_file, args.fwhm, args.sigma, args.box, args.max)
        print(f"Detected {len(table)} stars in {args.fits_file}")
        if not args.no_crossmatch:
            table = identify_stars(
                args.fits_file, table, args.match_radius, args.maglim, args.min_plx_snr
            )
    except Exception as e:
        sys.exit(f"Error: {e}")

    print()
    table.pprint(max_lines=-1, max_width=-1)

    if "reliable" in table.colnames:
        n_match = int(np.sum(table["gaia_id"] != ""))
        n_rel = int(np.sum(table["reliable"]))
        print(
            f"\n{n_match}/{len(table)} detections matched to Gaia; "
            f"{n_rel} have a reliable absolute magnitude."
        )
        print("Note: abs_mag_G is not corrected for interstellar extinction.")

    if args.csv:
        table.write(args.csv, format="csv", overwrite=True)
        print(f"Saved to {args.csv}")
    if args.plot:
        plot_result(args.fits_file, table)


if __name__ == "__main__":
    main()
