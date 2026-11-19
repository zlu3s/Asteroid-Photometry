#!/usr/bin/env python3
"""select_sources.py - click on a light source in each image and collect positions.

Workflow
--------
1. Give the script one or more images (FITS, or anything Pillow can open).
2. Every detected light source is circled in green.
3. Left-click near the source you want; it turns red and its position is stored.
4. Move between images, then quit. Positions are returned as a list (one entry
   per image, None where nothing was selected) and written to a CSV file.

Controls
--------
  left click   select the source closest to the click (click again to change)
  right click  next image
  d / a        next / previous image
  u            clear the selection for the current image
  + / -        zoom in / out
  q            quit and save (closing the window does the same)

Coordinates are reported in image pixels using the 1-based convention of
FITS / SExtractor (the centre of the first pixel is (1, 1)).

Usage
-----
  python select_sources.py img1.fits img2.fits img3.fits
  python select_sources.py *.fits --nsigma 5 --min-pixels 4 -o positions.csv

Or from Python:

  from select_sources import select_sources
  positions = select_sources(['a.fits', 'b.fits'])
  # -> [(x, y), (x, y)]  or None entries for images with no selection
"""
import argparse
import csv
import os
import sys
import tkinter as tk

import numpy as np
from PIL import Image, ImageTk
from scipy import ndimage

try:
    from astropy.io import fits
except ImportError:  # FITS support is optional if only using PNG/JPG etc.
    fits = None


# ----------------------------------------------------------------------
# image loading and source detection
# ----------------------------------------------------------------------

def load_image(filename):
    """Return the image as a 2D float array (row 0 is displayed at the top)."""
    if filename.lower().endswith(('.fit', '.fits', '.fts', '.fit.gz', '.fits.gz')):
        if fits is None:
            raise RuntimeError('astropy is required to read FITS files')
        with fits.open(filename, ignore_missing_end=True) as hdul:
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim >= 2:
                    data = hdu.data
                    break
            else:
                raise ValueError('no image data found in %s' % filename)
        while data.ndim > 2:      # take the first plane of cubes
            data = data[0]
    else:
        data = np.asarray(Image.open(filename).convert('L'))
    return np.nan_to_num(data.astype(float))


def robust_stats(data):
    """Median and robust sigma (from the MAD) of an array."""
    med = np.median(data)
    sigma = 1.4826 * np.median(np.abs(data - med))
    if sigma <= 0:
        sigma = np.std(data) or 1.0
    return med, sigma


def detect_sources(data, nsigma=5.0, min_pixels=4, max_sources=2000):
    """Find bright blobs and return an array of (x, y, flux), 0-based, brightest first."""
    smooth = ndimage.gaussian_filter(data, 1.0)
    bkg, sigma = robust_stats(smooth)

    mask = smooth > bkg + nsigma * sigma
    labels, n = ndimage.label(mask)
    if n == 0:
        return np.empty((0, 3))

    index = np.arange(1, n + 1)
    sizes = ndimage.sum(mask, labels, index)
    index = index[sizes >= min_pixels]
    if len(index) == 0:
        return np.empty((0, 3))

    weights = np.clip(data - bkg, 0, None)
    centres = np.array(ndimage.center_of_mass(weights, labels, index))  # (y, x)
    flux = np.asarray(ndimage.sum(weights, labels, index))

    order = np.argsort(flux)[::-1][:max_sources]
    return np.column_stack([centres[order, 1], centres[order, 0], flux[order]])


def to_display(data):
    """Scale a float array to an 8-bit PIL image for display."""
    med, sigma = robust_stats(data)
    lo, hi = med - sigma, med + 5 * sigma
    scaled = (np.clip(data, lo, hi) - lo) / (hi - lo) * 255
    return Image.fromarray(scaled.astype(np.uint8))


# ----------------------------------------------------------------------
# GUI
# ----------------------------------------------------------------------

class SourcePicker(object):

    R = 6  # circle radius in screen pixels

    def __init__(self, master, filenames, zoom=None, nsigma=5.0, min_pixels=4):
        self.top = master
        self.files = filenames
        self.index = 0
        self.selected = [None] * len(filenames)   # index into sources, per image
        self.images, self.sources = [], []

        print('loading images and detecting sources...')
        for i, fn in enumerate(filenames):
            data = load_image(fn)
            src = detect_sources(data, nsigma, min_pixels)
            self.images.append(to_display(data))
            self.sources.append(src)
            print('  %s: %d sources' % (os.path.basename(fn), len(src)))

        if zoom is None:   # fit the first image to ~90% of the screen
            w, h = self.images[0].size
            sw, sh = master.winfo_screenwidth(), master.winfo_screenheight()
            zoom = min(1.0, 0.9 * sw / w, 0.85 * sh / h)
        self.zoom = zoom

        self.title = tk.Label(master)
        self.title.pack()
        self.canvas = tk.Canvas(master)
        self.canvas.pack(side='top', expand=1)
        self.canvas.focus_set()
        self.canvas.bind('<Key>', self.on_key)
        self.canvas.bind('<Button-1>', self.on_left_click)
        self.canvas.bind('<Button-3>', lambda e: self.show(self.index + 1))
        self.canvas.bind('<Button-2>', lambda e: self.show(self.index + 1))  # macOS

        self.show(0)

    # -- drawing -------------------------------------------------------

    def show(self, index):
        self.index = index % len(self.files)
        img = self.images[self.index]
        w, h = int(img.width * self.zoom), int(img.height * self.zoom)
        self.tkimage = ImageTk.PhotoImage(img.resize((w, h), Image.BILINEAR))

        self.canvas.config(width=w, height=h)
        self.canvas.delete('all')
        self.canvas.create_image(0, 0, anchor='nw', image=self.tkimage)

        src = self.sources[self.index]
        sel = self.selected[self.index]
        for i, (x, y, _) in enumerate(src):
            if i == sel:
                continue
            self._circle(x, y, 'green', 1)
        if sel is not None:
            self._circle(src[sel, 0], src[sel, 1], 'red', 3)

        done = sum(s is not None for s in self.selected)
        self.title.config(text='%s  (%d/%d)   selected: %d/%d' % (
            os.path.basename(self.files[self.index]), self.index + 1,
            len(self.files), done, len(self.files)))

    def _circle(self, x, y, colour, width):
        cx, cy = (x + 0.5) * self.zoom, (y + 0.5) * self.zoom
        self.canvas.create_oval(cx - self.R, cy - self.R, cx + self.R,
                                cy + self.R, outline=colour, width=width)

    # -- events --------------------------------------------------------

    def on_left_click(self, event):
        src = self.sources[self.index]
        if len(src) == 0:
            return
        x, y = event.x / self.zoom - 0.5, event.y / self.zoom - 0.5
        self.selected[self.index] = int(np.argmin(np.hypot(src[:, 0] - x, src[:, 1] - y)))
        self.show(self.index)

    def on_key(self, event):
        c = event.char
        if c == 'd':
            self.show(self.index + 1)
        elif c == 'a':
            self.show(self.index - 1)
        elif c == 'u':
            self.selected[self.index] = None
            self.show(self.index)
        elif c == '+' or c == '=':
            self.zoom = min(self.zoom * 2, 8.0)
            self.show(self.index)
        elif c == '-':
            self.zoom = max(self.zoom / 2, 0.05)
            self.show(self.index)
        elif c == 'q':
            self.top.quit()

    # -- results -------------------------------------------------------

    def positions(self):
        """List of (x, y) per image (1-based pixels), None if not selected."""
        out = []
        for src, sel in zip(self.sources, self.selected):
            out.append(None if sel is None
                       else (float(src[sel, 0]) + 1, float(src[sel, 1]) + 1))
        return out


def select_sources(filenames, zoom=None, nsigma=5.0, min_pixels=4):
    """Open the picker and return the list of selected (x, y) positions."""
    root = tk.Tk()
    app = SourcePicker(root, filenames, zoom, nsigma, min_pixels)
    root.mainloop()
    result = app.positions()
    root.destroy()
    return result


def save_positions(filenames, positions, outfile):
    with open(outfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['filename', 'x', 'y'])
        for fn, pos in zip(filenames, positions):
            writer.writerow([fn] + (['', ''] if pos is None else ['%.3f' % pos[0], '%.3f' % pos[1]]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='click on a light source in each image')
    parser.add_argument('images', nargs='+', help='images to process')
    parser.add_argument('--zoom', type=float, default=None,
                        help='display zoom factor (default: fit to screen)')
    parser.add_argument('--nsigma', type=float, default=5.0,
                        help='detection threshold above background in sigma (default 5)')
    parser.add_argument('--min-pixels', type=int, default=4,
                        help='minimum connected pixels for a source (default 4)')
    parser.add_argument('-o', '--output', default='positions.csv',
                        help='output CSV file (default positions.csv)')
    args = parser.parse_args()

    positions = select_sources(args.images, args.zoom, args.nsigma, args.min_pixels)
    save_positions(args.images, positions, args.output)

    for fn, pos in zip(args.images, positions):
        print('%-40s %s' % (os.path.basename(fn),
                            'not selected' if pos is None else '%.3f  %.3f' % pos))
    print('saved to', args.output)
