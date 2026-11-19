import argparse
import glob
from images import MultiCalibrateImage, RawFitsImage
import sys, os

home = os.getenv('PIPELINE')

def gen_images(img_paths, bias_paths, flat_paths):
    images = []
    flat = MultiCalibrateImage(flat_paths)
    bias = MultiCalibrateImage(bias_paths)
    for path in img_paths:
        images.append(RawFitsImage(path, bias.master, flat.master))
    return images


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Total controller")
    parser.add_argument("filter", type=str, help="Images to process")
    parser.add_argument('-s', '--save', action='store_true', help="Save processed images")

    args = parser.parse_args()
    
    bias_files = glob.glob(f"{home}/dat/bias/r.ZERO*")
    flat_files = glob.glob(f"{home}/dat/flats/r.{args.filter}flat*")
    raw_images = glob.glob(f"{home}/dat/raw_images/niobe{args.filter}*")

    images = gen_images(raw_images, bias_files, flat_files)
    if args.save:
        for img in images:
            name = os.path.basename(img.img_path)
            path = f"{home}/dat/reduced_images/{name}"
            img.save_reduced_image(path)
