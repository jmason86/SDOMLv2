import os
import numpy as np
from scipy.io import readsav
import zarr
from numcodecs import Blosc, Delta


def eve_fits_to_zarr():
    lines_to_zarr() # "lines" here refers to the extracted emission lines data product
    spectra_to_zarr() # and this is the full spectra data product


def lines_to_zarr(lines_output_filename='eve_lines.zarr'):
    lines = load_lines_source_file()
    lines_zarr = store_lines_in_zarr(lines, lines_output_filename)


def load_lines_source_file(lines_source_file=os.path.expanduser('~') + '/Dropbox/Research/Data/SDO/EVE/eve_lines_2010121-2014146 MEGS-A Mission Bare Bones.sav'):
    return readsav(lines_source_file)


def store_lines_in_zarr(lines, lines_output_filename='eve.zarr'):
    store = zarr.DirectoryStore(lines_output_filename)
    compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
    root = zarr.group(store=store,overwrite=True)

    format_lines_for_zarr(lines, root)


def format_lines_for_zarr(lines, root):
    pass # TODO Populate based on Meng's code


def spectra_to_zarr():
    pass # TODO: Populate




if __name__ == "__main__":
    eve_fits_to_zarr()