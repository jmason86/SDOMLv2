import os
import numpy as np
from scipy.io import readsav
import zarr
from numcodecs import Blosc, Delta
from astropy.table import Table
from astropy.time import Time
import pandas as pd

delta_seconds_1958_to_unix = 378691200  # TAI's formal epoch starts at 1958-01-01 while astropy uses unix_tai, which is from 1970-01-01


def eve_fits_to_zarr():
    #lines_to_zarr() # "lines" here refers to the extracted emission lines data product
    spectra_to_zarr() # and this is the full spectra data product


def lines_to_zarr(lines_output_filename='eve_lines.zarr'):
    lines = read_lines_source_file()
    lines_zarr = store_lines_in_zarr(lines, lines_output_filename)


def read_lines_source_file(lines_source_file=os.path.expanduser('~') + '/Dropbox/Research/Data/SDO/EVE/eve_lines_2010121-2014146 MEGS-A Mission Bare Bones.sav'):
    return readsav(lines_source_file)


def store_lines_in_zarr(lines, lines_output_filename='eve.zarr'):
    store = zarr.DirectoryStore(lines_output_filename)
    compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
    root = zarr.group(store=store,overwrite=True)

    format_lines_for_zarr(lines, root)


def format_lines_for_zarr(lines, root):
    pass # TODO Populate based on Meng's code


def spectra_to_zarr():
    dataloc = '/Volumes/Roci Extension/sdo/eve/spectra/'
    filename = dataloc + 'EVS_L2_2011013_01_007_02.fit.gz'
    bla = read_spectra_source_file(filename)
    
    pass # TODO: Populate


def read_spectra_source_file(filename):
    meta = Table.read(filename, format='fits', hdu=2)
    data = Table.read(filename, format='fits', hdu=3)
    time = convert_time(data)
    pass


def convert_time(data):
    time = []
    for tai in data['TAI']:
        time.append(Time(tai - delta_seconds_1958_to_unix, format='unix_tai').iso)
    return time



if __name__ == "__main__":
    eve_fits_to_zarr()