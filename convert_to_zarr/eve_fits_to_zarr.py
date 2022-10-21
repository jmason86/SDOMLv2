import os
import numpy as np
from scipy.io import readsav
import zarr
from numcodecs import Blosc, Delta
from astropy.table import Table
from astropy.time import Time
import pandas as pd

delta_seconds_1958_to_unix = Time('1958-01-01').unix_tai  # TAI's formal epoch starts at 1958-01-01 while astropy uses unix_tai, which is from 1970-01-01


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
    meta, data = read_spectra_source_file(filename) 
    tmp = store_spectra_in_zarr(meta, data)

    pass # TODO: Populate


def read_spectra_source_file(filename):
    meta = Table.read(filename, format='fits', hdu=2)
    data = Table.read(filename, format='fits', hdu=3)
    
    time = convert_time(data)
    data.add_column(time, name='TIME_ISO', index=0)
    meta.add_column('ISO8601 UTC Time', name='TIME_ISO', index=0)
    data['IRRADIANCE'] = replace_fill_values(data)

    # Drop unnecessary columns (TAI, YYYYDOY, SOD, FLAGS, SC_FLAGS, INT_TIME, COUNT_RATE, PRECISION, BIN_FLAGS)
    data = data['TIME_ISO', 'IRRADIANCE']  
    meta = meta['TIME_ISO', 'IRRADIANCE']

    return meta, data  # TODO: Could have a test here to make sure that these are returned in the expected order (meta is a table with length 1, for example)


def convert_time(data):
    time = []
    for tai in data['TAI']:
        time.append(Time(tai + delta_seconds_1958_to_unix, format='unix_tai').iso)
    return time


def replace_fill_values(data):
    irradiance = data['IRRADIANCE']
    irradiance[irradiance == -1] = np.nan  # Required by Zarr spec https://zarr.readthedocs.io/en/stable/spec/v2.html#fill-value-encoding
    return irradiance


def store_spectra_in_zarr(meta, data, spectra_output_filename='eve_spectra.zarr'):
    store = zarr.DirectoryStore(spectra_output_filename)
    compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
    
    root = zarr.group(store=store, overwrite=True)
    spectra = root.create_group('Spectra')
    
    time = spectra.create_dataset('TIME_ISO', shape=(np.shape(data)[0]), dtype='<U23', compressor=compressor)
    time[:] = data['TIME_ISO'].value
    time.attrs['UNITS'] = meta['TIME_ISO'].value[0]
    
    irradiance = spectra.create_dataset('IRRADIANCE', shape=(np.shape(data['IRRADIANCE'])), dtype='>f4', compressor=compressor)
    irradiance[:] = data['IRRADIANCE'].value
    irradiance.attrs['UNITS'] = str(meta['IRRADIANCE'].value[0])


    pass


if __name__ == "__main__":
    eve_fits_to_zarr()