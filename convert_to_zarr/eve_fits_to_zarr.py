import os
import shutil
import glob
import numpy as np
from scipy.io import readsav
import zarr
from numcodecs import Blosc, Delta
from astropy.table import Table
from astropy.time import Time
import pandas as pd

delta_seconds_1958_to_unix = Time('1958-01-01').unix_tai  # TAI's formal epoch starts at 1958-01-01 while astropy uses unix_tai, which is from 1970-01-01


def eve_fits_to_zarr():
    eve_path = '/Volumes/Roci Extension/sdo/eve/'
    lines_to_zarr(path_to_input_data=eve_path + '/lines/', 
                  lines_output_filename=eve_path + '/lines/eve_lines.zarr') # "lines" here refers to the extracted emission lines data product
    spectra_to_zarr(path_to_input_data=eve_path + 'spectra/', 
                    spectra_output_filename=eve_path + '/spectra/eve_spectra.zarr') # and this is the full spectra data product


def lines_to_zarr(path_to_input_data, lines_output_filename='eve_lines.zarr'):
    remove_old_zarr_file(lines_output_filename)

    input_filenames = glob.glob(path_to_input_data + 'EVL_L2_*.fit*')
    for filename in input_filenames:
        meta, data = read_lines_source_file(filename) 
        store_lines_in_zarr(meta, data, lines_output_filename) 


def read_lines_source_file(filename):
    meta = Table.read(filename, format='fits', hdu=1)
    meta2 = Table.read(filename, format='fits', hdu=6)
    data = Table.read(filename, format='fits', hdu=5)
    
    time = convert_time(data)
    data.add_column(time, name='TIME_ISO', index=0)
    meta.add_column('ISO8601 UTC Time', name='TIME_ISO', index=0)
    meta.add_column(meta2['LINE_IRRADIANCE'], name='IRRADIANCE', index=2)
    data['IRRADIANCE'] = replace_fill_values(data)

    # Drop unnecessary columns (TAI, YYYYDOY, SOD, etc; everything in meta is good to keep)
    data = data['TIME_ISO', 'IRRADIANCE']  

    return meta, data


def store_lines_in_zarr(meta, data, lines_output_filename):
    store = zarr.DirectoryStore(lines_output_filename)
    root = zarr.group(store=store,overwrite=False)

    if 'Lines' not in root:
        lines = root.create_group('Lines')
        compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
        time = lines.create_dataset('TIME_ISO', shape=(np.shape(data['TIME_ISO'])), dtype='<U23', compressor=compressor)
        wavelength = lines.create_dataset('WAVELENGTH', shape=(np.shape(meta['WAVE_CENTER'])), dtype='>f4', compressor=compressor)
        logt = lines.create_dataset('LOG_T', shape=(np.shape(meta['LOGT'])), dtype='>f4', compressor=compressor)
        ion = lines.create_dataset('ION', shape=(np.shape(meta['NAME'])), dtype='<U23', compressor=compressor)
        blends = lines.create_dataset('BLENDS', shape=(np.shape(meta['BLENDS'])), dtype='<U23', compressor=compressor)
        irradiance = lines.create_dataset('IRRADIANCE', shape=(np.shape(data['IRRADIANCE'])), dtype='>f4', compressor=compressor)

        time.attrs['UNITS'] = 'Universal Coordinated Time (UTC) in ISO8601 format'
        wavelength.attrs['UNITS'] = 'nm'
        logt.attrs['UNITS'] = 'log_10 of temperature in Kelvin'
        ion.attrs['NOTE'] = 'This is the _dominate_ ion in the spectral range integrated to create this line. There may be other emission lines blended in (see "blends").'
        blends.attrs['NOTE'] = 'These are the other emission lines that are blended with the dominate emission listed in "ion".'
        irradiance.attrs['UNITS'] = 'W m^-2'
    else: 
        lines = root['Lines']
        time = lines['TIME_ISO']
        wavelength = lines['WAVELENGTH']
        logt = lines['LOG_T']
        ion = lines['ION']
        blends = lines['BLENDS']
        irradiance = lines['IRRADIANCE']

    time.append(data['TIME_ISO'].value)
    wavelength.append(meta['WAVE_CENTER'].value)
    logt.append(meta['LOGT'].value)
    ion.append(meta['NAME'].value)
    blends.append(meta['BLENDS'].value)
    irradiance.append(data['IRRADIANCE'].value)


def spectra_to_zarr(path_to_input_data, spectra_output_filename='eve_spectra.zarr'):
    remove_old_zarr_file(spectra_output_filename)
    
    input_filenames = glob.glob(path_to_input_data + 'EVS_L2_*.fit*')
    for filename in input_filenames:
        meta, data = read_spectra_source_file(filename) 
        store_spectra_in_zarr(meta, data, spectra_output_filename) 


def remove_old_zarr_file(output_filename): 
    if os.path.exists(output_filename):
        shutil.rmtree(output_filename)


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
    try:
        irradiance = data['IRRADIANCE']
    except: 
        irradiance = data['LINE_IRRADIANCE']
    irradiance[irradiance == -1] = np.nan  # Required by Zarr spec https://zarr.readthedocs.io/en/stable/spec/v2.html#fill-value-encoding
    return irradiance


def store_spectra_in_zarr(meta, data, spectra_output_filename):
    store = zarr.DirectoryStore(spectra_output_filename)
    root = zarr.group(store=store, overwrite=False)
    
    if 'Spectra' not in root: 
        spectra = root.create_group('Spectra')
        compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
        time = spectra.create_dataset('TIME_ISO', shape=(np.shape(data)[0]), dtype='<U23', compressor=compressor)
        irradiance = spectra.create_dataset('IRRADIANCE', shape=(np.shape(data['IRRADIANCE'])), dtype='>f4', compressor=compressor)

        time.attrs['UNITS'] = meta['TIME_ISO'].value[0]
        irradiance.attrs['UNITS'] = str(meta['IRRADIANCE'].value[0])
    else: 
        spectra = root['Spectra']
        time = spectra['TIME_ISO']
        irradiance = spectra['IRRADIANCE']
    
    time.append(data['TIME_ISO'].value)
    irradiance.append(data['IRRADIANCE'].value)


if __name__ == "__main__":
    eve_fits_to_zarr()