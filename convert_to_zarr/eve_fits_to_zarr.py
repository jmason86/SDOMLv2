import os
import shutil
import glob
import numpy as np
from scipy.io import readsav
import zarr
from numcodecs import Blosc, Delta
from astropy.table import Table, vstack
from astropy.time import Time
import pandas as pd
from datetime import datetime

eve_path = '/Volumes/Roci Extension/sdo/eve/'
eve_output_filename = eve_path + 'eve_L2B_2020-2023.zarr'
delta_seconds_1958_to_unix = Time('1958-01-01').unix_tai  # TAI's formal epoch starts at 1958-01-01 while astropy uses unix_tai, which is from 1970-01-01


def eve_fits_to_zarr(do_remove_old_zarr_file=False):    
    meta_lines, data_lines, evl_filenames_to_process = read_lines()
    meta_spectra, data_spectra, evs_filenames_to_process = read_spectra()

    if do_remove_old_zarr_file: 
        remove_old_zarr_file(eve_output_filename)
        setup_zarr_file(eve_output_filename, meta_lines, data_lines, meta_spectra, data_spectra)
    store_lines_in_zarr(meta_lines, data_lines, eve_output_filename)
    store_spectra_in_zarr(meta_spectra, data_spectra, eve_output_filename)

    log_processed_files(evl_filenames_to_process)
    log_processed_files(evs_filenames_to_process)


def remove_old_zarr_file(output_filename): 
    if os.path.exists(output_filename):
        shutil.rmtree(output_filename)
    

def read_lines(): 
    all_evl_filenames = sorted(glob.glob(eve_path + 'lines/EVL_L2B_*.fit*'))
    start_filename = get_start_filename_from_log(evl_or_evs='evl')
    start_index = all_evl_filenames.index(start_filename)
    end_index = min(start_index + 10, len(all_evl_filenames)) # batch to 10 files at a time
    filenames_to_process = all_evl_filenames[start_index:end_index]

    meta_all = Table()
    data_all = Table()
    for filename in filenames_to_process:
        meta, data = read_lines_source_file(filename)
        meta_all = vstack([meta_all, meta])
        data_all = vstack([data_all, data])
    return meta_all, data_all, filenames_to_process


def get_start_filename_from_log(evl_or_evs='evl'): 
    df = pd.read_csv(eve_path + '/processing_log.txt', names=['timestamp', 'processing', 'filename'])
    if evl_or_evs == 'evl':
        df = df[df['filename'].str.contains('EVL_L2B_')]
        input_filenames = sorted(glob.glob(eve_path + 'lines/EVL_L2B_*.fit*'))  
    else: 
        df = df[df['filename'].str.contains('EVS_L2B_')]
        input_filenames = sorted(glob.glob(eve_path + 'spectra/EVS_L2B_*.fit*')) 
    
    last_filename = df['filename'].iloc[-1]
    filtered_filenames = [filename for filename in input_filenames if os.path.basename(filename) > os.path.basename(last_filename)]
    return filtered_filenames[0]


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


def read_spectra():
    all_evs_filenames = sorted(glob.glob(eve_path + 'spectra/EVS_L2B_*.fit*'))
    start_filename = get_start_filename_from_log(evl_or_evs='evs')
    start_index = all_evs_filenames.index(start_filename)
    end_index = min(start_index + 10, len(all_evs_filenames)) # batch to 10 files at a time
    filenames_to_process = all_evs_filenames[start_index:end_index]

    meta_all = Table()
    data_all = Table()
    for filename in filenames_to_process:
        meta, data = read_spectra_source_file(filename)
        meta_all = vstack([meta_all, meta])
        data_all = vstack([data_all, data])
    return meta_all, data_all, filenames_to_process


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


def setup_zarr_file(eve_output_filename, meta_lines, data_lines, meta_spectra, data_spectra):
    store = zarr.DirectoryStore(eve_output_filename)
    root = zarr.group(store=store, overwrite=False)

    lines = root.create_group('Lines')
    compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
    time = lines.create_dataset('TIME_ISO', shape=(np.shape(data_lines['TIME_ISO'])), dtype='<U23', compressor=compressor)
    wavelength = lines.create_dataset('WAVELENGTH', shape=(np.shape(meta_lines['WAVE_CENTER'])), dtype='>f4', compressor=compressor)
    logt = lines.create_dataset('LOG_T', shape=(np.shape(meta_lines['LOGT'])), dtype='>f4', compressor=compressor)
    ion = lines.create_dataset('ION', shape=(np.shape(meta_lines['NAME'])), dtype='<U23', compressor=compressor)
    blends = lines.create_dataset('BLENDS', shape=(np.shape(meta_lines['BLENDS'])), dtype='<U23', compressor=compressor)
    irradiance = lines.create_dataset('IRRADIANCE', shape=(np.shape(data_lines['IRRADIANCE'])), dtype='>f4', compressor=compressor)

    time.attrs['UNITS'] = 'Universal Coordinated Time (UTC) in ISO8601 format'
    wavelength[:] = meta_lines['WAVE_CENTER'].value.data
    wavelength.attrs['UNITS'] = 'nm'
    logt[:] = meta_lines['LOGT'].value.data
    logt.attrs['UNITS'] = 'log_10 of temperature in Kelvin'
    ion[:] = meta_lines['NAME'].value.data
    ion.attrs['NOTE'] = 'This is the _dominate_ ion in the spectral range integrated to create this line. There may be other emission lines blended in (see "blends").'
    blends[:] = meta_lines['BLENDS'].value.data
    blends.attrs['NOTE'] = 'These are the other emission lines that are blended with the dominate emission listed in "ion".'
    irradiance.attrs['UNITS'] = 'W m^-2'

    spectra = root.create_group('Spectra')
    compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
    time = spectra.create_dataset('TIME_ISO', shape=(np.shape(data_spectra)[0]), dtype='<U23', compressor=compressor)
    irradiance = spectra.create_dataset('IRRADIANCE', shape=(np.shape(data_spectra['IRRADIANCE'])), dtype='>f4', compressor=compressor)

    time.attrs['UNITS'] = meta_spectra['TIME_ISO'].value[0]
    irradiance.attrs['UNITS'] = str(meta_spectra['IRRADIANCE'].value[0])


def store_lines_in_zarr(meta, data, lines_output_filename):
    store = zarr.DirectoryStore(lines_output_filename)
    root = zarr.group(store=store,overwrite=False)

    lines = root['Lines']
    time = lines['TIME_ISO']
    wavelength = lines['WAVELENGTH']
    logt = lines['LOG_T']
    ion = lines['ION']
    blends = lines['BLENDS']
    irradiance = lines['IRRADIANCE']

    time.append(data['TIME_ISO'].value)
    #wavelength.append(meta['WAVE_CENTER'].value)
    #logt.append(meta['LOGT'].value)
    #ion.append(meta['NAME'].value)
    #blends.append(meta['BLENDS'].value)
    irradiance.append(data['IRRADIANCE'].value)


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
    
    spectra = root['Spectra']
    time = spectra['TIME_ISO']
    irradiance = spectra['IRRADIANCE']
    
    time.append(data['TIME_ISO'].value)
    irradiance.append(data['IRRADIANCE'].value)


def log_processed_files(filenames): 
    file = open(eve_path + '/processing_log.txt', 'a')
    for filename in filenames: 
        file.write('\n' + datetime.now().strftime('%Y-%m-%dT%H:%M:%S') + ',Processed file,' + filename)
    file.close()


if __name__ == "__main__":
    while True:
        eve_fits_to_zarr()