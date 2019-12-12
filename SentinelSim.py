"""
Sentinel-2 Emulator
A Simulator for estimating sentinel 2 responses to arbitrary spectra

Created by Joseph T. Fennell
https://github.com/joe-fennell/sentinel2_sim/

"""

import numpy as np
import os
import pandas
import warnings
import xarray

def simulate(spectrum, force_out_of_range=False):
    """
    Simulate surface reflectance for every S2A band

    PARAMETERS
    ----------
    spectrum : xarray.DataArray
        spectrum data array with labelled wavelength coordinate

    force_out_of_range : bool, optional
        flag to force return of a value for out-of-range measurement.
        Only use if there is a very small difference

    RETURNS
    -------
    values : xarray.DataArray
        Surface Level Reflectance values labelled by Sentinel 2 Band and
        peak wavelength (nm)

    """
    s = SimpleSimulator()

    return s.simulate(spectrum, force_out_of_range)

class SimpleSimulator(object):
    """
    Simulator class for simulating a reflectance spectrum shape with the ASTM G173
    model solar spectrum. Values are not normalised and
    """
    def __init__(self):
        # read spectra and persist
        self.spectra, self.solar_spectrum = _read_spectra()

    def simulate(self, spectrum, force_out_of_range=False):
        """
        Simulate surface reflectance for every S2A band

        PARAMETERS
        ----------
        spectrum : xarray.DataArray
            spectrum data array with labelled wavelength coordinate

        force_out_of_range : bool, optional
            flag to force return of a value for out-of-range measurement.
            Only use if there is a very small difference

        RETURNS
        -------
        values : xarray.DataArray
            Surface Level Reflectance values labelled by Sentinel 2 Band and
            peak wavelength (nm)

        """
        if 'wavelength' not in spectrum.coords:
            raise ValueError('Spectrum must have a \'wavelength\' dimension')

        signal = spectrum * self.solar_spectrum

        out = []
        bands = []
        peaks = []
        for band, r in self.spectra.items():

            # check if out of range for a given band
            if _is_out_of_range(signal,r):
                warnings.warn('Input spectrum is out of range of {}'.format(band),stacklevel=2,)
                if force_out_of_range:
                    out.append((r*signal).integrate('wavelength')/(self.solar_spectrum*r).integrate('wavelength'))
                    warnings.warn('Value for {} may be inaccurate'.format(band),stacklevel=2)
                else:
                    out.append(np.nan)
                    warnings.warn('Value for {} not calculated'.format(band),stacklevel=2)
            else:
                out.append((r*signal).integrate('wavelength')/(self.solar_spectrum*r).integrate('wavelength'))
            bands.append(band)
            peaks.append(r.wavelength[r.argmax()])

        return xarray.DataArray(out,dims='x',coords={'peak_wavelength':('x',peaks),'band':('x',bands)}).sortby('peak_wavelength')

def read_file(fpath,header=None):
    """
    Reads a csv file into an xarray object. The first column should be a
    vector of wavelengths and the second a vector of corresponding values.

    PARAMETERS
    ----------
    fpath : str
        path/to/file.csv

    header : int, optional
        number of header lines to ignore. Default is None

    RETURNS
    -------
    spectrum : xarray.DataArray
        1D DataArray with a labelled 'wavelength' dimension
    """
    # expects all files in csv format, first column wavelength, second column response
    df = pandas.read_csv(fpath,header=None)
    return xarray.DataArray(df.iloc[:,1].astype('float64'),dims=['wavelength'],coords={'wavelength':df.iloc[:,0].astype('float64')})

# private funcs
def _read_spectra(s2_path='spectra/S2A/',astm_path='spectra/ASTMG173_spectrum.csv'):
    # reads all spectra into a dataset
    s2_paths = [os.path.join(s2_path,x) for x in os.listdir(s2_path) if x.endswith('.csv')]
    spectra = {}
    for path in s2_paths:
        band_name = os.path.basename(path).split('.')[0].split('_')[-1]
        # read as xarray and interpolate
        xr = read_file(path).interp(wavelength=np.arange(412,2320,1),kwargs={'fill_value':0})
        xr.attrs['Range'] = _min_max(xr)
        spectra[band_name] = xr

    return spectra, read_file(astm_path).interp(wavelength=np.arange(412,2320,1),kwargs={'fill_value':0})

def _min_max(spectrum):
    # max of spectrum
    max_ = np.max(spectrum.wavelength[spectrum.values>0]).values
    # min of spectrum
    min_ = np.min(spectrum.wavelength[spectrum.values>0]).values
    return (float(min_), float(max_))

def _is_out_of_range(signal,sensor):
    sig_min = signal.wavelength.min()
    sig_max = signal.wavelength.max()
    if sig_min > sensor.attrs['Range'][0]:
        return True
    if sig_max < sensor.attrs['Range'][1]:
        return True
    return False
