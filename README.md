# Sentinel Simulator
Module to provide estimates of the Bottom-of-Atmosphere (surface-level reflectance) using the ASTM G173 simulated spectrum and the sensor sensitivity of the S2A sensors.
Handles out-of-range values and returns in popular xarray.DataArray format.

## Installation
Not pip installable. Run scripts in root directory of this repository

## Requirements
- NumPy
- xarray
- pandas

## Example Useage

```python
from SentinelSim import read_file, simulate
# read the spectrum you want to simulate
my_spectrum = read_file('path/to/file')
# simulate
simulation = simulate(my_spectrum)
# plot a pseudo spectrum of band responses
simulation.plot(x='peak_wavelength')
```
