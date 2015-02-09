# Signal to Noise Calculator for CCD Photometry
Allows you to calculate the signal-to-noise ratio for a star observed with a telescope and CCD over a range of
magnitudes with a defined magnitude step.

This is an Flask app. The apps main code is `signal_to_noise.py`. The templates required as stored in the templates
directory and the css file is stored in the static directory.

All the functions required to perform the calculations can be found in `functions.py`. This module can be imported
and used with out creating the Flask app. E.g.

```python
import functions
results = functions.signal_to_noise(filter_name='V', mag_start=8.0, mag_end=9.0, dmag=0.5,
                                    tel_diam=10, qe=0.5, readnoise=10, pixsize=3, skymag=17,
                                    airmass=1.2, exptime=20.0, fwhm=6, aper_rad=8, gain=1.0,
                                    obs_diam=5.0, sat_level=60000)
for r in results:
    print r

Magnitude=8.0  star_adu=256108.292  sky_adu=16255.446  read_adu=2234.021  signal_to_noise=488.737
    saturated=False
Magnitude=8.5  star_adu=161593.408  sky_adu=16255.446  read_adu=2234.021  signal_to_noise=380.792
    saturated=False
Magnitude=9.0  star_adu=101958.548  sky_adu=16255.446  read_adu=2234.021  signal_to_noise=293.781
    saturated=False
```
