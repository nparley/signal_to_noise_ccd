"""
This module contains all the functions required to do the signal to noise calculations. Allows you
to calculate the signal-to-noise ratio for a star observed with a telescope and CCD. You can vary
the parameters and find out what sort of results to expect. Code adapted from original perl code
here: http://stupendous.rit.edu/richmond/signal.shtml

This module can be imported and used output of Flask for example:

> import functions
> results = functions.signal_to_noise(filter_name='V', mag_start=8.0, mag_end=9.0, dmag=0.5,
                                      tel_diam=10, qe=0.5, readnoise=10, pixsize=3, skymag=17,
                                      airmass=1.2, exptime=20.0, fwhm=6, aper_rad=8, gain=1.0,
                                      obs_diam=5.0, sat_level=60000)
> for r in results:
...    print r
...
Magnitude=8.0  star_adu=256108.292  sky_adu=16255.446  read_adu=2234.021  signal_to_noise=488.737
    saturated=False
Magnitude=8.5  star_adu=161593.408  sky_adu=16255.446  read_adu=2234.021  signal_to_noise=380.792
    saturated=False
Magnitude=9.0  star_adu=101958.548  sky_adu=16255.446  read_adu=2234.021  signal_to_noise=293.781
    saturated=False

"""
from __future__ import print_function
from __future__ import division
import collections

__author__ = 'Neil Parley'


class Result(collections.namedtuple('Result', 'magnitude star_adu sky_adu '
                                              'read_adu signal_to_noise saturated')):
    """
    Named tuple for the results of the calculation. The tuple as the parameters:
    Magnitude, star level (ADU), sky level (ADU), read level (ADU), signal to noise and boolean
    value for if the star will saturate or not.
    """
    def __str__(self):
            """
            Overload the printing of the named tuple to format the output as we like
            """
            return 'Magnitude={mag}  star_adu={star_e:.3f}  sky_adu={sky_e:.3f}  ' \
                   'read_adu={read_e:.3f}  signal_to_noise={s_to_n:.3f}  saturated={sat}'\
                .format(mag=self.magnitude, star_e=self.star_adu, sky_e=self.sky_adu,
                        read_e=self.read_adu, s_to_n=self.signal_to_noise, sat=self.saturated)


def get_extinct_coeff(filter_name):
    """
    Given a filter name, return an extinction coefficient. Numbers are approximate.
    V-band value is used if the filter name is 'none'

    :type filter_name: str
    :param filter_name: Filter name string
    :return: Extinction coefficient

    >>> get_extinct_coeff("none")
    0.2
    >>> get_extinct_coeff("U")
    0.6
    >>> get_extinct_coeff("other")
    Traceback (most recent call last):
        ...
    ValueError: Bad filter name: other
    """
    coeff_values = {"none": 0.2, "U": 0.6, "B": 0.4, "V": 0.2, "R": 0.1, "I": 0.08}
    try:
        return coeff_values[filter_name]
    except KeyError:
        raise ValueError("Bad filter name: {filter_name}".format(filter_name=filter_name))


def mag_zeropoint(filter_name):
    """
    Given a filter name, return the number of photons per square centimeter per second
    which a star of zero'th magnitude would produce above the atmosphere. We assume that the star
    has a spectrum like Vega's.

    The numbers are pre-calculated and we just pick the appropriate one for the given filter.

    :param filter_name: Filter name string
    :return: Number of photons per sq.cm. per sec from zero'th mag star

    >>> mag_zeropoint("none")
    4320000.0
    >>> mag_zeropoint("B")
    391000.0
    >>> get_extinct_coeff("other")
    Traceback (most recent call last):
        ...
    ValueError: Bad filter name: other
    """
    photons_per_sq_cm = {"none": 4.32e+06, "U": 5.50e+05, "B": 3.91e+05, "V": 8.66e+05,
                         "R": 1.10e+06, "I": 6.75e+05}
    try:
        return photons_per_sq_cm[filter_name]
    except KeyError:
        raise ValueError("Bad filter name: {filter_name}".format(filter_name=filter_name))


def fraction_inside_cython(fwhm, radius, pixsize, max_pix_rad=30, piece=20):
    """
    Figure out what fraction of a star's light falls within he aperture.  We assume that
    the starlight has a circular gaussian distribution with FWHM given by the first argument
    (with units of arcsec).

    This function goes to the trouble of calculating how much of the light falls within
    fractional pixels defined by the given radius of a synthetic aperture.

    Because this routine does a lot of loops the main code is written in cython for increased
    speed. See aperture_sum.pyx.

    :type piece: int
    :type max_pix_rad: int
    :param fwhm: FWHM in arcsec
    :param radius: radius of aperture, in arcsec
    :param pixsize: size of a pixel, in arcsec
    :param max_pix_rad: maximum pixel radius
    :param piece: how many pieces do we sub-divide pixels into
    :return: the fraction of light within aperture, the fraction of light in the maximum pixel
    """
    import pyximport
    pyximport.install()
    import aperture_sum

    if pixsize <= 0.0:
        raise ValueError("Radius must be greater than zero")

    fwhm /= pixsize
    radius /= pixsize

    if radius > max_pix_rad:
        raise ValueError("Radius exceeds limit of {limit}".format(limit=max_pix_rad))

    return aperture_sum.aperture_sum(fwhm, radius, max_pix_rad, piece)


def fraction_inside(fwhm, radius, pixsize, max_pix_rad=None, piece=None):
    """
    Figure out what fraction of a star's light falls within the aperture. We assume that
    the starlight has a circular gaussian distribution with FWHM given by the first argument
    (with units of arcsec).  We calculate the fraction of that light which falls within an aperture
    of radius given by second argument (with units of arcsec).

    This is only really needed if the cython version above is not available

    :param fwhm: FWHM in arcsec
    :param radius: radius of aperture, in arcsec
    :param pixsize: pixel size in arsec
    :param max_pix_rad: not used but included so that both fraction_inside functions
                        are called the same
    :param piece: not used but included so that both fraction_inside functions are called the same
    :return: The fraction of light within aperture and the fraction of light in the maximum pixel
    """
    import math

    sigma = fwhm / 2.35
    z = radius / (sigma * 1.414)
    x1 = math.erf(z)
    ratio = (x1 * x1)

    z_pixel = (pixsize * 0.5) / (sigma * 1.414)
    x1_pixel = math.erf(z_pixel)
    ratio_pixel = (x1_pixel * x1_pixel)

    max_pixel_fraction = ratio * ratio_pixel
    if max_pixel_fraction > 1.0:
        max_pixel_fraction = 1.0
    return ratio, max_pixel_fraction


def check_values(**kwargs):
    """
    Check the values of signal_to_noise. Values read from kwargs so the function can be easily
    extended if the parameters change.

    :param kwargs: Dictionary of key value pairs
    :raises: ValueError if any of the values are not within the correct ranges
    """
    for key, value in kwargs.iteritems():
        if value is None:
            raise ValueError("Missing value for {item}".format(item=key))

    try:
        if (kwargs['mag_start'] > kwargs['mag_end']) or kwargs['dmag'] < 0.0:
            raise ValueError("Invalid magnitude start/end or delta values")

        if kwargs['tel_diam'] <= 0.0:
            raise ValueError("Invalid telescope diameter, must be greater than 0")

        if kwargs['pixsize'] <= 0.0:
            raise ValueError("Invalid pixel size, must be greater than 0")

        if kwargs['readnoise'] < 0.0:
            raise ValueError("Invalid readnoise, must be >= 0")

        if kwargs['obs_diam'] < 0.0:
            raise ValueError("Invalid obscuration diameter, must be >= 0")

        if kwargs['exptime'] <= 0.0:
            raise ValueError("Invalid exposure time, must be greater than 0")

        if kwargs['fwhm'] <= 0.0:
            raise ValueError("Invalid FWHM, must be greater than 0")

        if kwargs['aper_rad'] <= 0.0:
            raise ValueError("Invalid aperture radius, must be greater than 0")

        if kwargs['sat_level'] <= 0.0:
            raise ValueError("Invalid saturation level, must be greater than 0")

        if kwargs['gain'] <= 0.0:
            raise ValueError("Invalid gain, must be greater than 0")

        if kwargs['airmass'] < 1.0:
            raise ValueError("Invalid airmass, can not be less than 1")

        if not 0.0 <= kwargs['qe'] < 1.0:
            raise ValueError("Invalid QE, must be between 0.0 and 1.0")

        if not 0.0 <= kwargs['skymag'] < 100.0:
            raise ValueError("Invalid skymag, must be between 0.0 and 100.0")

    except KeyError as missing:
        raise ValueError("Missing key {item}".format(item=missing))


def signal_to_noise(filter_name=None, mag_start=None, mag_end=None, dmag=None, tel_diam=None,
                    qe=None, readnoise=None, pixsize=None, skymag=None, airmass=None,
                    exptime=None, fwhm=None, aper_rad=None, obs_diam=0.0, sat_level=None, gain=None,
                    fraction_inside_func=fraction_inside_cython):
    """
    Calculates the signal to noise for a range of magnitudes supplied as mag_start, mag_end and
    dmag for magnitudes steps. The function performs the following steps:

    Checks for validity of values supplied by user
    Convert readnoise and saturation levels to electrons from ADU
    Prepare to perform calculations
    Calculates how many pixels are inside the aperture
    Calculates what fraction of the star's light falls within aperture
    Calculate the maximum fraction of light in a pixel

    Then for each magnitude step:
    Calculate the number of electrons collected on the CCD from star, total
    Decrease the number of electrons from star, due to extinction
    Calculate number of electrons collected on CCD from sky, per pixel
    Calculate the total number of electrons from star inside aperture
    Calculate if the star is saturated
    Calculate the total number of electrons from sky inside aperture
    Calculate the total variance from readout noise within the aperture
    Calculate signal-to-noise

    :param filter_name: Filter name (None, U, B, V, R, I)
    :param mag_start: Magnitude to start calculations
    :param mag_end: Magnitude to end calculations
    :param dmag: Magnitude steps for calculations
    :param tel_diam: Telescope diameter in cm
    :param obs_diam: Obscuration diameter in cm
    :param qe: QE primarily determined by CCD (0.3 cheap CCD, 0.7 expensive CCD)
    :param readnoise: CCD readout noise in ADU
    :param sat_level: Saturation level in ADU
    :param gain: Gain (electrons / ADU)
    :param pixsize: Pixel size (arcsec/pixel)
    :param skymag: Sky brightness: (mag per sq. arcsec) Suburbs 17, mountains 21
    :param airmass: Airmass
    :param exptime: Exposure time in seconds
    :param fwhm: FWHM in arcsec
    :param aper_rad: Radius for photometry area over which signal is measured in arcsec
    :param fraction_inside_func: Functions used for calculating star light inside the aperture
    :return: List of named Result tuples
    """
    import math

    check_values(mag_start=mag_start, mag_end=mag_end, dmag=dmag, tel_diam=tel_diam, qe=qe,
                 readnoise=readnoise, pixsize=pixsize, skymag=skymag, airmass=airmass,
                 exptime=exptime, fwhm=fwhm, aper_rad=aper_rad, obs_diam=obs_diam, gain=gain,
                 sat_level=sat_level)

    readnoise *= gain
    sat_level *= gain

    extinct_coeff = get_extinct_coeff(filter_name)
    nphoton = mag_zeropoint(filter_name)

    npix = (math.pi * aper_rad * aper_rad) / (pixsize * pixsize)

    fraction, max_pixel_fraction = fraction_inside_func(fwhm, aper_rad, pixsize)

    effective_area = math.pi * 0.25 * (tel_diam ** 2 - obs_diam ** 2)
    reduction_factor = effective_area * exptime * qe * nphoton
    extinction = airmass * extinct_coeff

    mag = mag_start

    output = []

    while mag <= mag_end:
        star_electrons = math.pow(10.0, -0.4 * mag) * reduction_factor
        star_electrons *= math.pow(10.0, -0.4 * extinction)

        sky_electrons_per_pix = math.pow(10.0, -0.4 * skymag) * reduction_factor
        sky_electrons_per_pix *= pixsize * pixsize

        max_pixel = star_electrons * max_pixel_fraction

        if max_pixel * gain > sat_level:
            saturated = True
        else:
            saturated = False

        star_electrons *= fraction
        sky_electrons = sky_electrons_per_pix * npix
        read_electrons = readnoise * readnoise * npix

        signal = star_electrons

        noise = math.sqrt(read_electrons + sky_electrons + star_electrons)
        signal_noise = signal / noise

        output.append(Result(magnitude=mag, star_adu=star_electrons / gain,
                             sky_adu=sky_electrons / gain, read_adu=read_electrons / gain,
                             signal_to_noise=signal_noise, saturated=saturated))

        mag += dmag

    return output