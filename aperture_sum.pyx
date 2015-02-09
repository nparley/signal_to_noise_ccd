cimport cython
from libc.math cimport exp

@cython.cdivision(True)
@cython.boundscheck(False)
def aperture_sum(double fwhm, double radius,int max_pix_rad, int piece):
    """
    Figure out what fraction of a star's light falls within he aperture.  We assume that
    the starlight has a circular gaussian distribution with FWHM given by the first argument
    (with units of arcsec).

    This function goes to the trouble of calculating how much of the light falls within
    fractional pixels defined by the given radius of a synthetic aperture.

    This function is called from functions.fraction_inside_cython which does a few checks on the
    data first

    :param fwhm: FWHM in pixel units
    :param radius: radius in pixel units
    :param max_pix_rad: maximum pixel radius
    :param piece: how many pieces do we sub-divide pixels into
    :return:the fraction of light within aperture, the fraction of light in the maximum pixel
    """
    cdef double psf_center_x = 0.5
    cdef double psf_center_y = 0.5

    cdef double sigma2 = (fwhm / 2.35) ** 2
    cdef double radius2 = radius ** 2
    cdef double bit = 1.0 / piece

    cdef double rad_sum = 0.0
    cdef double all_sum = 0.0
    cdef double max_pix_sum = 0.0

    cdef int start = int(0.0 - max_pix_rad)
    cdef int stop = int(max_pix_rad)

    cdef int i, j, k, l
    cdef double pix_sum, x, fx, y, fy, inten, this_bit, rad2

    for i in range(start, stop):
        for j in range(start, stop):
            pix_sum = 0.0

            for k in range(0, piece):
                x = (i - psf_center_x) + (k + 0.5) * bit
                fx = exp(-(x * x) / (2.0 * sigma2))

                for l in range(0, piece):
                    y = (j - psf_center_y) + (l + 0.5) * bit
                    fy = exp(-(y * y) / (2.0 * sigma2))
                    inten = fx * fy
                    this_bit = inten * bit * bit
                    pix_sum += this_bit
                    rad2 = x * x + y * y
                    if rad2 <= radius2:
                        rad_sum += this_bit

            all_sum += pix_sum
            if pix_sum > max_pix_sum:
                max_pix_sum = pix_sum

    cdef double ratio = rad_sum / all_sum
    cdef double max_pixel_fraction = max_pix_sum / all_sum
    return ratio, max_pixel_fraction


