# This  is Planck-specific code for HIRES
# Converts flux to RJ

import hires
# This points to development version of hires

import pyfits
import numpy as np

#====================================================================
# define detector specific information
# K_RJ = (K_CMB - Offset)/Scaling
det_offset_scale_rms = {}
det_offset_scale_rms['1'] = ( 33.858, 10351.1, 1.086 )
det_offset_scale_rms['2'] = ( 34.189,  9590.0, 1.097 )
det_offset_scale_rms['3'] = ( 33.840, 10640.6, 1.085 )

#====================================================================
# Override hires.read_one_IN_file method

hires.set_FLUX_UNITS('K_RJ')

def read_Planck_toast_file(filename, samples):
    ''' Read a Planck TOAST input file
        and append a SampleSet object to samples[] list
    '''
    detector_id = (filename.split('-')[-1])[0:-5]
    offset_scale_rms = det_offset_scale_rms[detector_id]
    offset = offset_scale_rms[0]
    scale  = offset_scale_rms[1]

    hdulist = pyfits.open(filename)
    phi = hdulist[1].data.field('PHI')
    theta = hdulist[1].data.field('THETA')
    signal = ( hdulist[1].data.field('SIGNAL') - offset ) / scale
    rtod = 180.0/3.14159265
    glon = phi * rtod
    glat = 90.0 - theta * rtod
    crval1 = hires.get_FITS_keyword('CRVAL1')
    crval2 = hires.get_FITS_keyword('CRVAL2')
    projection = hires.Gnomonic(crval1, crval2)

    det_id = '1'
    angle = 999.0 # Mark to be computed by HIRES

    for llist, blist, slist in zip(glon, glat, signal):
        xlist, ylist = projection.np_lonlat2xy(llist, blist)
        # negate xlist because CDELT1 is negative
        # angle = None forces estimation from delta x,y
        samp = hires.SampleSet(-xlist, ylist, slist, det_id, angle=None)
        samples.append(samp)

hires.read_one_IN_file = read_Planck_toast_file


#====================================================================
# Call hires
if __name__ == "__main__":
    hires.main()
