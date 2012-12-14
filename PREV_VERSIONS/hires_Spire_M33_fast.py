# This  is SPIRE-specific code for HIRES

import hires_v1_07e as hires


import pyfits


hires.set_FLUX_UNITS('Jy/beam')



#====================================================================
# Override hires.read_all_IN_files method

def read_all_Spire_files():
    ''' Read SPIRE input tables
        and return a list of Sample objects
    '''
    raHduList = pyfits.open('m33_raTable.fits')
    raHdu = raHduList[1]
    names = []
    for det in raHdu.data.names:
        if (det[:3] == 'PSW' and det[3] != 'T' and det[3] != 'R' and det[4] != 'P'):
        	names.append(det)
    names.sort()
    decHduList = pyfits.open('m33_decTable.fits')
    signalHduList = pyfits.open('m33_signalTable.fits')
    maskHduList = pyfits.open('m33_maskTable.fits')
    
    samples = []
    
    crval1 = hires.get_FITS_keyword('CRVAL1')
    crval2 = hires.get_FITS_keyword('CRVAL2')
    projection = hires.Gnomonic(crval1, crval2)

    for det in names:
    	#ra = raHduList[1].data[det]
    	ra = raHduList[1].data.field(det)[::17]
    	dec = decHduList[1].data.field(det)[::17]
    	sig = signalHduList[1].data.field(det)[::17]
    	mask = maskHduList[1].data.field(det)[::17]
    	inx = mask <= 1024
    	hires.log(3, 'Generating samples for detector %s',det)
    	x, y = projection.lonlat2xy(ra[inx], dec[inx])
    	for i in range(len(ra[inx])):
    	    #x, y = projection.lonlat2xy(ra[inx][i], dec[inx][i])
    	    samp = hires.Sample(-x[i], -y[i], sig[inx][i], 1, 1.0)
    	    if (i%3333 == 0):
    	        hires.log(3, 'Created %dth sample for detector %s, x=%f, y=%f, ra=%f, dec=%f',\
    	        	i,det,x[i], y[i], ra[inx][i],dec[inx][i]) 
    	    samples.append(samp)
    	
    return samples


#====================================================================
# Override hires.read_all_DRF_files method

def read_all_Spire_beams():
    ''' Read SPIRE beams
        and return a dictionary of response functions
    '''
    global DRF_SET_ID
    DRF_SET_ID = 'single'
    raHduList = pyfits.open('m33_raTable.fits')
    raHdu = raHduList[1]
    names = []
    for det in raHdu.data.names:
        if (det[:3] == 'PSW'):
        	names.append(det)
    detHduList = pyfits.open('../Spire/params/psw_beam_1arcsec_cutout.fits')
    drf_array = detHduList[1].data
    naxis1 = detHduList[1].header['NAXIS1']
    naxis2 = detHduList[1].header['NAXIS2']
    deg_per_pix = 1./3600.
    radius_pix = naxis1 / 2
    radius_degrees = radius_pix * deg_per_pix    	
    detectors = {}
    #print names
    #
    def xy2response_function(x, y):
    	''' interpolate in detector response array
       		  x, y in degrees (relative to DRF center
    	'''
    	iFloat = (x + radius_degrees) / deg_per_pix
    	jFloat = (y + radius_degrees) / deg_per_pix
    	# cheap "interpolation" -- just takes nearest one
        iInt = int(iFloat+0.5)
        jInt = int(jFloat+0.5)
        if iInt<0 or iInt>=naxis1 or jInt<0 or jInt>=naxis2 :
            response = 0.0
        else:
            response = drf_array[iInt, jInt]
        return response
    #
    #for id in names:
    #    detector = hires.Detector(id, radius_degrees, xy2response_function)
    #    detectors[id] = detector
    detectors[1] = hires.Detector(1, radius_degrees, xy2response_function, DRF_SET_ID)
    
    return detectors



hires.read_all_IN_files = read_all_Spire_files
hires.read_all_DRF_files = read_all_Spire_beams


#====================================================================
# Call hires
if __name__ == "__main__":
    hires.main()
