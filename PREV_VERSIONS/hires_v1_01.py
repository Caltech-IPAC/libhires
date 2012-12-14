# "HIRES" Python program to produce high resolution images
# Creates image using "maximum correlation method" (smoothest image that is
# consistent with the input data and the detector response functions)
# ==> Must override the "read_one_IN_file" method

PROGRAM = 'HIRES'
VERSION = 'v1_01'

#@======= (N) Notes =================
#   i,j  in pixels  0<=i<=NPIX
#   u,v in pixels from center of image (may not be center of pixel)
#   x,y in degrees from center of image

#   i = (x / DEG_PER_PIX) + NPIXi/2
#   j = (y / DEG_PER_PIX) + NPIXj/2

#   x = (i - NPIXi/2) * DEG_PER_PIX
#   y = (j - NPIXj/2) * DEG_PER_PIX

#====================================================================
# imports
#====================================================================
import numpy as np
import pyfits
import sys, os
import exceptions
import glob
import time
import logging
from math import radians, cos, sin

#====================================================================
# initialize logging
#====================================================================
logging.basicConfig(
  filename='log_hires.log',
  filemode='w',
  format='%(asctime)s %(message)s',
  datefmt='%Y/%m/%d %H:%M:%S',
  level=logging.INFO  )
#  level=logging.ERROR  )
#  level=logging.WARNING  )

LAST_TIME = time.clock()
def print_completed(msg):
   ''' print and log message, plus delta time(seconds) since
       last call (intended to be called after each step) '''
   global LAST_TIME
   t = time.clock()
   delta = t - LAST_TIME
   LAST_TIME = t
   full_msg =  str(delta) +' Completed '+ msg
   logging.info(full_msg)
   print full_msg


#====================================================================
class Detector:
#====================================================================
    ''' one instance for each detector, generally one DRF file '''

    all = {}

    def __init__(self, id, xy_coverage_degrees, response_function):
        self.id   = id # detector ID (integer)
        self.xy_coverage_degrees = xy_coverage_degrees
        self.response_function = response_function
        Detector.all[id] = self

#====================================================================
class Sample:
#====================================================================
    ''' one instance for each sample in the IN files '''

    def __init__(self, x, y, flux, detector_id, angle, noise=1.0):
        self.x           = x # delta from center (degrees)
        self.y           = y # delta from center (degrees)
        self.flux        = flux # signal in units that will go into image
        self.detector_id = detector_id # detector ID (string)
        self.angle       = angle # orientation degrees east of north
        self.noise       = noise # same units as flux

#====================================================================
class Footprint:
#====================================================================
    ''' A Footprint is a mapping of the detector response function
    on the image.  There is one footprint for each sample (unless out of range).
    The Footprint is basically an array of responses for each image pixel
    near the sample.
      Note: i,j "bounds" are ~Python style bounds e.g. i1 = max i + 1
    '''

    all_response_arrays = {} # indexed by (detID, iFracID, jFracID, angleID)

    def __init__(self, detector_id, iFloat, jFloat, angle, flux):
        self.flux = flux
        iWhole, iFrac = divmod(iFloat, 1.0)
        jWhole, jFrac = divmod(jFloat, 1.0)
        full_response_array = self.get_response_array(detector_id, iFrac, jFrac, angle)
        untrimmed_image_bounds = \
                          self.compute_image_bounds(full_response_array, iWhole, jWhole)
        trimmed_image_bounds = self.trim_image_bounds(untrimmed_image_bounds)
        foot_sum, trimmed_response_array = self.trim_response_array( \
                         full_response_array, untrimmed_image_bounds, trimmed_image_bounds)
        self.image_bounds = trimmed_image_bounds
        self.response_array = trimmed_response_array
        self.foot_sum = foot_sum

        
    #------------------------------------------------------------------------
    def get_response_array(self, detector_id, iFloat, jFloat, angle):
        ''' see if appropriate footprint array has already been generated
            If so, return it, otherwise call the method to generate it
        '''
        # compose key (to use previously computed footprint array)
        angle_round = round( angle / ANGLE_TOLERANCE )
        angleID = int( angle_round )
        angle = angle_round * ANGLE_TOLERANCE
        iFrac = iFloat % 1.0 # range:  0.0 <= frac < 1.0
        jFrac = jFloat % 1.0
        if FOOTPRINTS_PER_PIX <= 100:
            iFracID = int( iFrac * FOOTPRINTS_PER_PIX ) # range: 0 <= ID < FootPerPix
            jFracID = int( jFrac * FOOTPRINTS_PER_PIX )
            footArr_key =  ( detector_id, angleID, iFracID, jFracID )
            if footArr_key in Footprint.all_response_arrays:
                footArr = Footprint.all_response_arrays[footArr_key]
            else: 
                delta = 1.0 / float(FOOTPRINTS_PER_PIX)
                zero = (delta / 2.0 ) - 0.5
                iFracMod = zero + iFracID * delta # rounded for all within tolerance
                jFracMod = zero + jFracID * delta
                footArr =self.generate_response_array(detector_id, iFracMod,jFracMod, angle)
                Footprint.all_response_arrays[footArr_key] = footArr
        else:
            footArr = self.generate_response_array(detector_id, iFrac, jFrac, angle)

        return footArr

    #------------------------------------------------------------------------
    def generate_response_array(self, detector_id, iFrac, jFrac, angle):
        ''' generate full footprint array for given detector at
            given angle and centered at given x,y offset from pixel center
         detector_id
         angle (degrees) 
         iFrac, jFrac offset (frac of pixel) from center of pixel -.5 <= Frac < 0.5
        '''
        detector = Detector.all[detector_id]
        x0, x1, y0, y1 = detector.xy_coverage_degrees
        response_func = detector.response_function
        u0 = int( x0 / DEG_PER_PIX ) # (will do an implicit "floor", towards zero)
        u1 = int( x1 / DEG_PER_PIX )
        v0 = int( y0 / DEG_PER_PIX )
        v1 = int( y1 / DEG_PER_PIX )
        xStart = ( float(u0) - iFrac) * DEG_PER_PIX
        yStart = ( float(v0) - jFrac) * DEG_PER_PIX
        iNPIX = u1 - u0 + 1
        jNPIX = v1 - v0 + 1

        # create the response array
        foot_array = np.empty( (iNPIX, jNPIX), dtype=np.float32 )
        x = xStart
        for i in range(iNPIX):
            y = yStart
            for j in range(jNPIX):
                response_here = response_func(x,y)
                foot_array[i, j] = response_here
                y += DEG_PER_PIX
            x += DEG_PER_PIX

        # normalize to sum of 1.0 and return
        foot_sum = foot_array.sum()
        normalized_foot_array = foot_array / foot_sum # normalize
        return normalized_foot_array

    #------------------------------------------------------------------------
    def compute_image_bounds(self, full_reaponse_array, iWhole, jWhole):
        ''' Compute where this foorprint lands in the image.
            Could lie partially off the image but will be trimmed later.
        '''
        iSize, jSize = full_reaponse_array.shape
        i0 = int(iWhole) - iSize/2
        i1 = i0 + iSize
        j0 = int(jWhole) - jSize/2
        j1 = j0 + jSize
        return (i0, i1, j0, j1)

    #------------------------------------------------------------------------
    def trim_image_bounds(self, untrimmed_bounds):
        ''' Thim the i.j bounds to lie completely within the image. '''
        i0 = max(untrimmed_bounds[0], 0)
        i1 = min(untrimmed_bounds[1], NPIXi)
        j0 = max(untrimmed_bounds[2], 0)
        j1 = min(untrimmed_bounds[3], NPIXj)
        return (i0, i1, j0, j1)

    #------------------------------------------------------------------------
    def trim_response_array(self, full_array, untrimmed_bounds, trimmed_bounds):
        ''' Thim the actual footprint response array if image i,j bounds were trimmed.
            Just return full array if image bounds were not trimmed.
        '''
        if untrimmed_bounds== trimmed_bounds: return (1.0, full_array)
        iSize, jSize = full_array.shape
        i0 = trimmed_bounds[0] - untrimmed_bounds[0]
        i1 = trimmed_bounds[1] - untrimmed_bounds[1] + iSize 
        j0 = trimmed_bounds[2] - untrimmed_bounds[2]
        j1 = trimmed_bounds[3] - untrimmed_bounds[3] + iSize 
        trimmed_arr = full_array[i0:i1, j0:j1]
        arr_sum = trimmed_arr.sum()
        return ( arr_sum, trimmed_arr )


#====================================================================
def create_all_footprints(setOfSamples, setOfDetectors):
#====================================================================
    ''' Create footprints for all samples in setOfSamples
        Use detector (DRF definitions) info given in seOfDetectors
    '''
    all_footprints = [] ## 0=image_bounds, 1=array, 2=sum(array)  3=flux
    n_rejected = 0
    HALF_PIXi = float(NPIXi-1)/2.0
    HALF_PIXj = float(NPIXj-1)/2.0
    for samp in setOfSamples:
        # (i,j are image pixel coordinates (range 0 to NPIX-1))
        iFloat =  samp.x / DEG_PER_PIX + HALF_PIXi
        jFloat =  samp.y / DEG_PER_PIX + HALF_PIXj
        i = int(round( iFloat ))
        j = int(round( jFloat ))

        if i>=0 and i<NPIXi and j>=0 and j<NPIXj:
            footprint = Footprint(samp.detector_id, iFloat, jFloat, samp.angle, samp.flux)
            all_footprints.append(footprint)
        else: n_rejected += 1

    print_completed( 'footprint creation: ' +str(len(all_footprints))+ ' Footprints' )
    print 'Rejected', n_rejected, 'Samples and created', \
     len(Footprint.all_response_arrays), 'full response arrays'
    return all_footprints


#====================================================================
def read_start_image():
#====================================================================
    ''' Read in a starting image '''
    global ITER_START
    hduList = pyfits.open(STARTING_IMAGE)
    kwd = hduList[0].header
    if kwd.has_key('ITERNUM'):
        ITER_START = kwd['ITERNUM']
        print 'found ITERNUM in header =', ITER_START
    flux_image = hduList[0].data
    return flux_image

#====================================================================
def calc_wgt_array(footprints):
#====================================================================
    ''' Create an image of all footprints *1.0 '''
    wgt_array = np.zeros( (NPIXi, NPIXj), dtype=np.float32)
    wgt_array += 0.00001 # to avoid NaNs in results
    for foot in footprints:
        i0, i1, j0, j1 = foot.image_bounds
        wgt_array[i0:i1, j0:j1] += foot.response_array
    print_completed('calc_wgt_array')
    return wgt_array

#====================================================================
def calc_flux_wgt_array(footprints):
#====================================================================
    ''' Create an image of all footprints * flux '''
    flux_wgt_array = np.zeros( (NPIXi, NPIXj), dtype=np.float32)
    for foot in footprints:
        i0, i1, j0, j1 = foot.image_bounds
        flux_wgt_array[i0:i1, j0:j1] += foot.response_array * foot.flux
    print_completed('calc_flux_wgt_array')
    return flux_wgt_array

#====================================================================
def calc_corr_wgt_array(footprints, flux_array, iter_num):
#====================================================================
    ''' Create a correction image: For each footprint:
         1. Calculate "correcton" 
         2. Add into correction image: response_array * correction
    '''
    corr_wgt_array = np.zeros( (NPIXi, NPIXj), dtype=np.float32)
    for foot in footprints:
        i0, i1, j0, j1 = foot.image_bounds
        integration = foot.response_array * flux_array[i0:i1, j0:j1]
        flux_prime = integration.sum() / foot.foot_sum
        correction = foot.flux / flux_prime
        corr_wgt_array[i0:i1, j0:j1] += foot.response_array * correction

    min_c = np.nanmin(corr_wgt_array)
    max_c = np.nanmax(corr_wgt_array)
    std_c = np.std(corr_wgt_array)
    stats = 'iter=%3d min=%.6f max=%.6f stdev=%.6f' % (iter_num, min_c, max_c, std_c)
    print_completed('calc_corr_wgt_array ' + stats)

    return corr_wgt_array


#====================================================================
def get_paramaters(args):
#====================================================================
    ''' Process user parameters 
        Everything put into global variables (ALL CAPS).
        Global variables should generally only be set in this method. '''

    global INFILE_PREFIX # directory and prefix of files containing samples
    global OUTFILE_PREFIX # directory and prefix for output filenames
    global DRF_PREFIX # directory and prefix of files containing response func images
    global DEG_PER_PIX # image scale (CDELT) degrees per pixel
    global NPIXi, NPIXj # image size (pixels)
    global ITER_START # starting iteration number (0=initial, 1=after 1st correction)
    global ITER_MAX # max iteration nmber to compute
    global ITER_LIST # which iterations to write outupt files for
    global STARTING_IMAGE # user supplied starting image
    global ANGLE_TOLERANCE # footprints re-used if within this tollerance (degrees)
    global FOOTPRINTS_PER_PIX # generate footprints at 1/N pixel accuracy
    global FLUX_UNITS # for FITS keyword BUNIT in flux images
    global FITS_KEYWORDS # FITS keywords to be added to output files

    FITS_KEYWORDS = [] 
    ITER_LIST = {}
    STARTING_IMAGE = None
    ITER_START = 0
    FOOTPRINTS_PER_PIX = 1

    INFILE_PREFIX  = args[1]
    OUTFILE_PREFIX = args[2]
    param_files    = ( args[3:] )
    
    for filename in ( param_files ):
        for line in open(filename):
            par_and_comment = line.strip().split('#')
            par = par_and_comment[0]
            if len(par)==0: continue # just comments on this line
            words = par.split()
            if len(words)<2:
                print 'ERROR: No value given for parameter:', par, 'in file:', filename
                exit(-2)
            name = words[0]
            val = words[1]
            if   name == 'SIZE_NPIX':
               npix_split =  val.split(',')
               NPIXi = int(npix_split[0])
               if len(npix_split)==1: NPIXj = NPIXi
               else: NPIXj = int(npix_split[1])
            elif name == 'ARCSEC_PER_PIX': arcsec_per_pix = float(val)
            elif name == 'FLUX_UNITS': FLUX_UNITS = val
            elif name == 'ANGLE_TOLERANCE': ANGLE_TOLERANCE = float(val)
            elif name == 'FOOTPRINTS_PER_PIX': FOOTPRINTS_PER_PIX = int(val)
            elif name == 'DRF_PREFIX': DRF_PREFIX = val
            elif name == 'STARTING_IMAGE': STARTING_IMAGE = val
            elif name == 'ITER_MAX': ITER_MAX = int(val)
            elif name == 'ITER_LIST':
               iter_str =  val.split(',')
               for n_str in iter_str: ITER_LIST[int(n_str)] = 1
            elif name == 'KWD':
                if len(words) >=3: comment = ' '.join(words[3:])
                else: comment = None
                kwd = val
                val = numerify(words[2])
                FITS_KEYWORDS.append( ( kwd, val, comment ) )
            elif name == '': pass
            else:
                print 'ERROR: Unknown parameter name:', name, 'in file:', filename
                exit(-2)

    # compute globals from params
    DEG_PER_PIX = arcsec_per_pix / 3600.0
    ITER_LIST[ITER_MAX] = 1

    # print all params
    print 'hires.py VERSION:', VERSION
    print 'INFILE_PREFIX =', INFILE_PREFIX
    print 'OUTFILE_PREFIX =', OUTFILE_PREFIX
    print 'DRF_PREFIX =', DRF_PREFIX
    print 'ARCSEC_PER_PIX =', arcsec_per_pix
    print 'FLUX_UNITS =', FLUX_UNITS
    print 'ANGLE_TOLERANCE =', ANGLE_TOLERANCE
    print 'FOOTPRINTS_PER_PIX =', FOOTPRINTS_PER_PIX
    print 'ITER_START =', ITER_START
    print 'ITER_MAX =', ITER_MAX
    print 'ITER_LIST =', sorted( ITER_LIST.keys() )
    print 'NPIXi =', NPIXi
    print 'NPIXj =', NPIXj
    print 'STARTING_IMAGE =', STARTING_IMAGE
    print 'FITS keywords:', FITS_KEYWORDS


#====================================================================
def get_FITS_keyword(keyword):
#====================================================================
    ''' Get a FITS keyword from FITS_KEYWORDS.  Return None if not there. '''
    for k in FITS_KEYWORDS:
        kwd = k[0]
        if kwd == keyword: return k[1]
    return None # not found

#====================================================================
def numerify (s):
#====================================================================
    ''' convert string to int or float if possible '''
    try:
        return int(s)
    except exceptions.ValueError:
        try:
            return float(s)
        except exceptions.ValueError:
            return s

#====================================================================
def write_FITS_image(pixel_array, file_type, iter=None):
#====================================================================
    ''' Writes out one FITS image of specified type.
    '''

    filename = OUTFILE_PREFIX + '_' + file_type
    if iter is not None: filename += '_' + str(iter)
    filename += '.fits'

    directory = os.path.dirname(OUTFILE_PREFIX)
    if not os.path.exists(directory): os.makedirs(directory)

    hdu = pyfits.PrimaryHDU(pixel_array)
    prihdr = hdu.header

    # add keywords from user parameters
    for kwd_tuple in FITS_KEYWORDS:
        prihdr.update(kwd_tuple[0], kwd_tuple[1], comment=kwd_tuple[2])

    cdelt_rounded = float('%.7f' % DEG_PER_PIX)
    prihdr.update('CDELT1', -cdelt_rounded, comment='degrees per pixel')
    prihdr.update('CDELT2',  cdelt_rounded, comment='degrees per pixel')
    prihdr.update('CRPIX1', (NPIXi+1)/2 , comment='center pixel')
    prihdr.update('CRPIX2', (NPIXj+1)/2 , comment='center pixel')

    if (file_type=='flux'): prihdr.update('BUNIT', FLUX_UNITS)
    if   file_type=='flux': t_comment = 'HIRES flux image'
    elif file_type=='cov':  t_comment = 'HIRES coverage image'
    else:                   t_comment = 'HIRES image'
    prihdr.update('FILETYPE', file_type, comment=t_comment)

    if DRF_PREFIX is not None:
        drf_pref = DRF_PREFIX.split('/')[-1]
        prihdr.update('DRF'     , drf_pref  , comment='Detector Response Files used')  

    # add keywords from processing
    if iter is not None: prihdr.update('ITERNUM', iter, comment='HIRES iteration number')

    filename_only = filename.split('/')[-1] # strip off leading directory names
    prihdr.update('FILENAME', filename_only, comment='name of this file')  
    dt = time.strftime("%Y/%m/%d %H:%M:%S")
    prihdr.update('DATE'     , dt  , comment='when this file was created')  
    prihdr.update('CREATED', PROGRAM+' '+VERSION, comment='software version that created this file')

    hdulist = pyfits.HDUList([hdu])
    if os.path.isfile(filename): os.remove(filename)
    hdulist.writeto(filename)
    print_completed('write FITS file: '+ filename)

#====================================================================
def read_one_DRF_file(filename):
#====================================================================
    ''' Read one detector response function file.
        Return the detector ID and a Detector object.
        This method may be overiden for project-specific DRF files.
    '''
    hduList = pyfits.open(filename)
    kwd = hduList[0].header
    id = kwd['DETECTOR']
    naxis1 = kwd['NAXIS1']
    naxis2 = kwd['NAXIS2']
    crpix1 = kwd['CRPIX1']
    crpix2 = kwd['CRPIX2']
    deg_per_pix = kwd['CDELT1']
    x0 = (1.0 - crpix1) * deg_per_pix
    x1 = (float(naxis1) - crpix1) * deg_per_pix
    y0 = (1.0 - crpix2) * deg_per_pix
    y1 = (float(naxis2) - crpix2) * deg_per_pix
    drf_array = hduList[0].data

    def function_interpolate(x, y): # interpolate in detector response array
        iFloat = (x - x0) / deg_per_pix
        jFloat = (y - y0) / deg_per_pix
        iWhole, iFrac = divmod(iFloat, 1.0)
        jWhole, jFrac = divmod(jFloat, 1.0)
        iInt = int(iWhole)
        jInt = int(jWhole)
        p00 = drf_array[iInt  , jInt  ] * (1.0-iFrac) * (1.0-jFrac)
        p01 = drf_array[iInt  , jInt+1] * (1.0-iFrac) *      jFrac 
        p10 = drf_array[iInt+1, jInt  ] *      iFrac  * (1.0-jFrac)
        p11 = drf_array[iInt+1, jInt+1] *      iFrac  *      jFrac 
        response = p00 + p01 + p10 + p11
        return response

    print "Detector: ", id, 'file=', filename.split('/')[-1], \
            'x:', x0*60.0, x1*60.0, 'y:', y0*60.0, y1*60.0, '(arcmin)'
    detector = Detector(id, (x0, x1, y0, y1), function_interpolate)
    return ( id, detector )


#====================================================================
def read_all_DRF_files():
#====================================================================
    ''' Read all of the RDF (detector response function) files.
        (by calling the "read_one_DRF_file" method)
    '''
    all_detectors = {} # dictionay accessed by detector_id
    file_list = glob.glob(DRF_PREFIX + '*')
    for file in file_list:
        id, detector = read_one_DRF_file(file)
        all_detectors[id] = detector
    print_completed('setup of ' + str(len(all_detectors)) +' detector(s)')
    return all_detectors

#====================================================================
class Gnomonic:
#====================================================================
    ''' Compute gnomonic (tangent plane) projection of lon,lat to x,y
        lat,lon,x,y all in degrees
        This method intended to be used by overrides of "read_one_IN_file".
        Equations from: http://mathworld.wolfram.com/GnomonicProjection.html
        Example:
          gn = Gnomonic(lon_center_degrees, lat_center_degrees)
          x, y = lonlat2xy(lon_degrees, lat_degrees)
    '''

    def __init__(self, lon_degrees, lat_degrees):
        lat0 = radians(lat_degrees)
        self.sin_lat0  = sin(lat0)
        self.cos_lat0  = cos(lat0)
        self.lon0 = radians(lon_degrees)

    def lonlat2xy(self, lon_degrees, lat_degrees):
        lat = radians(lat_degrees)
        lon = radians(lon_degrees)
        cos_dlon = cos(lon - self.lon0)
        sin_dlon = sin(lon - self.lon0)
        sin_lat  = sin(lat)
        cos_lat  = cos(lat)
        cos_c = radians( self.sin_lat0 * sin_lat + self.cos_lat0 * cos_lat * cos_dlon )
        x = (cos_lat * sin_dlon) / cos_c
        y = (self.sin_lat0 * cos_lat * cos_dlon - self.cos_lat0 * sin_lat) / cos_c
        return (x, y)

#====================================================================
def read_one_IN_file(filename, all_samples):
#====================================================================
    ''' Read one IN (input data) file, appending Sample objects to all_samples.
        This method **must** be overiden for project-specific "IN" files.
    '''
    print 'PROGRAMMING ERROR: Must override read_one_IN_file method'
    exit(-1)

#====================================================================
def read_all_IN_files():
#====================================================================
    ''' Read all of the IN (input data) files.
        (by calling the "read_one_IN_file" method)
    '''
    all_samples = []
    prev_nsamps = 0
    file_list = sorted( glob.glob(INFILE_PREFIX + '*.fits') )
    for filename in file_list:
        read_one_IN_file(filename, all_samples)
        total_nsamps = len(all_samples)
        this_nsamps = total_nsamps - prev_nsamps
        prev_nsamps = total_nsamps
        print_completed('reading '+ str(this_nsamps) +' samples from ' + filename)
    print 'A total of', str(len(all_samples)), 'samples were read in.'
    return all_samples


#====================================================================
# MAIN PROGRAM
#====================================================================

def main():

    get_paramaters(sys.argv)

    all_detectors = read_all_DRF_files()

    all_samples = read_all_IN_files()

    all_footprints = create_all_footprints(all_samples, all_detectors)

    # Create initial image
    wgt_array = calc_wgt_array(all_footprints)
    if STARTING_IMAGE is None:
        flux_wgt_array = calc_flux_wgt_array(all_footprints)
        flux_array = flux_wgt_array / wgt_array
    else: 
        flux_array = read_start_image()
        print 'reset ITER_START = ', ITER_START
    if ITER_START in ITER_LIST: write_FITS_image(flux_array, 'flux', ITER_START)

    # Main Loop --- compute correction array and then apply it (iterate)
    for iter in range(ITER_START+1, ITER_MAX+1):
        corr_wgt_array = calc_corr_wgt_array(all_footprints, flux_array, iter)
        flux_array *= corr_wgt_array / wgt_array
        if iter in ITER_LIST: write_FITS_image(flux_array, 'flux', iter)

    # Write FITS files
    write_FITS_image( wgt_array, 'cov')

    # ----------  show histogram (as a quick check)
    histo = np.histogram(flux_array, bins=50, range=(0.0, 50.0) )
    print histo[0]
    print histo[1]

