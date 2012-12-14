# "HIRES" Python program to produce high resolution images
# Creates image using "maximum correlation method" (smoothest image that is
# consistent with the input data and the detector response functions)
# ==> Must override the "read_one_IN_file" method

PROGRAM = 'HIRES'
VERSION = 'v1_05'

#====================================================================
# imports
#====================================================================
import sys, os, exceptions, glob, time, logging
import numpy as np
import pyfits
from math import radians, cos, sin, pow

#====================================================================
# initialize logging
#====================================================================
logging.basicConfig(
  filename='log_hires.log',
  filemode='w',
  format='%(asctime)s %(message)s',
  datefmt='%Y/%m/%d %H:%M:%S',
  level=logging.INFO  )

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
    def __init__(self, id, radius_degrees, response_function):
        self.id   = id # detector ID (integer)
        self.radius_degrees = radius_degrees
        self.response_function = response_function
        Detector.all[id] = self

#====================================================================
class Sample:
#====================================================================
    ''' one instance for each sample in the IN files '''
    sum_noise = 0.0
    def __init__(self, x, y, flux, detector_id, angle, noise=1.0):
        self.x           = x # delta from center (degrees)
        self.y           = y # delta from center (degrees)
        self.flux        = flux # signal in units that will go into image
        self.detector_id = detector_id # detector ID (string)
        self.angle       = angle # orientation degrees east of north
        self.noise       = noise # same units as flux
        Sample.sum_noise +=  noise

#====================================================================
class Footprint:
#====================================================================
    ''' A Footprint is a mapping of the detector response function
    on the image.  There is one footprint for each sample (unless out of range).
    The Footprint is basically an array of responses for each image pixel
    near the sample.
      Note: i,j "bounds" are ~Python style bounds i.e . i0 = min(i); i1 = max(i) + 1
    '''

    all_response_arrays = {} # indexed by (detID, iFracID, jFracID, angleID)
    i0min =  9999; i1max = -9999; j0min =  9999; j1max = -9999

    def __init__(self, detector_id, iFloat, jFloat, angle, flux, weight):
        self.flux = flux
        self.weight = weight
        iWhole, iFrac = divmod(iFloat+0.5, 1.0)
        jWhole, jFrac = divmod(jFloat+0.5, 1.0)
        self.response_array = self._get_response_array(detector_id, iFrac, jFrac, angle)
        image_bounds = self._compute_image_bounds(self.response_array, iWhole, jWhole)

        Footprint.i0min = min(Footprint.i0min, image_bounds[0])
        Footprint.i1max = max(Footprint.i1max, image_bounds[1])
        Footprint.j0min = min(Footprint.j0min, image_bounds[2])
        Footprint.j1max = max(Footprint.j1max, image_bounds[3])

        self.image_bounds = image_bounds

    #------------------------------------------------------------------------
    def _get_response_array(self, detector_id, iFloat, jFloat, angle):
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
                footArr =self._generate_response_array(detector_id, iFracMod,jFracMod,angle)
                Footprint.all_response_arrays[footArr_key] = footArr
        else:
            footArr = self._generate_response_array(detector_id, iFrac, jFrac, angle)

        return footArr

    #------------------------------------------------------------------------
    def _generate_response_array(self, detector_id, iFrac, jFrac, angle):
        ''' generate full footprint array for given detector at
            given angle and centered at given x,y offset from pixel center
         detector_id
         angle (degrees) 
         iFrac, jFrac offset (frac of pixel) from center of pixel -.5 <= Frac < 0.5
        '''
        detector = Detector.all[detector_id]
        radius_degrees = detector.radius_degrees
        radius_pix = int( radius_degrees / DEG_PER_PIX )
        radius_degrees = radius_pix * DEG_PER_PIX # (now rounded to nearest pixel)
        response_func = detector.response_function
        u0 = -radius_pix; u1 = radius_pix; v0 = -radius_pix; v1 = radius_pix
        ijNPIX = radius_pix + radius_pix + 1
        cos_angle = cos(radians(angle))
        sin_angle = sin(radians(angle))
        duFrac = iFrac * DEG_PER_PIX
        dvFrac = jFrac * DEG_PER_PIX

        # print "==================   "
        # print iFrac, jFrac, duFrac, dvFrac
        # print angle, sin_angle, cos_angle

        # create the response array
        foot_array = np.empty( (ijNPIX, ijNPIX), dtype=np.float32 )
        dx = -radius_degrees
        for i in range(ijNPIX):
            dy = -radius_degrees
            for j in range(ijNPIX):
                du =  dx * cos_angle + dy * sin_angle - duFrac
                dv = -dx * sin_angle + dy * cos_angle - dvFrac
                response_here = response_func(du, dv)
                foot_array[i, j] = response_here
                dy += DEG_PER_PIX
            dx += DEG_PER_PIX

        # normalize to sum of 1.0 and return
        foot_sum = foot_array.sum()
        normalized_foot_array = foot_array / foot_sum # normalize
        return normalized_foot_array

    #------------------------------------------------------------------------
    def _compute_image_bounds(self, full_reaponse_array, iWhole, jWhole):
        ''' Compute where this foorprint lands in the image.
            Could lie partially off the image but will be trimmed later.
        '''
        iSize, jSize = full_reaponse_array.shape
        i0 = int(iWhole) - iSize/2
        i1 = i0 + iSize
        j0 = int(jWhole) - jSize/2
        j1 = j0 + jSize
        return (i0, i1, j0, j1)


#====================================================================
def create_all_footprints(setOfSamples, setOfDetectors):
#====================================================================
    ''' Create footprints for all samples in setOfSamples
        Use detector (DRF definitions) info given in seOfDetectors
    '''
    all_footprints = [] ## 0=image_bounds, 1=array, 2=sum(array)  3=flux
    mean_noise = Sample.sum_noise / float(len(setOfSamples))
    print 'mean noise=', mean_noise
    n_rejected = 0
    half_i = (NPIXi-1)/2
    half_j = (NPIXj-1)/2

    for samp in setOfSamples:
        iFloat =  samp.x / DEG_PER_PIX
        jFloat =  samp.y / DEG_PER_PIX
        i = int(round( iFloat ))
        j = int(round( jFloat ))

        # Create Footprint if inside image area 
        if i>=-half_i and i<=half_i and j>=-half_j and j<=half_j:
            relative_noise = samp.noise / mean_noise
            wgt = 1.0 / (relative_noise * relative_noise)
            footprint = Footprint(samp.detector_id, iFloat,jFloat,samp.angle,samp.flux,wgt)
            all_footprints.append(footprint)
        else:
            n_rejected += 1

    print_completed( 'footprint creation: ' +str(len(all_footprints))+ ' Footprints' )
    print 'Rejected', n_rejected, 'Samples and created', \
     len(Footprint.all_response_arrays), 'full response arrays'

    return all_footprints

#====================================================================
def adjust_for_padding(setOfFootprints):
#====================================================================
    ''' 1. Calculate padding adjustments to i,j bounds
        2. Apply to Footprint.image_bounds
    '''

    global UNPAD_bounds # ij bounds to get real image from padded image
    global NPIXiPadded # size of padded image arrrays
    global NPIXjPadded

    # print 'i0min=', Footprint.i0min
    # print 'i1max=', Footprint.i1max
    # print 'j0min=', Footprint.j0min
    # print 'j1max=', Footprint.j1max

    # Compute size of padded image
    NPIXiPadded = Footprint.i1max - Footprint.i0min
    NPIXjPadded = Footprint.j1max - Footprint.j0min
    # print 'NPIXiPadded=', NPIXiPadded 
    # print 'NPIXjPadded=', NPIXjPadded 

    # compute offsets to go FROM range -N/2:+N/2 TO  range 0:N
    delta_i = -Footprint.i0min
    delta_j = -Footprint.j0min
    # print "DELTAs", delta_i, delta_j

    # Compute bounds of "real" image within padded image
    half_i = (NPIXi-1)/2
    half_j = (NPIXj-1)/2
    unpad_i0 = delta_i - half_i
    unpad_i1 = unpad_i0 + NPIXi
    unpad_j0 = delta_j - half_j
    unpad_j1 = unpad_j0 + NPIXj
    # print 'UNPADs', unpad_i0, unpad_i1, unpad_j0, unpad_j1

    UNPAD_bounds =  ( unpad_i0, unpad_i1, unpad_j0, unpad_j1 )

    for foot in setOfFootprints:
        i0, i1, j0, j1 = foot.image_bounds
        new_bounds = (i0+delta_i, i1+delta_i, j0+delta_j, j1+delta_j)
        foot.image_bounds = new_bounds

#====================================================================
def make_start_image(filename):
#====================================================================
    ''' starting image - use flat image if 'flat'
        otherwise read in starting image '''
    iter_start = 0 # will be set non-zero if in starting image selected
    if filename == 'flat': # make a flat image to start
        image = np.zeros( (NPIXiPadded, NPIXjPadded), dtype=np.float32 )
        background = 1.0
        image += background
    else: # read image from file
        hduList = pyfits.open(filename)
        kwd = hduList[0].header
        # check that dimensions are correct
        if (kwd['NAXIS1'] != NPIXiPadded) or (kwd['NAXIS2'] != NPIXjPadded):
            print "ERROR: STARTING_IMAGE has incorrect dimensions"
            print '   (must be padded) ...', NPIXiPadded, 'by', NPIXjPadded
            exit(-1)
        if kwd.has_key('ITERNUM'):
            iter_start = kwd['ITERNUM']
            print 'ITERNUM in ', filename, ' is', iter_start
        image = hduList[0].data
    return (image, iter_start)

#====================================================================
def calc_wgt_image(footprints):
#====================================================================
    ''' Create an image of all footprints *1.0 '''
    wgt_image = np.zeros( (NPIXiPadded, NPIXjPadded), dtype=np.float32)
    wgt_image += 0.0000001 # to avoid NaNs in results
    for foot in footprints:
        i0, i1, j0, j1 = foot.image_bounds
        wgt_image[i0:i1, j0:j1] += foot.response_array * foot.weight
    print_completed('calc_wgt_image')
    return wgt_image

#====================================================================
def calc_corr_wgt_image(footprints, flux_image, iter_num, do_cfv):
#====================================================================
    ''' Create a correction image: For each footprint:
         1. Calculate "correcton" 
         2. Add into correction image: response_array * correction
    '''
    corr_wgt_image = np.zeros( (NPIXiPadded, NPIXjPadded), dtype=np.float32)
    if do_cfv: corr_sq_wgt_image = np.zeros((NPIXiPadded, NPIXjPadded),dtype=np.float32)
    else: corr_sq_wgt_image = None
    boost_mode = None
    for foot in footprints:
        i0, i1, j0, j1 = foot.image_bounds
        integration = foot.response_array * flux_image[i0:i1, j0:j1]
        flux_prime = integration.sum()
        correction = (foot.flux / flux_prime) * foot.weight
        if iter_num != 1 and iter_num <= BOOST_MAX_ITER:
            correction = BOOST_FUNC(correction)
            boost_mode = "   (BOOSTED correction)"
        corr_wgt_image[i0:i1, j0:j1] += foot.response_array * correction
        if do_cfv:
            corr_sq = correction * correction
            corr_sq_wgt_image[i0:i1, j0:j1] += foot.response_array * corr_sq

    min_c = np.nanmin(corr_wgt_image)
    max_c = np.nanmax(corr_wgt_image)
    std_c = np.std(corr_wgt_image)
    stats = 'iter=%3d min=%.6f max=%.6f stdev=%.6f' % (iter_num, min_c, max_c, std_c)
    print_completed('calc_corr_wgt_image ' + stats)
    if boost_mode is not None: print boost_mode

    return ( corr_wgt_image , corr_sq_wgt_image )

#====================================================================
def set_fluxes_to_sim_values(footprints, sim_image):
#====================================================================
    ''' Calculate simulated fluxes; replace foot.flux value '''
    for foot in footprints:
        i0, i1, j0, j1 = foot.image_bounds
        integration = foot.response_array * sim_image[i0:i1, j0:j1]
        sim_flux = integration.sum()
        foot.flux = sim_flux
    print_completed('set_fluxes_to_sim_values ')

#====================================================================
def create_spike_image(n, height):
#====================================================================
    ''' Create initial spike image
        (values are all small>0 except for nXn spikes of given height '''
    spike_image = np.zeros( (NPIXiPadded, NPIXjPadded), dtype=np.float32)
    spike_image += 0.000001
    i0, i1, j0, j1 = UNPAD_bounds
    iDelta = NPIXi / n
    jDelta = NPIXj / n
    for i in range(i0+iDelta/2, i0+NPIXi, iDelta):
        for j in range(j0+jDelta/2, j0+NPIXj, jDelta):
            spike_image[i,j] = height
    return spike_image

#====================================================================
def set_FLUX_UNITS(units):
#====================================================================
    ''' Set FLUX_UNITS global variable, the BUNIT keyword in output files.
        (intended  for use in "read_one_IN_file") '''
    global FLUX_UNITS
    FLUX_UNITS = units

#====================================================================
# Initialize parameters
#====================================================================
FITS_KEYWORDS = [] 
ITER_LIST = {}
OUTFILE_TYPES = ('flux', )
STARTING_IMAGE = 'flat'
FOOTPRINTS_PER_PIX = 1
BEAM_STARTING_IMAGE = 'flat'
BEAM_SPIKE_N = 5
BEAM_SPIKE_HEIGHT = 10.0
BOOST_MAX_ITER = 0
FLUX_UNITS = '??'

#====================================================================
def get_paramaters(args):
#====================================================================
    ''' Process user parameters 
        Everything put into global variables (ALL CAPS).
        Global variables should generally only be set in this method.
    '''
    global INFILE_PREFIX # directory and prefix of files containing samples
    global OUTFILE_PREFIX # directory and prefix for output filenames
    global OUTFILE_TYPES # 'flux' 'cov' beam' etc.
    global DRF_PREFIX # directory and prefix of files containing response func images
    global DEG_PER_PIX # image scale (CDELT) degrees per pixel
    global NPIXi, NPIXj # image size (pixels)
    global ITER_MAX # max iteration nmber to compute
    global ITER_LIST # which iterations to write outupt files for
    global BOOST_MAX_ITER # max iteration do accalerated correction
    global BOOST_FUNC # type of accelerated correction to do
    global STARTING_IMAGE # user supplied starting image or 'flat'
    global ANGLE_TOLERANCE # footprints re-used if within this tollerance (degrees)
    global FOOTPRINTS_PER_PIX # generate footprints at 1/N pixel accuracy
    global BEAM_SPIKE_N # NxN spikes in initial beam images
    global BEAM_SPIKE_HEIGHT # height of initial spikes for beam images
    global BEAM_STARTING_IMAGE # user supplied BEAM starting image
    global FLUX_UNITS # for FITS keyword BUNIT in flux images
    global FITS_KEYWORDS # FITS keywords to be added to output files

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
                NPIXi = int(val)
                if len(words)>=3: NPIXj = int(words[2])
                else: NPIXj = NPIXi
            elif name == 'ARCSEC_PER_PIX': arcsec_per_pix = float(val)
            elif name == 'FLUX_UNITS': FLUX_UNITS = val
            elif name == 'OUTFILE_TYPES': OUTFILE_TYPES = ( words[1:] )
            elif name == 'ANGLE_TOLERANCE': ANGLE_TOLERANCE = float(val)
            elif name == 'FOOTPRINTS_PER_PIX': FOOTPRINTS_PER_PIX = int(val)
            elif name == 'DRF_PREFIX': DRF_PREFIX = val
            elif name == 'STARTING_IMAGE': STARTING_IMAGE = val
            elif name == 'BEAM_SPIKE_N': BEAM_SPIKE_N = int(val)
            elif name == 'BEAM_SPIKE_HEIGHT': BEAM_SPIKE_HEIGHT = float(val)
            elif name == 'BEAM_STARTING_IMAGE': BEAM_STARTING_IMAGE = val
            elif name == 'ITER_MAX': ITER_MAX = int(val)
            elif name == 'ITER_LIST':
                for w in words[1:]:
                  iter_str =  w.split(',')
                  for n_str in iter_str:
                     if n_str != '': ITER_LIST[int(n_str)] = 1
            elif name == 'BOOST_CORRECTION':
                BOOST_MAX_ITER = int(val)
                if len(words) >= 3: boost_type = words[2]
                if   boost_type == 'TIMES_2': BOOST_FUNC = lambda x: x+x-1.0
                elif boost_type == 'TIMES_3': BOOST_FUNC = lambda x: x+x+x-2.0
                elif boost_type == 'SQUARED': BOOST_FUNC = lambda x: x*x
                elif boost_type == 'EXP_2.5': BOOST_FUNC = lambda x: pow(x, 2.5)
                elif boost_type == 'CUBED':   BOOST_FUNC = lambda x: x*x*x
                else:
                    print 'ERROR: Unknown BOOST type:', boost_type
                    exit(-2)
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

    # compute globals from pars
    DEG_PER_PIX = arcsec_per_pix / 3600.0
    ITER_LIST[ITER_MAX] = 1

    # print all params
    print 'hires.py VERSION:', VERSION
    print 'INFILE_PREFIX', INFILE_PREFIX
    print 'OUTFILE_PREFIX', OUTFILE_PREFIX
    print 'OUTFILE_TYPES', OUTFILE_TYPES
    print 'DRF_PREFIX', DRF_PREFIX
    print 'NPIXi', NPIXi
    print 'NPIXj', NPIXj
    print 'ARCSEC_PER_PIX', arcsec_per_pix
    print 'FLUX_UNITS', FLUX_UNITS
    print 'ITER_MAX', ITER_MAX
    print 'ITER_LIST', sorted( ITER_LIST.keys() )
    if BOOST_MAX_ITER > 0: print 'BOOST', boost_type, 'for ITER 2 to', BOOST_MAX_ITER
    else: print 'BOOST (none)'
    print 'STARTING_IMAGE', STARTING_IMAGE
    print 'ANGLE_TOLERANCE', ANGLE_TOLERANCE
    print 'FOOTPRINTS_PER_PIX', FOOTPRINTS_PER_PIX
    print 'BEAM_SPIKE_N', BEAM_SPIKE_N
    print 'BEAM_SPIKE_HEIGHT', BEAM_SPIKE_HEIGHT
    print 'BEAM_STARTING_IMAGE', BEAM_STARTING_IMAGE
    for k in FITS_KEYWORDS: print "KWD " + str(k)

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
    try: return int(s)
    except exceptions.ValueError:
        try: return float(s)
        except exceptions.ValueError: return s

#====================================================================
def write_FITS_image(image, file_type, iter=None, trim_padding=True):
#====================================================================
    ''' Writes out one FITS image of specified type.
         image = array with pixel values
         file_type = 'flux', 'cov', 'cfv', etc.
         iter = iteration number (put in FITS keyword)
         trim_padding = Should padding be removed ?
    '''
    directory = os.path.dirname(OUTFILE_PREFIX)
    if not os.path.exists(directory): os.makedirs(directory)
    filename = OUTFILE_PREFIX + '_' + file_type
    if iter is not None: filename += '_' + str(iter)
    i0, i1, j0, j1 = UNPAD_bounds
    if trim_padding:
        output_image = image[i0:i1, j0:j1] # trim off padding area
    else:
        output_image = image # keep the padding
        filename += '_padded'
    filename += '.fits'

    hdu = pyfits.PrimaryHDU(output_image)
    prihdr = hdu.header

    # add keywords
    if (file_type=='flux'): prihdr.update('BUNIT', FLUX_UNITS)
    for kwd_tuple in FITS_KEYWORDS: # user specified kwds
        prihdr.update(kwd_tuple[0], kwd_tuple[1], comment=kwd_tuple[2])
    cdelt_rounded = float('%.7f' % DEG_PER_PIX)
    prihdr.update('CDELT1', -cdelt_rounded, comment='left to right')
    prihdr.update('CDELT2',  cdelt_rounded, comment='scan direction at top')
    if trim_padding:
        prihdr.update('CRPIX1', (NPIXi+1)/2 , comment='center pixel')
        prihdr.update('CRPIX2', (NPIXj+1)/2 , comment='center pixel')
    else:
        prihdr.update('CRPIX1', (NPIXi+1)/2 + i0 , comment='center pixel')
        prihdr.update('CRPIX2', (NPIXj+1)/2 + j0 , comment='center pixel')
        prihdr.update('PADDING', 'KEPT', comment='Image pixel padding has been kept')
    if   file_type=='flux': t_comment = 'HIRES flux image'
    elif file_type=='cov':  t_comment = 'HIRES coverage image'
    elif file_type=='cfv':  t_comment = 'HIRES correction factor variance image'
    else:                   t_comment = 'HIRES image'
    prihdr.update('FILETYPE', file_type, comment=t_comment)
    if iter is not None: prihdr.update('ITERNUM', iter, comment='HIRES iteration number')
    filename_only = filename.split('/')[-1] # strip off leading directory names
    prihdr.update('FILENAME', filename_only, comment='name of this file')  
    in_pref = INFILE_PREFIX.split('/')[-1] +'*'
    prihdr.update('INFILES' , in_pref  , comment='INput data files')  
    if DRF_PREFIX is not None:
        drf_pref = DRF_PREFIX.split('/')[-1] +'*'
        prihdr.update('DRF'     , drf_pref  , comment='Detector Response Files used')  
    dt = time.strftime("%Y/%m/%d %H:%M:%S")
    prihdr.update('DATE'     , dt  , comment='when this file was created')  
    prihdr.update('CREATED', PROGRAM+' '+VERSION, comment='software version that created this file')

    # write file
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
    npix = kwd['NAXIS1']

    # naxis2 = kwd['NAXIS2']
    # crpix1 = kwd['CRPIX1']
    # crpix2 = kwd['CRPIX2']
    print 'check CRPIX1,2 are right; and NAXIS1=NAXIS2 and odd'
    print 'and signs of CDELTs are right (WARN)'

    deg_per_pix = kwd['CDELT1']
    radius_pix = npix / 2
    radius_degrees = radius_pix * deg_per_pix

    print "radius=", radius_pix, radius_degrees

    drf_array = hduList[0].data

    def xy2response_function(x, y):
        ''' interpolate in detector response array
         x, y in degrees (relative to DRF center
         '''
        iFloat = (x + radius_degrees) / deg_per_pix
        jFloat = (y + radius_degrees) / deg_per_pix
        # cheap "interpolation" -- just takes nearest one
        iInt = int(iFloat+0.5)
        jInt = int(jFloat+0.5)
        if iInt<0 or iInt>=npix or jInt<0 or jInt>=npix :
            response = 0.0
        else:
            response = drf_array[iInt, jInt]

        return response

    print "Detector: ", id, 'file=', filename.split('/')[-1], \
            'radius:', radius_pix, 'pixels', radius_degrees*60, ' arcmin'
    detector = Detector(id, radius_degrees, xy2response_function)
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

    #-----------------------------------------------------------
    # Initialize
    #-----------------------------------------------------------
    get_paramaters(sys.argv)
    all_detectors = read_all_DRF_files()
    all_samples = read_all_IN_files()
    all_footprints = create_all_footprints(all_samples, all_detectors)
    adjust_for_padding(all_footprints)
    wgt_image = calc_wgt_image(all_footprints)
    if 'cov' in OUTFILE_TYPES: write_FITS_image( wgt_image, 'cov')

    #-----------------------------------------------------------
    # Create FLUX image(s)
    #-----------------------------------------------------------
    if 'flux' in OUTFILE_TYPES:
        flux_image, iter_start = make_start_image(STARTING_IMAGE)
        for iter in range(iter_start+1, ITER_MAX+1):
            do_cfv_image =  ('cfv' in OUTFILE_TYPES) and (iter in ITER_LIST)
            corr_wgt_image, corr_sq_wgt_image = \
              calc_corr_wgt_image(all_footprints, flux_image, iter, do_cfv_image)
            correction_image = corr_wgt_image / wgt_image
            flux_image *= correction_image
            # print 'Mean flux in image =', flux_image.mean()
            if iter in ITER_LIST: write_FITS_image(flux_image, 'flux', iter)
            if do_cfv_image: 
                corr_sq_image = (corr_sq_wgt_image / wgt_image) - \
                                (correction_image * correction_image)
                write_FITS_image(corr_sq_image, 'cfv', iter)
        if 'padded' in OUTFILE_TYPES:
            write_FITS_image(flux_image, 'flux', iter, trim_padding=False)

    #-----------------------------------------------------------
    # Create BEAM image(s)
    #-----------------------------------------------------------
    if 'beam' in OUTFILE_TYPES:
        spike_image = create_spike_image(BEAM_SPIKE_N, BEAM_SPIKE_HEIGHT)
        set_fluxes_to_sim_values(all_footprints, spike_image) # reset Sample.flux
        beam_image, iter_start = make_start_image(BEAM_STARTING_IMAGE)
        for iter in range(iter_start+1, ITER_MAX+1):
            corr_wgt_image, corr_sq_wgt_image = \
              calc_corr_wgt_image(all_footprints, beam_image, iter, False)
            beam_image *= corr_wgt_image / wgt_image
            if iter in ITER_LIST: write_FITS_image(beam_image, 'beam', iter)
        if 'padded' in OUTFILE_TYPES:
            write_FITS_image(beam_image, 'beam', iter, trim_padding=False)

