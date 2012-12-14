# HIRES" Python program to produce high resolution images
# Creates image using "maximum correlation method" (smoothest image that is
# consistent with the input data and the detector response functions)
# ==> Must override the "read_one_IN_file" method for reading input data files

PROGRAM = 'HIRES'
VERSION = 'v1_07'

#====================================================================
# imports
#====================================================================
import sys, os, exceptions, glob, time, logging
import numpy as np
import pyfits
from math import radians, degrees, cos, sin, atan2, sqrt, pow

#====================================================================
# log message to print and/or file, and maybe exit
#====================================================================
LOG_FILE_LEVEL  = 1 # put message in log file only if at least this level
LOG_PRINT_LEVEL = 1 # print message only if at least this level
LOG_QUIT_LEVEL  = 6 # only quit if message is at least this level
LOG_FILE_file = open('logfile.log', 'w')
def log(level, message, *values):
    '''
    print/log a message
    level:
       1 = extra message
       2 = standard message
       3 = time tagged message (e.g. start, completed, end)
       4 = warning (something odd, but not really an error)
       5 = error (something wrong, but may wish to continue)
       6 = fatal (impossible to continue)
    '''
    if len(values)>0: message = message % values
    if   level ==1: prelude = '       ... '
    elif level ==2: prelude = '       '
    elif level ==3: prelude = '%5.2f' % time.clock() + ' '
    elif level ==4: prelude = '** WARNING: '
    elif level ==5: prelude = '**** ERROR: '
    elif level ==6: prelude = '**** FATAL: '
    if level >= LOG_FILE_LEVEL: LOG_FILE_file.write(prelude+message +'\n')
    if level >= LOG_PRINT_LEVEL:
        if level >= 4: print ' '
        print prelude+message
        if level >= 4: print ' '
    if level >= LOG_QUIT_LEVEL:
        log(3, 'End PROCESSING *** EARLY TERMINATION ***')
        exit(-1)

#====================================================================
class Sample:
#====================================================================
    ''' one instance for each sample in the IN files '''
    flux_n_reset = 0
    flux_n_nonpositive = 0
    sum_noise = 0.0
    def __init__(self, x, y, flux, detector_id, angle, noise=1.0):
        self.x           = x # delta from center (degrees)
        self.y           = y # delta from center (degrees)
        if flux < FLUX_MIN:
            flux = FLUX_MIN
            Sample.flux_n_reset += 1
        if flux <= 0.0: Sample.flux_n_nonpositive += 1
        self.flux        = flux # signal in units that will go into image
        self.detector_id = detector_id # detector ID (string)
        self.angle       = angle # orientation degrees east of north
        self.noise       = noise # same units as flux
        Sample.sum_noise +=  noise

#====================================================================
class Detector:
#====================================================================
    ''' one instance for each detector, generally one DRF file '''
    all = {}
    def __init__(self, id, radius_degrees, response_function, drf_set_id):
        self.id   = id # detector ID (integer)
        self.radius_degrees = radius_degrees
        self.response_function = response_function
        self.drf_set_id = drf_set_id
        Detector.all[id] = self

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
        self.response_array = Footprint._get_response_array(detector_id, iFrac,jFrac, angle)
        bounds = Footprint._compute_footprint_bounds(self.response_array, iWhole, jWhole)

        Footprint.i0min = min(Footprint.i0min, bounds[0])
        Footprint.i1max = max(Footprint.i1max, bounds[1])
        Footprint.j0min = min(Footprint.j0min, bounds[2])
        Footprint.j1max = max(Footprint.j1max, bounds[3])

        self.footprint_bounds = bounds

    #------------------------------------------------------------------------
    @staticmethod
    def _get_response_array(detector_id, iFloat, jFloat, angle):
        ''' see if appropriate footprint array has already been generated
            If so, return it, otherwise call the method to generate it
        '''
        # compose key (to use previously computed footprint array)
        angle_round = round( angle / ANGLE_TOLERANCE )
        angleID = int( angle_round )
        angle = angle_round * ANGLE_TOLERANCE
        iFrac = iFloat % 1.0 # range:  0.0 <= frac < 1.0
        jFrac = jFloat % 1.0
        if FOOTPRINTS_PER_PIX <= 10:
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
                footArr =Footprint._generate_response_array(detector_id, iFracMod,jFracMod,angle)
                Footprint.all_response_arrays[footArr_key] = footArr

                print "DEBUG: angle:", angle, "full key:", footArr_key

        else:
            footArr = Footprint._generate_response_array(detector_id, iFrac, jFrac, angle)

        return footArr

    #------------------------------------------------------------------------
    @staticmethod
    def _generate_response_array(detector_id, iFrac, jFrac, angle):
        ''' generate full footprint array for given detector at
            given angle and centered at given x,y offset from pixel center
         detector_id
         angle (degrees) 
         iFrac, jFrac offset (frac of pixel) from center of pixel -.5 <= Frac < 0.5
        '''
        detector = Detector.all[detector_id]
        radius_degrees = detector.radius_degrees
        radius_pix = int( radius_degrees / DEG_PER_PIX )
        response_func = detector.response_function
        duFrac = iFrac * DEG_PER_PIX
        dvFrac = jFrac * DEG_PER_PIX
        duArray , dvArray =  \
         Footprint._generate_response_array_coords(radius_pix, angle, duFrac, dvFrac)
        response_array =  \
         Footprint._fill_in_response_array(duArray, dvArray, response_func)
        return response_array

    #------------------------------------------------------------------------
    @staticmethod
    def _generate_response_array_coords(radius_pix, angle, duFrac, dvFrac):
        ''' generate du,dv coordinate arrays of footprint
            returns: duArray, dvArray two np float arrays with coords in response func
        '''
        u0 = -radius_pix; u1 = radius_pix; v0 = -radius_pix; v1 = radius_pix
        ijNPIX = radius_pix + radius_pix + 1
        radius_degrees = radius_pix * DEG_PER_PIX # (now rounded to nearest pixel)
        cos_angle = cos(radians(angle))
        sin_angle = sin(radians(angle))

        # create the response array
        duArray = np.empty( (ijNPIX, ijNPIX), dtype=np.float32 )
        dvArray = np.empty( (ijNPIX, ijNPIX), dtype=np.float32 )
        dx = -radius_degrees
        for i in range(ijNPIX):
            dy = -radius_degrees
            for j in range(ijNPIX):
                du =  dx * cos_angle - dy * sin_angle - duFrac
                dv =  dy * cos_angle + dx * sin_angle - dvFrac
                duArray[i, j] = du
                dvArray[i, j] = dv
                dy += DEG_PER_PIX
            dx += DEG_PER_PIX
        return (duArray , dvArray)
       

    #------------------------------------------------------------------------
    @staticmethod
    def _fill_in_response_array(duArray, dvArray, response_func):
        ''' generate footprint response array, given u & v coordinate arrays
            input: duArray, dvArray: two np float arrays with coords in response func
            returns: response_array, an np array of floats, same dimensions as input arrays
               normalized to a sum of 1.0
        '''
        iNPIX, jNPIX = duArray.shape # iNpix should equal jNPIX
        response_array = np.empty( (iNPIX, jNPIX), dtype=np.float32 )
        for i in range(iNPIX):
            for j in range(jNPIX):
                response_here = response_func(duArray[i,j], dvArray[i,j])
                response_array[i, j] = response_here
        # normalize to sum of 1.0 and return
        sum = response_array.sum()
        normalized_response_array = response_array / sum # normalize
        return normalized_response_array

    #------------------------------------------------------------------------
    @staticmethod
    def _compute_footprint_bounds(full_reaponse_array, iWhole, jWhole):
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
FLUX_MIN = -1.0e200
FLUX_UNITS = '??'


# *decrease* du,dv to simulate *increased* axis length
# drf_fwhm = 4.23 # arcmin
# delta_axis = 2.0 # move between samples
# delta_axis = 2.80 # adjusted a bit for start/stop

# mean FWHM = 4.23
# kludge_factor = sqrt( (drf_fwhm + delta_axis) / drf_fwhm )
# KLUDGE_DU =      kludge_factor # cross-scan: inrease du to effectvely deccrease axis size
# KLUDGE_DV = 1.0 /kludge_factor # in-scan: decrease dv to effectvely increase axis size

# print "DEBUG: *** Using kludged du, dv (elongated gaussian) du*", KLUDGE_DU
# print "DEBUG: *** Using kludged du, dv (elongated gaussian) dv*", KLUDGE_DV
# print "DEBUG: *** Using kludged du, dv (elongated gaussian)"
# print "DEBUG: *** Using kludged du, dv (elongated gaussian)"

#====================================================================
def get_paramaters(args):
#====================================================================
    ''' Process user parameters 
        Everything put into global variables (ALL CAPS).
        Global variables should generally only be set in this method.
    '''
    global INFILE_PREFIX # directory and prefix of files containing samples
    global OUTFILE_PREFIX # directory and prefix for output filenames
    global OUTFILE_TYPES # 'flux' 'cov' beam' cfv padded
    global DRF_PREFIX # directory and prefix of files containing response func images
    global NPIXi, NPIXj # image size (pixels)
    global DEG_PER_PIX # image scale (CDELT) degrees per pixel
    global CRVAL1, CRVAL2 # coordinates of center image --- lon, lat in degrees
    global CTYPE1, CTYPE2 # coordinates type e.g. RA---TAN, DEC--TAN
    global ITER_MAX # max iteration nmber to compute
    global ITER_LIST # which iterations to write outupt files for
    global BOOST_MAX_ITER # max iteration do accalerated correction
    global BOOST_TYPE # type of accelerated correction to do
    global BOOST_FUNC # function to use for accelerated correction
    global STARTING_IMAGE # user supplied starting image or 'flat'
    global ANGLE_TOLERANCE # footprints re-used if within this tollerance (degrees)
    global FOOTPRINTS_PER_PIX # generate footprints at 1/N pixel accuracy
    global BEAM_SPIKE_N # NxN spikes in initial beam images
    global BEAM_SPIKE_HEIGHT # height of initial spikes for beam images
    global BEAM_STARTING_IMAGE # user supplied BEAM starting image
    global FLUX_MIN # lowest value allowed for flux in putut data (reset if less)
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
            if len(words)<2: log(6, 'parameter: %s has no value in file: %s', par, filename)
            name = words[0]
            val = words[1]
            if   name == 'SIZE_NPIX':
                NPIXi = int(val)
                if len(words)>=3: NPIXj = int(words[2])
                else: NPIXj = NPIXi
            elif name == 'ARCSEC_PER_PIX': arcsec_per_pix = float(val)
            elif name == 'CRVAL1': CRVAL1 = float(val)
            elif name == 'CRVAL2': CRVAL2 = float(val)
            elif name == 'CTYPE1': CTYPE1 = val
            elif name == 'CTYPE2': CTYPE2 = val
            elif name == 'OUTFILE_TYPES': OUTFILE_TYPES = ( words[1:] )
            elif name == 'ANGLE_TOLERANCE': ANGLE_TOLERANCE = float(val)
            elif name == 'FOOTPRINTS_PER_PIX': FOOTPRINTS_PER_PIX = int(val)
            elif name == 'DRF_PREFIX': DRF_PREFIX = val
            elif name == 'STARTING_IMAGE': STARTING_IMAGE = val
            elif name == 'BEAM_SPIKE_N': BEAM_SPIKE_N = int(val)
            elif name == 'BEAM_SPIKE_HEIGHT': BEAM_SPIKE_HEIGHT = float(val)
            elif name == 'BEAM_STARTING_IMAGE': BEAM_STARTING_IMAGE = val
            elif name == 'FLUX_UNITS': FLUX_UNITS = val
            elif name == 'FLUX_MIN': FLUX_MIN = float(val)
            elif name == 'ITER_MAX': ITER_MAX = int(val)
            elif name == 'ITER_LIST':
                for w in words[1:]:
                  iter_str =  w.split(',')
                  for n_str in iter_str:
                     if n_str != '': ITER_LIST[int(n_str)] = 1
            elif name == 'BOOST_CORRECTION':
                BOOST_MAX_ITER = int(val)
                if len(words) >= 3: BOOST_TYPE = words[2]
                if   BOOST_TYPE == 'TIMES_2': BOOST_FUNC = lambda x: x+x-1.0
                elif BOOST_TYPE == 'TIMES_3': BOOST_FUNC = lambda x: x+x+x-2.0
                elif BOOST_TYPE == 'SQUARED': BOOST_FUNC = lambda x: x*x
                elif BOOST_TYPE == 'EXP_2.5': BOOST_FUNC = lambda x: pow(x, 2.5)
                elif BOOST_TYPE == 'CUBED':   BOOST_FUNC = lambda x: x*x*x
                else: log(6, 'Unknown BOOST type: %s', BOOST_TYPE)
            elif name == 'KWD':
                if len(words) >=3: comment = ' '.join(words[3:])
                else: comment = None
                kwd = val
                val = numerify(words[2])
                need_warn = True
                if   kwd == 'CRVAL1': CRVAL1 = val;
                elif kwd == 'CRVAL2': CRVAL2 = val;
                elif kwd == 'CTYPE1': CTYPE1 = val;
                elif kwd == 'CTYPE2': CTYPE2 = val;
                else:
                    FITS_KEYWORDS.append( ( kwd, val, comment ) )
                    need_warn = False
                if need_warn: log(4, 'Using KWD before %s parameter is deprecated', kwd)
            elif name == '': pass
            else:
                print 'ERROR: Unknown parameter name:', name, 'in file:', filename
                exit(-2)

    # check for errors
    for t in OUTFILE_TYPES:
        if t not in ('flux', 'cov', 'beam', 'cfv', 'padded' ):
            log(6,'Illegal OUTFILE_TYPE: %s', t)

    # compute globals from pars
    DEG_PER_PIX = arcsec_per_pix / 3600.0
    ITER_LIST[ITER_MAX] = 1 # make sure to output files for final iteration 

#====================================================================
def print_paramaters():
#====================================================================
    log(2, '\nInput data file options:')
    log(2, '  INFILE_PREFIX %s', INFILE_PREFIX)
    log(2, '  FLUX_MIN %.2e', FLUX_MIN)
    log(2, '  STARTING_IMAGE %s', STARTING_IMAGE)

    log(2, '\nDRF (detector response files) to use:')
    log(2, '  DRF_PREFIX %s', DRF_PREFIX)

    log(2, '\nOutput image geometry:')
    log(2, '  NPIX %d %d', NPIXi, NPIXj)
    log(2, '  DEG_PER_PIX %.6f', DEG_PER_PIX)
    log(2, '  CRVAL1 %.5f', CRVAL1)
    log(2, '  CRVAL2 %.5f', CRVAL2)
    log(2, '  CTYPE1 %s', CTYPE1)
    log(2, '  CTYPE2 %s', CTYPE2)

    log(2, '\nOutput file options:')
    log(2, '  OUTFILE_PREFIX %s', OUTFILE_PREFIX)
    log(2, '  OUTFILE_TYPES %s', OUTFILE_TYPES)
    log(2, '  ITER_MAX %d', ITER_MAX)
    log(2, '  ITER_LIST %s', str(sorted(ITER_LIST.keys())) )
    log(2, '  FLUX_UNITS %s', FLUX_UNITS)

    if 'beam' in OUTFILE_TYPES:
        log(2, '\nBeam image file options:')
        log(2, '  BEAM_SPIKE_N %d', BEAM_SPIKE_N)
        log(2, '  BEAM_SPIKE_HEIGHT %f', BEAM_SPIKE_HEIGHT)
        log(2, '  BEAM_STARTING_IMAGE %s', BEAM_STARTING_IMAGE)

    log(2, '\nAdditional FITS keywords:')
    for k in FITS_KEYWORDS: log(2, "  KWD %s", str(k) )

    log(2, '\nAccelerated correction option:')
    if BOOST_MAX_ITER > 0: log(2, '  BOOST %s for ITER 2 to %d', BOOST_TYPE, BOOST_MAX_ITER)
    else: log(2, '  BOOST (none)')

    log(2, '\nFootprint accuracy options:')
    log(2, '  ANGLE_TOLERANCE %.2f', ANGLE_TOLERANCE)
    log(2, '  FOOTPRINTS_PER_PIX %d', FOOTPRINTS_PER_PIX)

    print " "
    log(3, 'Start PROCESSING  (' + PROGRAM  +' '+ VERSION +')')

#====================================================================
def read_all_DRF_files():
#====================================================================
    ''' Read all of the RDF (detector response function) files.
        (by calling the "read_one_DRF_file" method)
    '''
    global DRF_SET_ID
    all_detectors = {} # dictionay accessed by detector_id
    file_list = glob.glob(DRF_PREFIX + '*')
    if len(file_list)==0: log(6, 'No DRF files for: ' + DRF_PREFIX) # QUIT
    for file in file_list:
        id, detector = read_one_DRF_file(file)
        all_detectors[id] = detector
        DRF_SET_ID = detector.drf_set_id
    log(3, 'DRF file reading complete')
    log(2, '%d DRF files read', len(all_detectors) )
    log(2, 'DRF file set ID = %s', DRF_SET_ID)
    return all_detectors

#====================================================================
def read_all_IN_files():
#====================================================================
    ''' Read all of the IN (input data) files.
        (by calling the "read_one_IN_file" method)
    '''
    all_samples = []
    prev_nsamps = 0
    file_list = sorted( glob.glob(INFILE_PREFIX + '*.fits') )
    if len(file_list)==0: log(6, 'No input files for: ' + INFILE_PREFIX) # QUIT
    for filename in file_list:
        read_one_IN_file(filename, all_samples)
        total_nsamps = len(all_samples)
        this_nsamps = total_nsamps - prev_nsamps
        prev_nsamps = total_nsamps
        log(1,'read '+ str(this_nsamps) +' data samples from ' + filename)
    log(3, 'Input data sample reading complete')
    log(2, str(len(all_samples)) + ' input data samples read')
    return all_samples

#====================================================================
def create_all_footprints(setOfSamples, setOfDetectors):
#====================================================================
    ''' Create footprints for all samples in setOfSamples
        Use detector (DRF definitions) info given in seOfDetectors
    '''
    # do messages about samples read 
    if Sample.flux_n_reset != 0:
        log(2, 'Reset %d flux values in input data (to %e)', Sample.flux_n_reset, FLUX_MIN)
    if Sample.flux_n_nonpositive != 0: 
        log(4, 'Using %d negative fluxes in input data samples', Sample.flux_n_nonpositive)

    all_footprints = [] 
    mean_noise = Sample.sum_noise / float(len(setOfSamples))
    log(1, 'mean noise = %f', mean_noise)
    n_rejected = 0
    half_i = (NPIXi-1)/2
    half_j = (NPIXj-1)/2

    # DEBUG only needed for computing angle
    prev_x = None
    angle_calc = None

    for samp in setOfSamples:

        if prev_x is not None:
            delta_x = samp.x - prev_x
            if abs(delta_x) >= 0.05: # can't compute angle if coordinates jump
                angle_calc = None
            else:
                delta_y = samp.y - prev_y
                # ORIG: # angle_calc = degrees( atan2(delta_y, delta_x) )
                angle_calc = degrees( atan2(-delta_y, -delta_x) )
                distance_degrees = sqrt(delta_y*delta_y + delta_x*delta_x)
                # slope = delta_y / delta_x
                # print "DEBUG: dist=", "%.5f" % distance_degrees
                # print "DEBUG: slope=", "%.5f %.5f %.5f %.5f" % \
                #                     (angle_calc, slope, delta_x, delta_y)
        prev_x = samp.x
        prev_y = samp.y

        iFloat =  samp.x / DEG_PER_PIX
        jFloat =  samp.y / DEG_PER_PIX
        i = int(round( iFloat ))
        j = int(round( jFloat ))

        # Create Footprint if inside image area 
        # PREVIOUS:
        # if i>=-half_i and i<=half_i and j>=-half_j and j<=half_j:
        if i>=-half_i and i<=half_i and j>=-half_j and j<=half_j and angle_calc is not None:
            relative_noise = samp.noise / mean_noise
            wgt = 1.0 / (relative_noise * relative_noise)
            footprint = Footprint(samp.detector_id, iFloat,jFloat,angle_calc,samp.flux,wgt)
            # PREVIOUS:
            #footprint = Footprint(samp.detector_id, iFloat,jFloat,samp.angle,samp.flux,wgt)
            all_footprints.append(footprint)
        else:
            n_rejected += 1

    log(3, 'Footprint creation complete')
    if len(all_footprints)==0: log(6, 'All data samples rejected (probably out of image)' )
    log(2, '%d data samples accepted (footprints created)', len(all_footprints) )
    log(2, '%d data samples rejected', n_rejected)
    log(1, '(%d full response arrays created)',  len(Footprint.all_response_arrays) ) 
    pct_accepted = float(len(all_footprints)) / float(len(setOfSamples))
    if pct_accepted<0.50: log(4, 'More than 50% of data samples rejected' )

    return all_footprints

#====================================================================
def adjust_for_padding(setOfFootprints):
#====================================================================
    ''' 1. Calculate padding adjustments to i,j bounds
        2. Apply to Footprint.footprint_bounds
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
        i0, i1, j0, j1 = foot.footprint_bounds
        new_bounds = (i0+delta_i, i1+delta_i, j0+delta_j, j1+delta_j)
        foot.footprint_bounds = new_bounds

#====================================================================
def calc_wgt_image(footprints):
#====================================================================
    ''' Create an image of all footprints *1.0 '''
    wgt_image = np.zeros( (NPIXiPadded, NPIXjPadded), dtype=np.float32)
    wgt_image += 0.0000001 # to avoid NaNs in results
    for foot in footprints:
        i0, i1, j0, j1 = foot.footprint_bounds
        wgt_image[i0:i1, j0:j1] += foot.response_array * foot.weight
    log(3, 'Weight array computed')
    return wgt_image

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
            log(6, "STARTING_IMAGE %s has incorrect dimensions\n" + \
             '   Must be: %d %d  (requires padding)', filename, NPIXiPadded, NPIXjPadded)
        if (kwd['CRVAL1'] != CRVAL1) or (kwd['CRVAL2'] != CRVAL2):
            log(5, "STARTING_IMAGE %s has inconsistent CRVAL value(s)\n" + \
             '   Should be: %.5f %.5f', filename, CRVAL1, CRVAL2)
        if (kwd['BUNIT'] != FLUX_UNITS):
            log(5, "STARTING_IMAGE %s has inconsistent BUNIT\n" + \
             "   Should be: '%s'", filename, FLUX_UNITS)
        if kwd.has_key('ITERNUM'):
            iter_start = kwd['ITERNUM']
            log(2, 'ITERNUM in %s is %d', filename, iter_start)
        else: log(4, 'No ITERNUM keyword in START_IMAGE %s', filename)
        image = hduList[0].data
    return (image, iter_start)

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
        i0, i1, j0, j1 = foot.footprint_bounds
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
    log(3, 'Correction array computed, iter ' + str(iter_num) )
    if boost_mode is not None: log(2, boost_mode)
    return ( corr_wgt_image , corr_sq_wgt_image )

#====================================================================
def set_fluxes_to_sim_values(footprints, sim_image):
#====================================================================
    ''' Calculate simulated fluxes; replace foot.flux value '''
    for foot in footprints:
        i0, i1, j0, j1 = foot.footprint_bounds
        integration = foot.response_array * sim_image[i0:i1, j0:j1]
        sim_flux = integration.sum()
        foot.flux = sim_flux
    log(3, 'Fluxes have been reset to simulated values')

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
def get_FITS_keyword(keyword):
#====================================================================
    ''' Get a FITS keyword from FITS_KEYWORDS.  Return None if not there. '''
    if keyword == 'CRVAL1': return CRVAL1
    if keyword == 'CRVAL2': return CRVAL2
    if keyword == 'CTYPE1': return CTYPE1
    if keyword == 'CTYPE2': return CTYPE2
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
    prihdr.update('CRVAL1', CRVAL1)
    prihdr.update('CRVAL2', CRVAL2)
    prihdr.update('CTYPE1', CTYPE1)
    prihdr.update('CTYPE2', CTYPE2)
    cdelt_rounded = float('%.7f' % DEG_PER_PIX)
    prihdr.update('CDELT1', -cdelt_rounded, comment='left(+) to right(-)')
    prihdr.update('CDELT2',  cdelt_rounded, comment='bottom(-) to top(+)')
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
    prihdr.update('INFILES' , in_pref  , comment='Name of input data files')  
    if DRF_PREFIX is not None:
        drf_pref = DRF_PREFIX.split('/')[-1] +'*'
        prihdr.update('DRF_IN'  , drf_pref  , comment='Name of Detector Response Files')  
    prihdr.update('DRF_ID', DRF_SET_ID  , comment='ID of Detector Response Files')  
    dt = time.strftime("%Y/%m/%d %H:%M:%S")
    prihdr.update('DATE'     , dt  , comment='when this file was created')  
    prihdr.update('CREATED', PROGRAM+' '+VERSION, comment='software version that created this file')

  # prihdr.update('KLUGE_DU', KLUDGE_DU, comment='DRF cross-scan axis size divided by this')
  # prihdr.update('KLUGE_DV', KLUDGE_DV, comment='DRF in-scan axis size divided by this')

    # write file
    hdulist = pyfits.HDUList([hdu])
    if os.path.isfile(filename): os.remove(filename)
    hdulist.writeto(filename)
    log(3, 'Output file written: %s', filename)

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
    cdelt1 = kwd['CDELT1']
    cdelt2 = kwd['CDELT2']
    if naxis1 != naxis2: log(5, 'in ' + filename + ' NAXIS1 must equal NAXIS2')
    if naxis1 %2==0 : log(5, 'in ' + filename + ' NAXIS1 must be odd')
    if abs(cdelt1) != abs(cdelt2):
        log(5, 'in ' + filename + ' CDELT1 and CDELT2 must be same size')
    if cdelt1 >= 0: log(4, 'in ' + filename + ' CDELT1 must be negative')
    if cdelt2 <= 0: log(4, 'in ' + filename + ' CDELT2 must be positive')
    deg_per_pix = cdelt2

    radius_pix = naxis1 / 2
    radius_degrees = radius_pix * deg_per_pix

    if ('DRF_ID') in kwd: drf_set_id = kwd['DRF_ID']
    else:                 drf_set_id = '(unknown)'

    drf_array = hduList[0].data

    def duv2response_function(du, dv):
        ''' interpolate in detector response array
         du, dv in degrees relative to DRF center
         '''
        ## DEBUG KLUDGE to transform circular Gaussian into elliptical one
        # multiplying by factor effectively makes that axis *shorter* by that factor
        # du *= KLUDGE_DU
        # dv *= KLUDGE_DV

        iFloat = (du + radius_degrees) / deg_per_pix
        jFloat = (dv + radius_degrees) / deg_per_pix
        # cheap "interpolation" -- just takes nearest one
        iInt = int(iFloat+0.5)
        jInt = int(jFloat+0.5)
        if iInt<0 or iInt>=naxis1 or jInt<0 or jInt>=naxis2 :
            response = 0.0
        else:
            response = drf_array[iInt, jInt]
        return response

    log(3, 'detector: %d  file=%s radius= %d pixels = %f arcmin', \
      id, filename.split('/')[-1], radius_pix, radius_degrees*60 ) 
    detector = Detector(id, radius_degrees, duv2response_function, drf_set_id)
    return ( id, detector )

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
# MAIN PROGRAM
#====================================================================

def main():

    #-----------------------------------------------------------
    # Initialize
    #-----------------------------------------------------------
    get_paramaters(sys.argv)
    print_paramaters()
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

    log(3, 'End PROCESSING')
