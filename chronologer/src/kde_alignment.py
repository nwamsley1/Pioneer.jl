'''
KDE-based retention time alignment functions
'''

import itertools
import numpy as np
from scipy.interpolate import interp1d

from src.peptide import modseq_to_seq


class cosine_gaussian:
    def __init__(self, mean=0.0, sd=1.0, prior=1.0):
        self.mean = mean
        self.sd = sd
        self.prior = prior
        self.A = np.pi / (8.0*sd)
        self.M = np.pi / (4.0*sd)
        self.min = mean - 2.0*sd
        self.max = mean + 2.0*sd
    def compute_pdf( self, x ):
        if x < self.min or x > self.max:
            return 0.0
        else:
            return self.A * np.cos( (self.mean - x) * self.M )
    def compute_cdf( self, x ):
        if x < self.min or x > self.max:
            return 0.0
        else:
            d = (x-self.mean) / self.sd
            return 1.0 / (1.0+np.exp(-0.07056*np.power(d,3)-1.5976*d))
    def return_sd( self ):
        return self.sd


def define_stamp( dist ):
    sd = dist.return_sd()
    stamp_radius = int( np.round( 2.0 * sd ) )
    stamp_width = 2*stamp_radius + 1
    stamp = np.zeros( (stamp_width, stamp_width) )
    for i,j in itertools.product( range(stamp_width), range(stamp_width) ):
        distance = np.sqrt( (i-stamp_radius)**2 + (j-stamp_radius)**2 )
        stamp[i,j] = dist.compute_pdf( distance )
    return stamp

def KDE_align( x_vals, y_vals, n=3000, min_x=None, max_x=None, min_y=None, max_y=None ):
    # Transform input values to 0-1 scale
    if min_x == None: min_x = np.min( x_vals )
    if max_x == None: max_x = np.max( x_vals )
    range_x = max_x - min_x
    if min_y == None: min_y = np.min( y_vals )
    if max_y == None: max_y = np.max( y_vals )
    range_y = max_y - min_y

    # Discretize the normalized values to the number of array points
    xs = np.array( [ int( np.round( (n-1) * (x-min_x)/range_x ) )
                     for x in np.array(x_vals) ] )
    ys = np.array( [ int( np.round( (n-1) * (y-min_y)/range_y ) )
                     for y in np.array(y_vals) ] )

    # Compute Silverman kernel and define the "stamp"
    bandwidth = np.power( float(len(xs)), -1/6 ) * ( np.std(xs) + np.std(ys) ) / 2
    kernel_sd = bandwidth / 2 / np.sqrt( 2 * np.log(2) )
    distribution = cosine_gaussian( sd=kernel_sd )
    stamp = define_stamp( distribution )
    stamp_width = stamp.shape[0]
    stamp_radius = (stamp_width-1)/2

    # Generate the KDE map
    array = np.zeros( (n,n) )
    for i,j in zip( xs, ys ):
        l_bound = int( np.max( [ 0, i-stamp_radius ] ) )
        r_bound = int( np.min( [ n-1, i+stamp_radius ] ) )
        b_bound = int( np.max( [ 0, j-stamp_radius ] ) )
        t_bound = int( np.min( [ n-1, j+stamp_radius ] ) )

        if l_bound == 0:
            l_trim = stamp_width - r_bound - 1
            r_trim = stamp_width
        elif r_bound == n-1:
            l_trim = 0
            r_trim = n - l_bound
        else:
            l_trim = 0
            r_trim = stamp_width
        if b_bound == 0:
            b_trim = stamp_width - t_bound - 1
            t_trim = stamp_width
        elif t_bound == n-1:
            b_trim = 0
            t_trim = n - b_bound
        else:
            b_trim = 0
            t_trim = stamp_width

        local_stamp = stamp[ l_trim:r_trim, b_trim:t_trim ]
        array[ l_bound:r_bound+1, b_bound:t_bound+1 ] += local_stamp

    # Identify the KDE apex and perform forward ridge walk
    peak_i, peak_j = np.unravel_index( array.argmax(), array.shape)
    trace_i = peak_i
    trace_j = peak_j
    steps = [ (0,1), (1,0), (1,1) ]
    fit_points = [ (peak_i, peak_j) ]
    while trace_i < n-1 and trace_j < n-1:
        best_step_index = np.argmax( [ array[ trace_i+step_i, trace_j+step_j ]
                                     for step_i, step_j in steps  ] )
        trace_i += steps[best_step_index][0]
        trace_j += steps[best_step_index][1]
        if best_step_index > 0:
            fit_points.append( [ trace_i, trace_j ] )
    
    # Begin reverse ridge walk
    trace_i = peak_i
    trace_j = peak_j
    while trace_i > 0 and trace_j > 0:
        best_step_index = np.argmax( [ array[ trace_i-step_i, trace_j-step_j ]
                                     for step_i, step_j in steps  ] )
        trace_i -= steps[best_step_index][0]
        trace_j -= steps[best_step_index][1]
        if best_step_index > 0:
            fit_points.append( [ trace_i, trace_j ] )
            
    # Sort ridge walk along x-axis, reshape + rescale values, and compute spline
    fit_points.sort( key=lambda x:x[0] )
    if fit_points[0][0] != 0: fit_points = [ [0, fit_points[0][1]] ] + fit_points
    if fit_points[-1][0] != n-1: fit_points = fit_points + [ [n-1, fit_points[-1][1]] ]
    fit_x, fit_y = np.array(fit_points).T / float(n-1)
    fit_x = fit_x*range_x + min_x
    fit_y = fit_y*range_y + min_y
    spline = interp1d( fit_x, fit_y, kind='slinear' )
    interpolation_range = ( np.min( fit_x ), np.max( fit_x ) )
    return spline, interpolation_range, array, fit_x, fit_y


def align_rt_dicts( rt_obs, rt_ref, n=3000, min_len=7, max_len=30,  ):
    # rt_obs and rt_ref both need to be dicts of peptide_modseq to RT
    shared_modseqs = sorted( set( rt_obs ) & set( rt_ref ) )
    shared_modseqs = [ p for p in shared_modseqs 
                       if len( modseq_to_seq( p ) ) >= min_len
                       and len( modseq_to_seq( p ) ) <= max_len ]

    shared_rt_obs = np.array( [ rt_obs[x] for x in shared_modseqs ] )
    shared_rt_ref = np.array( [ rt_ref[x] for x in shared_modseqs ] )

    spline, bounds, array, fit_x, fit_y = KDE_align( shared_rt_obs,
                                                     shared_rt_ref,
                                                     n=n, )
    clipped_rt_obs = np.clip( list( rt_obs.values() ), *bounds, )
    aligned_rt_obs = [ float( spline( rt ) ) for rt in clipped_rt_obs ]
    aligned_rt_dict = dict( zip( list(rt_obs), aligned_rt_obs ) )

    return aligned_rt_dict, spline, bounds

