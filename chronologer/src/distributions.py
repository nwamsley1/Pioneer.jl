import math
from scipy.special import erfinv
import torch

class gaussian:
    def __init__( self, data=None, center=torch.zeros(1), scale=torch.ones(1), ):
        if data is None:
            assert type(center) in [ float, torch.Tensor ], 'Invalid data type for center'
            self.mu = center if type(center) == torch.Tensor else torch.full( (1,), center )
            assert type(scale) in [ float, torch.Tensor ], 'Invalid data type for scale'
            self.sigma = scale if type(scale) == torch.Tensor else torch.full( (1,), scale )
        else:
            self.mu = torch.mean( data, )
            self.sigma = torch.sqrt( ( data - self.mu )**2 / 
                                        ( data.size(0) - 1 ) )
    def cdf( self, x, ):
        return 0.5 * ( 1 + torch.erf( ( x - self.mu ) / self.sigma / math.sqrt(2) ) )
    def ppf( self, q, ):
        #q = q if type(q) == torch.Tensor else torch.full( (1,), q )
        return self.mu + self.sigma * math.sqrt(2) * erfinv( 2*q - 1 )
    def logL( self, x, ):
        return -1 * ( torch.log( self.sigma*math.sqrt(2*math.pi) ) + 0.5*( (x-self.mu)**2/self.sigma  ) )    
    
class laplace:
    def __init__( self, data=None, center=torch.zeros(1), scale=torch.ones(1), ):
        # MLE of mu is median, but use mean instead so differentiable 
        # Class can be initialized either by fitting data or setting mannually
        if data is None:
            assert type(center) in [ float, torch.Tensor ], 'Invalid data type for center'
            self.mu = center if type(center) == torch.Tensor else torch.full( (1,), center )
            assert type(scale) in [ float, torch.Tensor ], 'Invalid data type for scale'
            self.b = scale if type(scale) == torch.Tensor else torch.full( (1,), scale )
        else:
            self.mu = torch.mean( data, )
            self.b = torch.mean( torch.abs( data - self.mu ), )
    def cdf( self, x, ):
        if x <= self.mu:
            return 0.5 * torch.exp( ( x - self.mu ) / self.b )
        else:
            return 1 - 0.5 * torch.exp( -( x - self.mu ) / self.b )
    def ppf( self, q, ):
        if q <= 0.5:
            return self.mu + self.b * math.log( 2*q )
        else:
            return self.mu - self.b * math.log( 2 - 2*q )
    def logL( self, x, ):
        return -1 * ( torch.log( 2*self.b ) + torch.abs( x - self.mu ) / self.b )
        
class gumbel:
    def __init__( self, data, ):
        mean = torch.mean( data, )
        std = torch.sum( ( data - mean )**2, ) / ( data.size(0) - 1 )
        self.beta = std * math.sqrt(6) / math.pi
        self.mu = mean - 0.57721 * self.beta
        
    def cdf( self, x ):
        return torch.exp( -torch.exp( -( x - self.mu ) / self.beta ) )
    
    def ppf( self, q, ):
        return self.mu - self.beta * math.log( -math.log( q ) )
        
def return_distribution_dict():
    return { 'gaussian' : gaussian,
             'laplace'  : laplace,
             'gumbel'   : gumbel, }