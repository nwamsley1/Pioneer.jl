import numpy as np
import torch
import torch.nn as nn

import chronologer.src.constants as constants
from chronologer.src.distributions import return_distribution_dict


def generate_outlier_mask( array, dist_family, fdr, ):
    dist_dict = return_distribution_dict()
    assert dist_family in dist_dict, 'UNDEFINED DISTRIBUTION ' + dist_family
    dist = dist_dict[ dist_family ]( data = array ) # Initialize with data
    threshold = dist.ppf( 1.0 - fdr, )
    return array < threshold

def fdr_masked_mean( array, dist_family, fdr, eps=constants.epsilon, ):
    outlier_mask = generate_outlier_mask( array, dist_family, fdr, )
    masked_mean = ( torch.sum( array * outlier_mask ) + eps ) /\
                  ( torch.sum( outlier_mask ) + eps )
    return masked_mean


def l2_norm( vector, eps = constants.epsilon, ):
    sum_sqs = torch.sum( vector**2, dim=(1,-1), keepdims=True, )
    denom = torch.sqrt( sum_sqs.clamp( eps, ) )
    #print( vector.shape, sum_sqs.shape, denom.shape, )
    return vector / denom

def neg_logit( array, eps = constants.epsilon, ):
    return torch.log( (1 - array).clamp(eps,) ) - torch.log( array.clamp(eps,) ) 


class MAE_Loss( nn.Module ):
    def forward( self, pred, true, weight, eps = constants.epsilon, ):
        return torch.mean( torch.abs( pred - true ) )
    
class MSE_Loss( nn.Module ):
    def forward( self, pred, true, weight, eps = constants.epsilon, ):
        return torch.mean( torch.square( pred - true ) )
        

class LogL_Loss( nn.Module ):
    def __init__( self, n_sources=1, family='laplace', fdr=constants.train_fdr ):
        super().__init__()
        self.source_scale = nn.Linear( n_sources, 1, bias=False, )
        nn.init.constant_( self.source_scale.weight, 10.0, )
        self.family = family
        self.dist = return_distribution_dict()[ family ]
        self.fdr = fdr
        
    def forward( self, pred, true, source, eps = constants.epsilon, ):
        scale = self.source_scale( source ).clamp( eps, )
        dist = self.dist( center=true, scale=scale, )
        logL_loss = -1 * dist.logL( pred )
        return fdr_masked_mean( logL_loss, self.family, self.fdr, )
        


class RT_masked_negLogL( nn.Module ):
    def __init__( self, n_sources, fdr = constants.train_fdr, ):
        super().__init__()
        self.source_b = nn.Linear( n_sources, 1, bias=False, )
        nn.init.constant_( self.source_b.weight, 10.0, )
        self.fdr = fdr
    
    def forward( self, pred, true, source, eps = constants.epsilon, ):
        abs_error = torch.abs( pred - true )
        exp_abs_error = self.source_b( source ).clamp( eps, )
        norm_abs_error = abs_error / exp_abs_error
        neg_logL = norm_abs_error + torch.log( 2 * exp_abs_error )
        
        outlier_mask = generate_outlier_mask( norm_abs_error, 'laplace', self.fdr, )
        mean_neg_logL = ( torch.sum( neg_logL * outlier_mask ) + eps ) /\
                        ( torch.sum( outlier_mask ) + eps )

        return mean_neg_logL
    
    
class Spectrum_masked_negLogit( nn.Module ):
    def __init__( self, fdr = constants.train_fdr, ):
        super().__init__()
        self.fdr = fdr

    def forward( self, pred, true, weights, eps = constants.epsilon, ):
        ion_mask =  ( true + 1.0 ) / ( true + 1.0 + eps ) # Incompatible ions = -1

        pred_masked = pred * ion_mask
        true_masked = true * ion_mask
        pred_norm = l2_norm( pred_masked, eps, )
        true_norm = l2_norm( true_masked, eps, )
        product = torch.sum( pred_norm * true_norm, dim=(1,-1), )
        score = neg_logit( product, )
        
        outlier_mask = generate_outlier_mask( score, 'gumbel', 0.001, )
        mean_neg_logit = ( torch.sum( score * outlier_mask ) + eps ) /\
                         ( torch.sum( outlier_mask ) + eps )
                         
        return mean_neg_logit
        


class Masked_spectral_angle_loss( nn.Module ):
    def __init__( self, fdr = constants.train_fdr, ):
        super().__init__()
        self.fdr = fdr

    def forward( self, pred, true, weights, eps = constants.epsilon, ):
        ion_mask =  ( true + 1.0 ) / ( true + 1.0 + eps ) # Incompatible ions = -1

        pred_masked = pred * ion_mask
        true_masked = true * ion_mask
        pred_norm = l2_norm( pred_masked, eps, )
        true_norm = l2_norm( true_masked, eps, )
        product = torch.sum( pred_norm * true_norm, dim=(1,-1), )
        arccos = torch.arccos( product )
        score = 2 * arccos / np.pi
        
        
        outlier_mask = generate_outlier_mask( score, 'gumbel', self.fdr, )
        mean_neg_logit = ( torch.sum( score * outlier_mask ) + eps ) /\
                         ( torch.sum( outlier_mask ) + eps )
                         
        return mean_neg_logit