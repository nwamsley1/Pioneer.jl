
import torch
import torch.nn as nn


def activation_func(activation):
    return  nn.ModuleDict( [ ['relu', nn.ReLU(inplace=True)],
                             ['leaky_relu', nn.LeakyReLU(negative_slope=0.01, inplace=True)],
                             ['selu', nn.SELU(inplace=True)],
                             ['none', nn.Identity()] ] )[ activation ]


def convolution_unit( in_chs, out_chs, kernel, d_rate, ):
    conv_layer = nn.Conv1d( in_chs, 
                            out_chs, 
                            kernel, 
                            dilation=d_rate, 
                            padding='same' )
    torch.nn.init.kaiming_normal_( conv_layer.weight, 
                                   mode='fan_out', 
                                   nonlinearity='relu' )
    norm_layer = nn.BatchNorm1d( out_chs, )
    torch.nn.init.constant_( norm_layer.weight, 1.0, )
    torch.nn.init.constant_( norm_layer.bias, 0.0, )
    return nn.Sequential( conv_layer,
                          norm_layer, )


def resnet_unit( in_chs, out_chs,  kernel, d_rate, act_fx, ):
    return nn.Sequential( convolution_unit( in_chs, 
                                            out_chs, 
                                            kernel, 
                                            d_rate, ),
                          activation_func( act_fx ) )

class resnet_block( nn.Module ):
    def __init__( self, in_channels, out_channels, kernel, d_rate, act_fx, ):
        super().__init__()
        self.in_channels, self.out_channels = in_channels, out_channels
        self.process_blocks = nn.Sequential( resnet_unit( in_channels, 
                                                          out_channels, 
                                                          1, 
                                                          d_rate,
                                                          act_fx, ),
                                             resnet_unit( out_channels, 
                                                          out_channels, 
                                                          kernel, 
                                                          d_rate,
                                                          act_fx,) )
        self.activate = activation_func( act_fx )
        self.shortcut = convolution_unit( in_channels, 
                                          out_channels, 
                                          1, 
                                          d_rate, )

        self.term_block = nn.Identity()

    def forward( self, x ):
        residual = x
        if self.in_channels != self.out_channels: 
            residual = self.shortcut( x )
        x = self.process_blocks( x )
        x = self.term_block( x )
        x = x + residual
        x = self.activate( x )
        return x


