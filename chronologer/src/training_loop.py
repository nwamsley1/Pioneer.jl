
import numpy as np
from copy import deepcopy
import time, datetime

import torch
from torch.utils.data import DataLoader, RandomSampler
from torch.optim.lr_scheduler import ReduceLROnPlateau

from encyclopydia.src.utils import timer

def train_model( model, 
                 datasets, 
                 initial_batch_size, 
                 max_batch_size, 
                 epochs_to_double_batch,
                 loss_fx, 
                 optimizer,
                 num_epochs, 
                 train_device, 
                 other_device,
                 file_name,
                 start_time,
                 initial_loss = 1e40,
                 score_fxs = {}, ):
    
    s_time = time.time()
    
    phases = list( datasets )
    # Store scores as nested dicts: phase -> fx -> list by epoch
    all_scores = dict( [ ( phase, dict( [ ( name, [] ) 
                                          for name in score_fxs ] ) ) 
                         for phase in phases ] )
    
    batch_sizes = dict( [ ( p, 0 ) if p == 'train' 
                          else ( p, max_batch_size ) for p in phases ] )
    devices = dict( [ ( p, train_device ) if p == 'train' 
                      else ( p, other_device ) for p in phases ] )
    
    best_epoch = 0
    best_loss = initial_loss
    tolerance = 1e-4

    for epoch in range( 1, num_epochs+1 ):
        print( 'Epoch ' + str(epoch) + ' of ' + str(num_epochs) )
        print( '-' * 50 )
        
        float_batch_scaler = initial_batch_size * np.exp( np.log(2) * (epoch-1) / epochs_to_double_batch  ) / 8
        train_batch_size = int( round( float_batch_scaler ) ) * 8
        if train_batch_size != batch_sizes['train']:
            batch_sizes['train'] = train_batch_size
            dataloaders = dict( [ ( p, DataLoader( datasets[p], batch_sizes[p], shuffle=True, ) ) 
                                for p in phases ] )
            
        for phase in phases:
            model.to( devices[ phase ] )
            loss_fx.to( devices[ phase ] )
            # score_fxs is a dict with (name,fx)
            for name in score_fxs:
                score_fxs[ name ].to( devices[ phase ] )
            
            if phase == 'train':
                model.train()  # Set model to training mode
                print( 'Batch size = ' + str(batch_sizes['train']) )
            else:
                model.eval()   # Set model to evaluate mode

            running_loss = 0.0
            running_scores = np.zeros( len(score_fxs) )
            total_samples = 0

            data = dataloaders[phase]

            # Iterate over data.
            for i, batch in enumerate( data ):
                batch_size = batch[0].size(0)
                batch = [ b.to( devices[phase] ) for b in batch ]
                inputs = batch[:-2]
                outputs = batch[-2:] # y and weight/source
                
                pred = model( *inputs )                    
                loss = loss_fx( pred, *outputs, )
                scores = [ score_fxs[name]( pred, *outputs, ).item() for name in score_fxs ]
                
                if phase == 'train':
                    # zero the parameter gradients
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                    
                # statistics
                running_loss += loss.item() * batch_size
                running_scores += np.array( scores ) * batch_size
                total_samples += batch_size

            epoch_loss = running_loss / total_samples
            epoch_scores = running_scores / total_samples
            for i, name in enumerate( score_fxs ):
                all_scores[phase][name].append( epoch_scores[i] )
            runtime = time.time() - s_time
            print( phase.capitalize() + format( epoch_loss, '.4f' ).rjust(8) )

            if phase == 'test':
                #MAEs = loss_fx.source_b.weight.cpu().detach().numpy().tolist()[0]
                #for t, learned_mae in enumerate( MAEs ):
                #    print( '\t' + unique_sources[t].ljust(25) + format(learned_mae,'.3f') )
                if epoch_loss < best_loss-tolerance:
                    print( 'New best weights! Copying and saving model to\n\t' + file_name )
                    best_epoch = epoch
                    best_loss = epoch_loss
                    torch.save( model.state_dict(), file_name )
                    epochs_wo_improv = 0
                else:
                    print( 'Did not improve, best performance was epoch ' + 
                           str(best_epoch) + ' (' + format(best_loss,'.4f') + ')' )
        runtime = time.time() - s_time
        print( 'Runtime: ' + timer( start_time ) + '\n' )
    
    return best_loss, all_scores


def nce_optimization_routine( cartographer_beam_model,
                              charge_to_nce_model, 
                              dataset,
                              batch_size,
                              loss_fx, 
                              optimizer,
                              device,
                              start_time,
                              score_threshold = 1e-3,
                              patience = 20,
                              epoch_per_update = 5, ):
    ## Sample data until max change in NCE is < threshold
    cartographer_beam_model.to(device)
    cartographer_beam_model.train()
    charge_to_nce_model.to(device)
    charge_to_nce_model.train()
    loss_fx.to(device)
    
    data_loader = DataLoader( dataset, batch_size=batch_size, )
    
    lr_scheduler = ReduceLROnPlateau( optimizer,
                                      factor = 0.5,
                                      patience = 3, 
                                      threshold = score_threshold,
                                      threshold_mode = 'abs', )
    
        
    best_loss = 1e7
    epochs_wo_improve = 0
    epoch = 0
    while epochs_wo_improve < patience:
        epoch += 1
        total_loss = 0.0
        total_n = 0
        for _, data in enumerate( data_loader ):
            peptide_tensor, charge_tensor, intensity_tensor = [ x.to(device) for x in data ]
            nce_tensor = charge_to_nce_model( charge_tensor )
            pred_ladder_intensities = cartographer_beam_model( peptide_tensor, 
                                                charge_tensor, 
                                                nce_tensor )
            loss = loss_fx( pred_ladder_intensities, intensity_tensor, torch.ones((1,1)) )
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            n = peptide_tensor.size(0)
            total_loss += loss.item() * n
            total_n += n
        
        lr_scheduler.step( total_loss, )
        mean_loss = total_loss / total_n
        if mean_loss < best_loss - score_threshold:
            best_loss = mean_loss
            best_state_dict = deepcopy( charge_to_nce_model.state_dict() )
            best_nces = charge_to_nce_model.weight.cpu().detach().numpy()[0]
            epochs_wo_improve = 0
        else:
            charge_to_nce_model.load_state_dict( best_state_dict )
            epochs_wo_improve += 1
        if epoch%epoch_per_update == 0:
            print( str('\tEpoch ' + str( epoch ) + ':').ljust(12) + 
                   str('Score = ' + format( best_loss, '.3f' )).ljust(16) +
                   '('+timer(start_time)+')' )
    
    print( str('Epoch ' + str( epoch ) + ':').ljust(12) + 
           'Score = ' + format( best_loss, '.3f' ) )     
    return best_nces
            
        
        

