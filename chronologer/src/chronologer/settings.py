import torch


hyperparameters = { 'embed_dimension' : 64,
                    'n_resnet_blocks' : 3,
                    'kernel_size' : 7,
                    'activation_function' : 'relu',
                  }

training_parameters = { 'n_epochs' : 100,
                        'learning_rate' : 1e-3,
                        'dropout_rate' : 0.1,
                        'initial_batch_size' : 64,
                        'epochs_to_2x_batch' : 30,
                        'max_batch_size' : 1024,
                        'optimizer' : torch.optim.Adam,
                        'train_device' : 'cuda',
                        'eval_device' : 'cpu', }

