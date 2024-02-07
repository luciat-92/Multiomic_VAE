import sys, os
import torch.autograd as autograd
from itertools import chain
from collections import defaultdict
from models.model_definitions import *
from train.helper_eval import *
from train.helper_utils import *
from utils.logger import *

def dsn_ae_train_step(s_dsnae, t_dsnae, s_batch, t_batch, device, optimizer, history, scheduler=None):
    """
    The function `dsn_ae_train_step` performs a training step for two autoencoders, `s_dsnae` and
    `t_dsnae`, using source and target batches, and updates the optimizer and scheduler.
    
    :param s_dsnae: s_dsnae is the source domain-specific neural autoencoder model
    :param t_dsnae: The parameter `t_dsnae` is an instance of a target domain-specific neural
    autoencoder model
    :param s_batch: s_batch is a batch of data from the source domain. It contains the input data and
    the corresponding labels or targets
    :param t_batch: The `t_batch` parameter represents the batch of data from the target domain. It is
    used in the training step of a domain-specific neural autoencoder (DSNAE) model
    :param device: The "device" parameter is used to specify the device (e.g., CPU or GPU) on which the
    computations should be performed. It is typically a torch.device object
    :param optimizer: The optimizer is an optimization algorithm that updates the parameters of the
    models (s_dsnae and t_dsnae) based on the computed gradients. It is used to minimize the loss
    function
    :param history: The "history" parameter is a dictionary that keeps track of the loss values during
    training. It is updated with the loss values for each training step
    :param scheduler: The scheduler parameter is an optional learning rate scheduler that adjusts the
    learning rate during training. It can be used to dynamically change the learning rate based on the
    training progress. If a scheduler is provided, it will be called after each training step to update
    the learning rate
    :return: the updated history dictionary.
    """
    s_dsnae.zero_grad()
    t_dsnae.zero_grad()
    s_dsnae.train() # set the model to training mode
    t_dsnae.train() # set the model to training mode

    s_x = s_batch[0].to(device) 
    t_x = t_batch[0].to(device)

    s_loss_dict = s_dsnae.loss_function(*s_dsnae(s_x)) # forward pass + loss computation
    t_loss_dict = t_dsnae.loss_function(*t_dsnae(t_x))  # forward pass + loss computation

    optimizer.zero_grad() # clear the gradients
    loss = s_loss_dict['loss'] + t_loss_dict['loss'] # compute the total loss
    loss.backward() # backward pass

    optimizer.step() # update the model parameters
    if scheduler is not None:
        scheduler.step()
    loss_dict = {k: v.cpu().detach().item() + t_loss_dict[k].cpu().detach().item() for k, v in s_loss_dict.items()}

    for k, v in loss_dict.items():
        history[k].append(v)

    return history

def compute_gradient_penalty(critic, real_samples, fake_samples, device):
    """
    The function `compute_gradient_penalty` calculates the gradient penalty loss for the Wasserstein GAN
    with gradient penalty (WGAN-GP) algorithm.
    
    :param critic: The critic is a neural network model that is used to evaluate the quality of
    generated samples. It takes in samples (real or fake) as input and outputs a score indicating how
    realistic the samples are
    :param real_samples: The real_samples parameter is a tensor containing the real samples from the
    dataset
    :param fake_samples: The `fake_samples` parameter represents the generated samples from the
    generator network in a GAN (Generative Adversarial Network). These samples are usually generated by
    passing random noise through the generator network
    :param device: The `device` parameter is used to specify whether the computation should be performed
    on a CPU or a GPU. It is typically a `torch.device` object, such as `torch.device("cuda")` for GPU
    or `torch.device("cpu")` for CPU
    :return: the gradient penalty, which is a scalar value representing the penalty for the gradient of
    the critic network.
    """
    # Random weight term for interpolation between real and fake samples
    alpha = torch.rand((real_samples.shape[0], 1)).to(device)
    # Get random interpolation between real and fake samples
    interpolates = (alpha * real_samples + ((1 - alpha) * fake_samples)).requires_grad_(True)
    critic_interpolates = critic(interpolates)
    fakes = torch.ones((real_samples.shape[0], 1)).to(device)
    # Get gradient w.r.t. interpolates
    gradients = autograd.grad(
        outputs=critic_interpolates,
        inputs=interpolates,
        grad_outputs=fakes,
        create_graph=True,
        retain_graph=True,
        only_inputs=True,
    )[0]
    gradients = gradients.view(gradients.size(0), -1)
    gradient_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean()
    return gradient_penalty

def critic_dsn_train_step(critic, s_dsnae, t_dsnae, s_batch, t_batch, device, optimizer, 
                          history, scheduler=None, clip=None, gp=None, only_shared=False):
    """
    The function `critic_dsn_train_step` trains a critic model using the DSN-AE architecture, 
    with optional gradient penalty and gradient clipping.
    
    :param critic: The critic is a neural network model that is used to evaluate the quality of
    generated samples compared to real samples. It takes as input the encoded representations of samples
    from the source and target domains
    :param s_dsnae: s_dsnae is an autoencoder model for the source domain
    :param t_dsnae: The parameter "t_dsnae" is a neural network model that represents the target domain
    source-to-target domain autoencoder. It is used to encode the input data from the target domain
    :param s_batch: s_batch is a tuple containing the source data batch. The first element of the tuple
    is the source data
    :param t_batch: The `t_batch` parameter represents the batch of target data samples. It is used in
    the `critic_dsn_train_step` function to compute the code representation (`t_code`) of the target
    data using the `t_dsnae` model
    :param device: The "device" parameter is used to specify whether the computations should be
    performed on a CPU or a GPU. It is typically a torch.device object, and you can pass either "cuda"
    or "cuda:0" to use the first available GPU, or "cpu" to use the CPU
    :param optimizer: The optimizer is an instance of the optimizer class that is used to update the
    parameters of the critic model during training. It is responsible for computing and applying the
    gradients to the model parameters based on the loss function. Examples of optimizers include Adam,
    SGD, and RMSprop
    :param history: The `history` parameter is a dictionary that keeps track of the training progress.
    It is updated with the critic loss value after each training step
    :param scheduler: The `scheduler` parameter is an optional learning rate scheduler that adjusts the
    learning rate during training. It can be used to change the learning rate based on the number of
    training steps or epochs
    :param clip: The `clip` parameter is used to clip the gradients during training. It specifies the
    maximum value for the gradient norm. If the norm of the gradients exceeds this value, the gradients
    are scaled down to ensure they do not explode
    :param gp: The parameter "gp" stands for gradient penalty. It is used to control the strength of the
    gradient penalty term in the loss function. The gradient penalty is a regularization term that
    encourages the gradients of the critic network to have a norm of 1. It helps to stabilize the
    training of the critic network
    :param only_shared: The `only_shared` parameter is a boolean flag that specifies whether to use only
    the shared part of the encoder for the adversarial training. If set to True, only the shared part of
    the encoder is used for the adversarial training. If set to False, the entire encoder (private 
    concatenated with shared) is used
    :return: the updated history dictionary.
    """
    critic.zero_grad()
    s_dsnae.zero_grad()
    t_dsnae.zero_grad()
    s_dsnae.eval()
    t_dsnae.eval()
    critic.train()

    s_x = s_batch[0].to(device)
    t_x = t_batch[0].to(device)


    if only_shared:
        s_code = s_dsnae.s_encode(s_x)
        t_code = t_dsnae.s_encode(t_x)
    else:
        s_code = s_dsnae.encode(s_x)
        t_code = t_dsnae.encode(t_x)
    
    loss = torch.mean(critic(t_code)) - torch.mean(critic(s_code))

    if gp is not None:
        gradient_penalty = compute_gradient_penalty(critic,
                                                    real_samples=s_code,
                                                    fake_samples=t_code,
                                                    device=device)
        loss = loss + gp * gradient_penalty

    optimizer.zero_grad()
    loss.backward()
    #     if clip is not None:
    #         torch.nn.utils.clip_grad_norm_(model.parameters(), clip)
    optimizer.step()

    if clip is not None:
        for p in critic.parameters():
            p.data.clamp_(-clip, clip)
    if scheduler is not None:
        scheduler.step()

    history['critic_loss'].append(loss.cpu().detach().item())

    return history


def gan_dsn_gen_train_step(critic, s_dsnae, t_dsnae, s_batch, t_batch, device, optimizer, alpha, history, scheduler=None, only_shared=False):
    """
    The function `gan_dsn_gen_train_step` performs a single training step for a generative adversarial
    network with domain-specific neural autoencoders.
    
    :param critic: The critic is a neural network model that evaluates the quality of generated samples.
    It takes as input the encoded representation of the target samples and outputs a scalar value
    indicating the quality
    :param s_dsnae: s_dsnae is an instance of the source domain domain-specific autoencoder (s_dsnae)
    model
    :param t_dsnae: t_dsnae is an instance of a target domain domain-specific autoencoder
    :param s_batch: s_batch is a batch of samples from the source domain. It contains the input data and
    possibly other information such as labels or target values
    :param t_batch: The `t_batch` parameter represents a batch of data from the target domain. It is
    used in the training step of the GAN-DSN model
    :param device: The "device" parameter is used to specify the device (CPU or GPU) on which the
    computations should be performed. It is typically a torch.device object
    :param optimizer: The optimizer is an instance of the optimizer class that is used to update the
    parameters of the models during training. It is responsible for computing and applying the gradients
    to the model parameters based on the loss function. Examples of optimizers include Adam, SGD, and
    RMSprop
    :param alpha: The parameter "alpha" is a weighting factor that determines the importance of the
    generator loss (from GAN) compared to the reconstruction loss (from AE). It is used to balance the trade-off between
    generating realistic samples (closest to cell lines) and preserving the original information during training. 
    A higher value of alpha will prioritize the generator loss, while a lower value will prioritize
    :param history: The `history` parameter is a dictionary that keeps track of the loss values during
    training. It is updated with the loss values for each training step
    :param scheduler: The "scheduler" parameter is an optional learning rate scheduler that can be used
    to adjust the learning rate during training. It is used to update the learning rate based on a
    predefined schedule. If a scheduler is provided, it will be called after each training step to
    update the learning rate
    :param only_shared: The `only_shared` parameter is a boolean flag that specifies whether to use only
    the shared part of the encoder for the adversarial training. If set to True, only the shared part of
    the encoder is used for the adversarial training. If set to False, the entire encoder (private 
    concatenated with shared) is used
    :return: the updated history dictionary.
    """
    critic.zero_grad()
    s_dsnae.zero_grad()
    t_dsnae.zero_grad()
    critic.eval()
    s_dsnae.train()
    t_dsnae.train()

    s_x = s_batch[0].to(device)
    t_x = t_batch[0].to(device)

    if only_shared:
        t_code = t_dsnae.s_encode(t_x)
    else:
        t_code = t_dsnae.encode(t_x)
    
    optimizer.zero_grad()
    gen_loss = -torch.mean(critic(t_code))
    s_loss_dict = s_dsnae.loss_function(*s_dsnae(s_x))
    t_loss_dict = t_dsnae.loss_function(*t_dsnae(t_x))
    recons_loss = s_loss_dict['loss'] + t_loss_dict['loss'] # called recon_loss BUT formed as recons_loss + self.alpha * ortho_loss (for each data source)
    loss = recons_loss + alpha * gen_loss
    optimizer.zero_grad()

    loss.backward()
    optimizer.step()
    if scheduler is not None:
        scheduler.step()

    loss_dict = {k: v.cpu().detach().item() + t_loss_dict[k].cpu().detach().item() for k, v in s_loss_dict.items()}

    for k, v in loss_dict.items():
        history[k].append(v)
    history['gen_loss'].append(gen_loss.cpu().detach().item())

    return history


def train_code_adv(s_dataloaders, t_dataloaders, **kwargs):
    """
    The function `train_code_adv` trains a deep learning model using adversarial training with domain
    adaptation.
    
    :param s_dataloaders: A list containing the source domain training and testing dataloaders
    :param t_dataloaders: The `t_dataloaders` parameter is a list containing two dataloaders:
    `t_train_dataloader` and `t_test_dataloader`. These dataloaders are used for training and testing
    the target domain model
    :return: The function `train_code_adv` returns the trained shared encoder and a tuple containing the
    training history for the DSNAE model, validation history for the DSNAE model, training history for
    the critic model, and training history for the generator model.
    """
   
    s_train_dataloader = s_dataloaders[0]
    s_test_dataloader = s_dataloaders[1]

    t_train_dataloader = t_dataloaders[0]
    t_test_dataloader = t_dataloaders[1]

    shared_encoder = MLP(input_dim=kwargs['input_dim'],
                         output_dim=kwargs['latent_dim'],
                         hidden_dims=kwargs['encoder_hidden_dims'],
                         dop=kwargs['dop']).to(kwargs['device'])

    shared_decoder = MLP(input_dim=2 * kwargs['latent_dim'],
                         output_dim=kwargs['input_dim'],
                         hidden_dims=kwargs['encoder_hidden_dims'][::-1],
                         dop=kwargs['dop']).to(kwargs['device'])

    s_dsnae = DSNAE(shared_encoder=shared_encoder,
                    decoder=shared_decoder,
                    alpha=kwargs['alpha'],
                    input_dim=kwargs['input_dim'],
                    latent_dim=kwargs['latent_dim'],
                    hidden_dims=kwargs['encoder_hidden_dims'],
                    dop=kwargs['dop'],
                    norm_flag=kwargs['norm_flag']).to(kwargs['device'])

    t_dsnae = DSNAE(shared_encoder=shared_encoder,
                    decoder=shared_decoder,
                    alpha=kwargs['alpha'],
                    input_dim=kwargs['input_dim'],
                    latent_dim=kwargs['latent_dim'],
                    hidden_dims=kwargs['encoder_hidden_dims'],
                    dop=kwargs['dop'],
                    norm_flag=kwargs['norm_flag']).to(kwargs['device'])

    # used as critic in the GAN training
    critic_input_dim = kwargs['latent_dim'] if kwargs['only_shared'] else kwargs['latent_dim'] * 2
    confounding_classifier = MLP(input_dim=critic_input_dim,
                                 output_dim=1,
                                 hidden_dims=kwargs['classifier_hidden_dims'],
                                 dop=kwargs['dop']).to(kwargs['device'])

    ae_params = [t_dsnae.private_encoder.parameters(),
                 s_dsnae.private_encoder.parameters(),
                 shared_decoder.parameters(),
                 shared_encoder.parameters()
                 ]
    t_ae_params = [t_dsnae.private_encoder.parameters(),
                   s_dsnae.private_encoder.parameters(),
                   shared_decoder.parameters(),
                   shared_encoder.parameters()
                   ]

    ae_optimizer = torch.optim.AdamW(chain(*ae_params), lr=kwargs['lr']) # used in the pretraining
    classifier_optimizer = torch.optim.RMSprop(confounding_classifier.parameters(), lr=kwargs['lr']) # used in the GAN training (why not adam?)
    t_ae_optimizer = torch.optim.RMSprop(chain(*t_ae_params), lr=kwargs['lr']) # used in the GAN training (why not adam?)

    dsnae_train_history = defaultdict(list)
    dsnae_val_history = defaultdict(list)
    critic_train_history = defaultdict(list)
    gen_train_history = defaultdict(list)

    if kwargs['retrain_flag']:

        # start dsnae pre-training (without adversarial training)
        for epoch in range(int(kwargs['pretrain_num_epochs'])):
            if epoch % 50 == 0:
                print(f'AE training epoch {epoch}')
            for step, s_batch in enumerate(s_train_dataloader):
                t_batch = next(iter(t_train_dataloader))
                dsnae_train_history = dsn_ae_train_step(s_dsnae=s_dsnae,
                                                        t_dsnae=t_dsnae,
                                                        s_batch=s_batch,
                                                        t_batch=t_batch,
                                                        device=kwargs['device'],
                                                        optimizer=ae_optimizer,
                                                        history=dsnae_train_history)
            # NOTE: currently, the validation is part of the training data, NOT PROPER!
            dsnae_val_history = eval_dsnae_epoch(model=s_dsnae,
                                                 data_loader=s_test_dataloader,
                                                 device=kwargs['device'],
                                                 history=dsnae_val_history
                                                 )
            dsnae_val_history = eval_dsnae_epoch(model=t_dsnae,
                                                 data_loader=t_test_dataloader,
                                                 device=kwargs['device'],
                                                 history=dsnae_val_history
                                                 )
            # TODO: remove for loop?
            for k in dsnae_val_history:
                if k != 'best_index':
                    dsnae_val_history[k][-2] += dsnae_val_history[k][-1]
                    dsnae_val_history[k].pop()
            if kwargs['es_flag']:
                save_flag, stop_flag = model_save_check(dsnae_val_history, metric_name='loss', tolerance_count=20)
                if save_flag: # save models parameters
                    torch.save(s_dsnae.state_dict(), os.path.join(kwargs['model_save_folder'], 'a_s_dsnae.pt')) 
                    torch.save(t_dsnae.state_dict(), os.path.join(kwargs['model_save_folder'], 'a_t_dsnae.pt'))
                if stop_flag:
                    break
        if kwargs['es_flag']: # load model if pretraining was already computed
            s_dsnae.load_state_dict(torch.load(os.path.join(kwargs['model_save_folder'], 'a_s_dsnae.pt')))
            t_dsnae.load_state_dict(torch.load(os.path.join(kwargs['model_save_folder'], 'a_t_dsnae.pt')))

        # start GAN training
        for epoch in range(int(kwargs['train_num_epochs'])):
            if epoch % 50 == 0:
                print(f'confounder wgan training epoch {epoch}')
            for step, s_batch in enumerate(s_train_dataloader):
                t_batch = next(iter(t_train_dataloader))
                critic_train_history = critic_dsn_train_step(critic=confounding_classifier,
                                                             s_dsnae=s_dsnae,
                                                             t_dsnae=t_dsnae,
                                                             s_batch=s_batch,
                                                             t_batch=t_batch,
                                                             device=kwargs['device'],
                                                             optimizer=classifier_optimizer,
                                                             history=critic_train_history,
                                                             # clip=0.1,
                                                             gp=10.0, # gradient penalty set from wGAN with gradient penality paper
                                                             only_shared=kwargs['only_shared'])
                if (step + 1) % 5 == 0: # every 5 steps, update of the encoder/decoder param
                    gen_train_history = gan_dsn_gen_train_step(critic=confounding_classifier,
                                                               s_dsnae=s_dsnae,
                                                               t_dsnae=t_dsnae,
                                                               s_batch=s_batch,
                                                               t_batch=t_batch,
                                                               device=kwargs['device'],
                                                               optimizer=t_ae_optimizer,
                                                               alpha=1.0,
                                                               history=gen_train_history, 
                                                               only_shared=kwargs['only_shared'])

        torch.save(s_dsnae.state_dict(), os.path.join(kwargs['model_save_folder'], 'a_s_dsnae.pt'))
        torch.save(t_dsnae.state_dict(), os.path.join(kwargs['model_save_folder'], 'a_t_dsnae.pt'))

    else:
        try:
            t_dsnae.load_state_dict(torch.load(os.path.join(kwargs['model_save_folder'], 'a_t_dsnae.pt'))) # load already computed model

        except FileNotFoundError:
            raise Exception("No pre-trained encoder")

    return t_dsnae.shared_encoder, (dsnae_train_history, dsnae_val_history, critic_train_history, gen_train_history)
