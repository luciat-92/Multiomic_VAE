from typing import TypeVar, List, Any
import torch
from torch import nn
from torch.nn import functional as F
from torch.autograd import Function
from abc import abstractmethod
from train.helper_eval import compute_mmd, compute_kernel

Tensor = TypeVar('torch.tensor') # variable of type Tensor is expected to be a PyTorch tensor

# The `RevGrad` class is a PyTorch module that implements a gradient reversal layer, which reverses
# the gradient in the backward pass.
class RevGrad(Function):
    @staticmethod
    def forward(ctx, input_, alpha_):
        ctx.save_for_backward(input_, alpha_)
        output = input_
        return output

    @staticmethod
    def backward(ctx, grad_output):  # pragma: no cover
        grad_input = None
        _, alpha_ = ctx.saved_tensors
        if ctx.needs_input_grad[0]:
            grad_input = -grad_output * alpha_
        return grad_input, None

revgrad = RevGrad.apply

class RevGrad(torch.nn.Module):
    def __init__(self, alpha=1., *args, **kwargs):
        """
        A gradient reversal layer.
        This layer has no parameters, and simply reverses the gradient
        in the backward pass.
        """
        super().__init__(*args, **kwargs)

        self._alpha = torch.tensor(alpha, requires_grad=False)

    def forward(self, input_):
        return revgrad(input_, self._alpha)


# The MLP class is a multi-layer perceptron model with customizable input and output dimensions,
# hidden layer dimensions, activation functions, dropout rate, and output function.
class MLP(nn.Module):

    def __init__(self, input_dim: int, output_dim: int, hidden_dims: List = None, dop: float = 0.1, act_fn=nn.SELU, out_fn=None, gr_flag=False, **kwargs) -> None:
        super(MLP, self).__init__()
        self.output_dim = output_dim
        self.dop = dop

        if hidden_dims is None:
            hidden_dims = [32, 64, 128, 256, 512]

        modules = []
        if gr_flag:
            modules.append(RevGrad())

        modules.append(
            nn.Sequential(
                nn.Linear(input_dim, hidden_dims[0], bias=True),
                #nn.BatchNorm1d(hidden_dims[0]),
                act_fn(),
                nn.Dropout(self.dop)
            )
        )

        for i in range(len(hidden_dims) - 1):
            modules.append(
                nn.Sequential(
                    nn.Linear(hidden_dims[i], hidden_dims[i + 1], bias=True),
                    #nn.BatchNorm1d(hidden_dims[i + 1]),
                    act_fn(),
                    nn.Dropout(self.dop)
                )
            )

        self.module = nn.Sequential(*modules)

        if out_fn is None:
            self.output_layer = nn.Sequential(
                nn.Linear(hidden_dims[-1], output_dim, bias=True)
            )
        else:
            self.output_layer = nn.Sequential(
                nn.Linear(hidden_dims[-1], output_dim, bias=True),
                out_fn()
            )



    def forward(self, input: Tensor) -> Tensor:
        embed = self.module(input)
        output = self.output_layer(embed)

        return output



# The above class is a base class for an autoencoder model with methods for encoding, decoding,
# sampling, generating, and defining the forward pass and loss function.
class BaseAE(nn.Module):
    def __init__(self) -> None:
        super(BaseAE, self).__init__()

    def encode(self, input: Tensor) -> List[Tensor]:
        raise NotImplementedError

    def decode(self, input: Tensor) -> List[Tensor]:
        raise NotImplementedError

    def sample(self, batch_size: int, current_device: int, **kwargs) -> Tensor:
        raise RuntimeWarning()

    def generate(self, x: Tensor, **kwargs) -> Tensor:
        raise NotImplementedError

    @abstractmethod
    def forward(self, *inputs: Tensor) -> Tensor:
        pass

    @abstractmethod
    def loss_function(self, *inputs: Any, **kwargs) -> Tensor:
        pass


    
# The DSNAE class is a neural network model that performs unsupervised learning by encoding input data
# into a latent space and reconstructing the input data from the latent space.
# The model has a shared encoder and decoder, and a private encoder.
# The model also has a loss function that includes a reconstruction loss and an orthogonality loss.
class DSNAE(BaseAE):

    def __init__(self, shared_encoder, decoder, input_dim: int, latent_dim: int, alpha: float = 1.0,
                 hidden_dims: List = None, dop: float = 0.1, noise_flag: bool = False, norm_flag: bool = False,
                 act_fn=nn.SELU, **kwargs) -> None:
        super(DSNAE, self).__init__()
        self.latent_dim = latent_dim
        self.alpha = alpha
        self.noise_flag = noise_flag
        self.dop = dop
        self.norm_flag = norm_flag
        self.act_fn = act_fn

        if hidden_dims is None:
            hidden_dims = [32, 64, 128, 256, 512]

        self.shared_encoder = shared_encoder # shared encoder
        self.decoder = decoder # shared decoder

        # build private encoder
        modules = []

        modules.append(
            nn.Sequential(
                nn.Linear(input_dim, hidden_dims[0], bias=True),
                # nn.BatchNorm1d(hidden_dims[0]),
                self.act_fn(),
                nn.Dropout(self.dop)
            )
        )

        for i in range(len(hidden_dims) - 1):
            modules.append(
                nn.Sequential(
                    nn.Linear(hidden_dims[i], hidden_dims[i + 1], bias=True),
                    # nn.BatchNorm1d(hidden_dims[i + 1]),
                    self.act_fn(),
                    nn.Dropout(self.dop)
                )
            )
        modules.append(nn.Dropout(self.dop))
        modules.append(nn.Linear(hidden_dims[-1], latent_dim, bias=True))
        # modules.append(nn.LayerNorm(latent_dim, eps=1e-12, elementwise_affine=False))

        self.private_encoder = nn.Sequential(*modules)

    def p_encode(self, input: Tensor) -> Tensor:
        if self.noise_flag and self.training:
            latent_code = self.private_encoder(input + torch.randn_like(input, requires_grad=False) * 0.1)
        else:
            latent_code = self.private_encoder(input)

        if self.norm_flag:
            return F.normalize(latent_code, p=2, dim=1)
        else:
            return latent_code

    def s_encode(self, input: Tensor) -> Tensor:
        if self.noise_flag and self.training:
            latent_code = self.shared_encoder(input + torch.randn_like(input, requires_grad=False) * 0.1)
        else:
            latent_code = self.shared_encoder(input)
        if self.norm_flag:
            return F.normalize(latent_code, p=2, dim=1)
        else:
            return latent_code

    def encode(self, input: Tensor) -> Tensor:
        p_latent_code = self.p_encode(input)
        s_latent_code = self.s_encode(input)

        return torch.cat((p_latent_code, s_latent_code), dim=1)

    def decode(self, z: Tensor) -> Tensor:
        outputs = self.decoder(z)

        return outputs

    def forward(self, input: Tensor, **kwargs) -> List[Tensor]:
        z = self.encode(input)
        return [input, self.decode(z), z]

    def loss_function(self, *args, **kwargs) -> dict:
        input = args[0]
        recons = args[1]
        z = args[2]

        p_z = z[:, :z.shape[1] // 2]
        s_z = z[:, z.shape[1] // 2:]

        recons_loss = F.mse_loss(input, recons)

        s_l2_norm = torch.norm(s_z, p=2, dim=1, keepdim=True).detach()
        s_l2 = s_z.div(s_l2_norm.expand_as(s_z) + 1e-6)

        p_l2_norm = torch.norm(p_z, p=2, dim=1, keepdim=True).detach()
        p_l2 = p_z.div(p_l2_norm.expand_as(p_z) + 1e-6)

        ortho_loss = torch.mean((s_l2.t().mm(p_l2)).pow(2)) # frobenioius norm
        loss = recons_loss + self.alpha * ortho_loss
        return {'loss': loss, 'recons_loss': recons_loss, 'ortho_loss': ortho_loss}

    def sample(self, num_samples: int, current_device: int, **kwargs) -> Tensor:
        z = torch.randn(num_samples, self.latent_dim)
        z = z.to(current_device)
        samples = self.decode(z)

        return samples

    def generate(self, x: Tensor, **kwargs) -> Tensor:
        return self.forward(x)[1]
    
########### VAE ###########
class VAE_Encoder(MLP):
    def __init__(self, input_dim: int, latent_dim: int, hidden_dims: List = None, dop: float = 0.1, act_fn=nn.SELU, out_fn=None, gr_flag=False) -> None:
        super(VAE_Encoder, self).__init__(input_dim, latent_dim * 2, hidden_dims, dop, act_fn, out_fn, gr_flag)  # Double latent_dim for mean and log_var

    def forward(self, input: Tensor) -> List[Tensor]:
        output = super().forward(input)
        mu, log_var = output.chunk(2, dim=1)  # Split output into mean and log_var
        # Reparameterization trick
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        z = mu + eps * std
        return [z, mu, log_var]

class VAE_Decoder(MLP):
    def __init__(self, latent_dim: int, output_dim: int, hidden_dims: List = None, dop: float = 0.1, act_fn=nn.SELU, out_fn=None, gr_flag=False) -> None:
        super(VAE_Decoder, self).__init__(latent_dim, output_dim, hidden_dims, dop, act_fn, out_fn, gr_flag)

    def forward(self, input: Tensor) -> Tensor:
        return super().forward(input)

#class VAE_Encoder(nn.Module):
#    def __init__(self, input_dim, latent_dim):
#        super(VAE_Encoder, self).__init__(input_dim, latent_dim * 2, hidden_dims=[512, 256], dop=0.1)
#        self.latent_dim = latent_dim
#
#    def __init__(self, input_dim, latent_dim):
#        super(VAE_Encoder, self).__init__()
#        self.input_dim = input_dim
#        self.latent_dim = latent_dim
#        self.encoder = nn.Sequential(
#            nn.Linear(input_dim, 512),
#            nn.ReLU(),
#            nn.Linear(512, 256),
#            nn.ReLU()
#        )
#        self.z_mean = nn.Linear(256, latent_dim)
#        self.z_log_var = nn.Linear(256, latent_dim)
#
#    def reparameterize(self, mu, log_var):
#        std = torch.exp(0.5 * log_var)
#        eps = torch.randn_like(std)
#        return mu + eps * std
#
#    def forward(self, x):
#        x = self.encoder(x)
#        mu = self.z_mean(x)
#        log_var = self.z_log_var(x)
#        z = self.reparameterize(mu, log_var)
#        return z, mu, log_var

#class VAE_Decoder(nn.Module):
#    def __init__(self, latent_dim, output_dim):
#        super(VAE_Decoder, self).__init__()
#        self.latent_dim = latent_dim
#        self.output_dim = output_dim
#        self.decoder = nn.Sequential(
#            nn.Linear(latent_dim, 256),
#            nn.ReLU(),
#            nn.Linear(256, 512),
#            nn.ReLU(),
#            nn.Linear(512, output_dim),
#            nn.Sigmoid()  # Output between 0 and 1
#        )
#
#    def forward(self, z):
#        x_recon = self.decoder(z)
#        return x_recon

# The BaseVAE class is an abstract base class for Variational Autoencoders (VAEs) that defines the
# basic structure and methods that need to be implemented by subclasses.
class BaseVAE(nn.Module):
    def __init__(self) -> None:
        super(BaseVAE, self).__init__()

    def encode(self, input: Tensor) -> List[Tensor]:
        raise NotImplementedError

    def reparameterize(self, mu: Tensor, log_var: Tensor) -> Tensor:
        raise NotImplementedError

    def decode(self, input: Tensor) -> List[Tensor]:
        raise NotImplementedError

    def sample(self, batch_size: int, current_device: int, **kwargs) -> Tensor:
        raise RuntimeWarning()

    def generate(self, x: Tensor, **kwargs) -> Tensor:
        raise NotImplementedError

    @abstractmethod
    def forward(self, *inputs: Tensor) -> Tensor:
        pass

    @abstractmethod
    def loss_function(self, *inputs: Any, **kwargs) -> Tensor:
        pass

class DSNVAE(BaseVAE):

    def __init__(self, shared_encoder, decoder, input_dim: int, latent_dim: int, alpha: float = 1.0, beta: float = 1.0,
                 hidden_dims: List = None, dop: float = 0.1, noise_flag: bool = False, norm_flag: bool = False, 
                 act_fn=nn.SELU, device = 'cpu', **kwargs) -> None:
        super(DSNVAE, self).__init__()
        self.latent_dim = latent_dim
        self.alpha = alpha
        self.beta = beta
        self.noise_flag = noise_flag
        self.dop = dop
        self.norm_flag = norm_flag
        self.act_fn = act_fn
        self.device = device

        if hidden_dims is None:
            hidden_dims = [32, 64, 128, 256, 512]

        self.shared_encoder = shared_encoder # shared encoder
        self.decoder = decoder # shared decoder

        # build private encoder
        self.private_encoder = VAE_Encoder(input_dim, 
                                           latent_dim, 
                                           hidden_dims, 
                                           dop, 
                                           act_fn, 
                                           out_fn=None, 
                                           gr_flag=False)
        
    # remove the normalization flag (?)
    def p_encode(self, input: Tensor) -> List[Tensor]:
        latent_code, mu, log_var = self.private_encoder(input)
        if self.norm_flag:
            return [F.normalize(latent_code, p=2, dim=1) , mu, log_var]
        else:
            return [latent_code, mu, log_var]

    # remove the normalization flag (?)
    def s_encode(self, input: Tensor) -> List[Tensor]:
        latent_code, mu, log_var = self.shared_encoder(input)
        if self.norm_flag:
            return [F.normalize(latent_code, p=2, dim=1) , mu, log_var]
        else:
            return [latent_code, mu, log_var]

    def encode(self, input: Tensor) -> List[Tensor]:
        p_latent_code = self.p_encode(input)[0] # get only the encoded input
        s_latent_code = self.s_encode(input)[0] # get only the encoded input
        return torch.cat((p_latent_code, s_latent_code), dim=1)
    
    def decode(self, z: Tensor) -> Tensor:
        outputs = self.decoder(z)
        return outputs

    def forward(self, input: Tensor, **kwargs) -> List[Tensor]:
        z = self.encode(input)
        mu_p = self.p_encode(input)[1] 
        log_var_p = self.p_encode(input)[2]
        mu_s = self.s_encode(input)[1]
        log_var_s = self.s_encode(input)[2]
        return [input, self.decode(z), z, mu_p, log_var_p, mu_s, log_var_s]

    def loss_function(self, *args, **kwargs) -> dict:
        input = args[0]
        recons = args[1]
        z = args[2]
        mu_p = args[3]
        log_var_p = args[4]
        mu_s = args[5]
        log_var_s = args[6]

        # reconstruction loss
        p_z = z[:, :z.shape[1] // 2]
        s_z = z[:, z.shape[1] // 2:]
        recons_loss = F.mse_loss(input, recons)

        # orthogonality loss
        s_l2_norm = torch.norm(s_z, p=2, dim=1, keepdim=True).detach()
        s_l2 = s_z.div(s_l2_norm.expand_as(s_z) + 1e-6)
        p_l2_norm = torch.norm(p_z, p=2, dim=1, keepdim=True).detach()
        p_l2 = p_z.div(p_l2_norm.expand_as(p_z) + 1e-6)
        ortho_loss = torch.mean((s_l2.t().mm(p_l2)).pow(2)) # frobenioius norm

        # MMD loss
        #true_samples_p=torch.randn(p_z.shape[0], p_z.shape[1]).to(self.device)
        #MMD_p=torch.sum(compute_mmd(true_samples_p, p_z))
        #true_samples_s=torch.randn(s_z.shape[0], s_z.shape[1]).to(self.device)
        #MMD_s=torch.sum(compute_mmd(true_samples_s, s_z))
        #MMD = MMD_p + MMD_s
        ## total loss
        # loss = recons_loss + self.alpha * ortho_loss + MMD
        
        # KL divergence
        kl_div_p = -0.5 * torch.sum(1 + log_var_p - mu_p.pow(2) - log_var_p.exp(), dim=1)
        kl_div_s = -0.5 * torch.sum(1 + log_var_s - mu_s.pow(2) - log_var_s.exp(), dim=1)
        #print(kl_div_p.mean()) # explodes to infinte, why?
        #print(kl_div_s.mean())
        kl_div = kl_div_p.mean() + kl_div_s.mean()

        # total loss
        loss = recons_loss + self.alpha * ortho_loss + self.beta * kl_div

        return {'loss': loss, 'recons_loss': recons_loss, 'ortho_loss': ortho_loss, 'kl_div': kl_div}

    # is this necessary? If not remove it
    def sample(self, num_samples: int, current_device: int, **kwargs) -> Tensor:
        z = torch.randn(num_samples, self.latent_dim)
        z = z.to(current_device)
        samples = self.decode(z)
        return samples

    def generate(self, x: Tensor, **kwargs) -> Tensor:
        return self.forward(x)[1]



