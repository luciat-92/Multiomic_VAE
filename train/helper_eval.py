import pandas as pd
from collections import defaultdict
import numpy as np
import torch
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, f1_score, \
    log_loss, auc, precision_recall_curve


def auprc(y_true, y_score):
    """
    The function calculates the area under the precision-recall curve for a binary classification
    problem.
    
    :param y_true: The true labels or ground truth values for the binary classification problem. It is a
    1-dimensional array-like object (e.g., list, numpy array) containing the true labels for each sample
    in the dataset. The labels should be binary (0 or 1) indicating the positive and negative classes
    :param y_score: The predicted scores or probabilities for the positive class. These scores indicate
    the likelihood of an instance belonging to the positive class
    :return: the area under the precision-recall curve (AUPRC) for the given true labels and predicted
    scores.
    """
    lr_precision, lr_recall, _ = precision_recall_curve(y_true=y_true, probas_pred=y_score)
    return auc(lr_recall, lr_precision)

def compute_kernel(x, y):
    x_size = x.size(0)
    y_size = y.size(0)
    dim = x.size(1)
    x = x.unsqueeze(1) # (x_size, 1, dim)
    y = y.unsqueeze(0) # (1, y_size, dim)
    tiled_x = x.expand(x_size, y_size, dim)
    tiled_y = y.expand(x_size, y_size, dim)
    kernel_input = (tiled_x - tiled_y).pow(2).mean(2)/float(dim)
    return torch.exp(-kernel_input) # (x_size, y_size)

def compute_mmd(x, y):
    x_kernel = compute_kernel(x, x)
    y_kernel = compute_kernel(y, y)
    xy_kernel = compute_kernel(x, y)
    mmd = x_kernel.mean() + y_kernel.mean() - 2*xy_kernel.mean()
    return mmd

def eval_dsn_epoch(model, data_loader, device, history):
    """
    The function `eval_dsnae_epoch` evaluates the performance of a model on a given dataset for one
    epoch and updates the loss history.
    
    :param model: The model is the neural network model that you want to evaluate. It should be an
    instance of a PyTorch model class
    :param data_loader: The `data_loader` parameter is an iterator that provides batches of data to the
    model for evaluation. It is typically created using the `torch.utils.data.DataLoader` class and is
    used to load data in parallel during training or evaluation
    :param device: The `device` parameter is used to specify whether the model should be run on a CPU or
    a GPU. It is a string that can take values like "cpu" or "cuda:0" (if you have a GPU). This
    parameter is used to move the input data to the specified device
    :param history: The `history` parameter is a dictionary that stores the average loss values for each
    epoch. It is updated with the loss values calculated for each batch in the `data_loader`. The keys
    in the `history` dictionary represent different loss components, and the values represent the
    average loss values for each component
    :return: the updated history dictionary.
    """
    model.eval()
    avg_loss_dict = defaultdict(float)
    for x_batch in data_loader:
        x_batch = x_batch[0].to(device)
        with torch.no_grad():
            loss_dict = model.loss_function(*(model(x_batch)))
            for k, v in loss_dict.items():
                avg_loss_dict[k] += v.cpu().detach().item() / len(data_loader)

    for k, v in avg_loss_dict.items():
        history[k].append(v)
    return history


