import sys, os
import pandas as pd
import torch

def model_save_check(history, metric_name, tolerance_count=5, reset_count=1):
    """
    The function `model_save_check` checks if a given metric value in a history dictionary is better
    than the previous best value, and determines whether to save the model or stop training based on a
    tolerance count and reset count.
    
    :param history: The `history` parameter is a dictionary that contains the training history of a
    model. It typically includes the values of different metrics (e.g., loss, accuracy) at each epoch
    during training
    :param metric_name: The metric_name parameter is the name of the metric that you want to use for
    saving and stopping the model. It could be any metric that you are tracking during the training
    process, such as accuracy, loss, precision, recall, etc
    :param tolerance_count: The tolerance_count parameter is the number of consecutive epochs where the
    metric does not improve before considering stopping the training process, defaults to 5 (optional)
    :param reset_count: The `reset_count` parameter is used to determine how many times the
    `tolerance_count` should be exceeded before stopping the training process. If the difference between
    the current index and the best index is greater than `tolerance_count * reset_count`, the
    `stop_flag` will be set, defaults to 1 (optional)
    :return: two values: `save_flag` and `stop_flag`.
    """
    save_flag = False
    stop_flag = False
    if 'best_index' not in history:
        history['best_index'] = 0
    if metric_name.endswith('loss'):
        if history[metric_name][-1] <= history[metric_name][history['best_index']]:
            save_flag = True
            history['best_index'] = len(history[metric_name]) - 1
    else:
        if history[metric_name][-1] >= history[metric_name][history['best_index']]:
            save_flag = True
            history['best_index'] = len(history[metric_name]) - 1

    if len(history[metric_name]) - history['best_index'] > tolerance_count * reset_count and history['best_index'] > 0:
        stop_flag = True

    return save_flag, stop_flag

# different from the original code
def generate_encoded_features(encoder, dataloader, normalize_flag=False):
    encoder.eval()
    raw_feature_tensor = dataloader[0].dataset.tensors[0].cpu()
    sample_info = dataloader[2]
    with torch.no_grad():
        encoded_feature_tensor = encoder.cpu()(raw_feature_tensor)
    # in the training, we use the normalization, but in the evaluation, we don't. Correct?
    if normalize_flag:
        encoded_feature_tensor = torch.nn.functional.normalize(encoded_feature_tensor, p=2, dim=1)
    encoded_feature_tensor = encoded_feature_tensor.cpu().detach().numpy()
    enc_out = pd.DataFrame(encoded_feature_tensor, index = sample_info.index)
    return enc_out

