import matplotlib.pyplot as plt
import numpy as np
import torch as tc
from torch import nn
from torch import optim
import torch.nn.functional as F
import torchvision as tv

# List Labels
label_names = [
    'bl',
    'da',
    'dt',
    'oc',
    'sq',
    'ut',
    'xx',
]


# Helper utility
def plot_images(images, cls_true, cls_pred=None):
    """
    Adapted from https://github.com/Hvass-Labs/TensorFlow-Tutorials/
    """
    fig, axes = plt.subplots(3, 3)

    for i, ax in enumerate(axes.flat):
        # plot img
        ax.imshow(images[i, :, :, :], interpolation='spline16')

        # show true & predicted classes
        cls_true_name = label_names[cls_true[i]]
        if cls_pred is None:
            xlabel = "{0} ({1})".format(cls_true_name, cls_true[i])
        else:
            cls_pred_name = label_names[cls_pred[i]]
            xlabel = "True: {0}\nPred: {1}".format(
                cls_true_name, cls_pred_name
            )
        ax.set_xlabel(xlabel)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.show()


# Develop training_data sets and batches

def load_dataset(data_path='/users/spyder/Desktop/training/',
                 batch_size=64,
                 valid_size=0.2,
                 random_seed = 0,
                 shuffle=True,
                 show_sample=False,
                 num_workers=0,
                 max_num = None
        ):    

    # Normalize Data
    normalize = tv.transforms.Normalize(
            mean=[0.4914, 0.4822, 0.4465],
            std=[0.2023, 0.1994, 0.2010],
        )

    valid_transform = tv.transforms.Compose([
            tv.transforms.Grayscale(num_output_channels=3),
            tv.transforms.ToTensor(),
            normalize
            ])
    train_transform = tv.transforms.Compose([
            tv.transforms.Grayscale(num_output_channels=3),
            tv.transforms.ToTensor(),
            normalize,
        ])

    # load the dataset
    train_dataset = tv.datasets.ImageFolder(
                root=data_path, 
                transform=train_transform
            )
    
    valid_dataset = tv.datasets.ImageFolder(
                root=data_path,
                transform=valid_transform
            )

    # Shuffle indices of data
    num_train = len(train_dataset)
    indices = list(range(num_train))
    split = int(np.floor(valid_size*num_train))

    if shuffle:
        np.random.seed(random_seed)
        np.random.shuffle(indices)
    
    train_idx, valid_idx = indices[split:], indices[:split]
    if(max_num):
        train_sampler = tc.utils.data.sampler.SubsetRandomSampler(train_idx[:max_num])
        valid_sampler = tc.utils.data.sampler.SubsetRandomSampler(valid_idx[:int(max_num*valid_size)])
    else:
        train_sampler = tc.utils.data.sampler.SubsetRandomSampler(train_idx)
        valid_sampler = tc.utils.data.sampler.SubsetRandomSampler(valid_idx)

    # Load data
    train_loader = tc.utils.data.DataLoader(
                train_dataset,
                batch_size=batch_size,
                num_workers=num_workers,
                sampler=train_sampler
        )

    valid_loader = tc.utils.data.DataLoader(
                valid_dataset,
                batch_size=batch_size,
                num_workers=num_workers,
                sampler=valid_sampler
        )
    
    train_size = len(train_dataset)
    valid_size = len(valid_dataset)

    if show_sample:
        sample_loader = tc.utils.data.DataLoader(
            train_dataset, batch_size=9, shuffle=shuffle,
            num_workers=num_workers
        )
        data_iter = iter(sample_loader)
        images, labels = data_iter.next()
        X = images.numpy().transpose([0, 2, 3, 1])
        plot_images(X, labels)


    return train_loader, valid_loader, train_size, valid_size



    





