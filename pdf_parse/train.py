import matplotlib.pyplot as plt
import numpy as np
import torch as tc
from torch import nn
from torch import optim
import torch.nn.functional as F
import torchvision as tv
from torch.optim import Adam, lr_scheduler
import time
import os
import copy

import model as m
import set_train_data as std

device = tc.device("cpu")

# Load Training Data & Validation Data for split
seed = 1
b_size = 100
v_size = 0.2
data_path='/users/spyder/Desktop/training/'

train_data, valid_data = std.load_dataset(data_path=data_path,
                                      random_seed = seed, 
                                      batch_size = b_size,
                                      valid_size = v_size
        )

data_loader = {"train": train_data, "val":valid_data}

# Hyper Parameters for tuning
lr = 0.001
w_decay = 0.0001
n_class = 7

num_epoch = 3

# 
model = m.ParseNet(num_classes=7)

optimizer = Adam(model.parameters(), lr=lr, weight_decay= w_decay)
loss_fn = nn.CrossEntropyLoss()

print("Data Finished Loading -------")


def save_models(epoch):
    tc.save(model.state_dict(), "parsenet_{}.model".format(epoch))
    print("Checkpoint Saved")

def train(model, criterion, optimizer, num_epochs, scheduler):
    since = time.time()

    best_model_wts = copy.deepcopy(model.state_dict())
    best_acc = 0.0

    for epoch in range(num_epoch):
        print('Epoch {}/{}'.format(epoch, num_epochs - 1))
        print('-' * 10)

        # Each epoch has a training and validation phase
        for phase in ['train', 'val']:
            if phase == 'train':
                scheduler.step()
                model.train()  # Set model to training mode
            else:
                model.eval()   # Set model to evaluate mode

            running_loss = 0.0
            running_corrects = 0

            # Iterate over data.
            for inputs, labels in data_loader[phase]:
                inputs = inputs.to(device)
                labels = labels.to(device)

                # zero the parameter gradients
                optimizer.zero_grad()

                # forward
                # track history if only in train
                with tc.set_grad_enabled(phase == 'train'):
                    outputs = model(inputs)
                    _, preds = tc.max(outputs, 1)
                    loss = criterion(outputs, labels)

                    # backward + optimize only if in training phase
                    if phase == 'train':
                        loss.backward()
                        optimizer.step()

                # statistics
                running_loss += loss.item() * inputs.size(0)
                running_corrects += tc.sum(preds == labels.data)

            epoch_loss = running_loss / dataset_sizes[phase]
            epoch_acc = running_corrects.double() / dataset_sizes[phase]

            print('{} Loss: {:.4f} Acc: {:.4f}'.format(
                phase, epoch_loss, epoch_acc))

            # deep copy the model
            if phase == 'val' and epoch_acc > best_acc:
                best_acc = epoch_acc
                best_model_wts = copy.deepcopy(model.state_dict())
    

exp_lr_scheduler = lr_scheduler.StepLR(optimizer, step_size=7, gamma=0.1)
train(model, loss_fn, optimizer, num_epoch, exp_lr_scheduler)
