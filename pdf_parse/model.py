import torch
import torch.nn as nn

class ParseNet(nn.Module):
    def __init__(self, num_classes:int=7):
        super(ParseNet, self).__init__()

        self.conv1 = nn.Conv2d(in_channels = 3, out_channels=9, kernel_size=3,
                               stride=1, padding=1);
        self.relu1 = nn.ReLU()

        self.pool = nn.MaxPool2d(kernel_size=2, stride=2)

        self.conv2 = nn.Conv2d(in_channels=9, out_channels=9, kernel_size=3,
                               stride=1, padding=1)
        self.relu2 = nn.ReLU()

        self.fc = nn.Linear(in_features=15 * 15 * 9,  out_features = num_classes)
    
    def forward(self, input):
        output = self.conv1(input)
        output = self.relu1(output)

        output = self.pool(output)

        output = self.conv2(output)
        output = self.relu2(output)
        output = output.view(output.size(0), -1)
        
        output = self.fc(output)
        return output
