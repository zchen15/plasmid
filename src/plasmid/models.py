#!/usr/bin/env python
# Functions to 

import torch

class NeuralNetwork(torch.nn.Module):
    def __init__(self, args=None):
        super(NeuralNetwork, self).__init__()
        self.params = ftpbase.loadargs(args)

    def load_data(self):
        '''
        Load the data set
        '''
        fname = self.params['infile']

    def model1(self, x):
        '''
        Defines a basic dense network
        '''
        self.model = torch.nn.Sequential(
            torch.nn.Linear(28*28, 512),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(512),
            torch.nn.Linear(512, 512),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(512),
            torch.nn.Linear(512, 10),
            torch.nn.Sigmoid()
        )

    def set_loss_fn(self):
        self.loss_fn = torch.nn.KLDivLoss()

    def set_optimizer(self):
        self.optimizer = torch.optim.SGD()

    def forward(self, x):
        return self.model(x)

    def train(self, dataloader):
        size = len(dataloader.dataset)
        self.model.train()
        self.optimizer = self.optimizer(self.model.paramters(), lr=0.01)
        for batch, (X, y) in enumerate(dataloader):
            X, y = X.to(device), y.to(device)

            # Compute prediction error
            pred = self.model(X)
            loss = self.loss_fn(pred, y)

            # Backpropagation
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()

            if batch % 100 == 0:
                loss, current = loss.item(), batch * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

    def test(dataloader, model, loss_fn):
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        self.model.eval()
        test_loss, correct = 0, 0
        with torch.no_grad():
            for X, y in dataloader:
                X, y = X.to(device), y.to(device)
                pred = model(X)
                test_loss += self.loss_fn(pred, y).item()
                correct += (pred.argmax(1) == y).type(torch.float).sum().item()
        test_loss /= num_batches
        correct /= size
        print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

    def save_model(self):
        torch.save(model.state_dict(), "model.pth")
