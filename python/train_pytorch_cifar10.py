import argparse
import os
import torch
from torch import nn, optim
from torch.nn import functional as F
from torch.utils.data import DataLoader
from torchvision import datasets, transforms





# arguments to parse
parser = argparse.ArgumentParser(description='PyTorch CIFAR10')
parser.add_argument('--batch-size', type=int, default=64, metavar='N',
                    help='input batch size for training (default: 64)')
parser.add_argument('--epochs', type=int, default=10, metavar='N',
                    help='number of epochs to train (default: 10)')
parser.add_argument('--lr', type=float, default=0.01, metavar='LR',
                    help='learning rate (default: 0.01)')
parser.add_argument('--momentum', type=float, default=0.5, metavar='M',
                    help='SGD momentum (default: 0.5)')
parser.add_argument('--no-cuda', action='store_true', default=False,
                    help='disables CUDA training')
parser.add_argument('--seed', type=int, default=1, metavar='S',
                    help='random seed (default: 1)')
parser.add_argument('--save-interval', type=int, default=10, metavar='N',
                    help='how many batches to wait before checkpointing')
parser.add_argument('--resume', action='store_true', default=False,
                    help='resume training from checkpoint')
args = parser.parse_args()


# using cuda or not
use_cuda = torch.cuda.is_available() and not args.no_cuda
device = torch.device('cuda' if use_cuda else 'cpu')
torch.manual_seed(args.seed)
if use_cuda:
  torch.cuda.manual_seed(args.seed)



# loading data 
transform = transforms.Compose(
    [transforms.ToTensor(),
     transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])

trainset = torchvision.datasets.CIFAR10(root='./data', train=True,
                                        download=True, transform=transform)
train_loader = torch.utils.data.DataLoader(trainset, batch_size=64,
                                          shuffle=True, num_workers=2)

testset = torchvision.datasets.CIFAR10(root='./data', train=False,
                                       download=True, transform=transform)
test_loader = torch.utils.data.DataLoader(testset, batch_size=10000,
                                         shuffle=False, num_workers=2)

classes = ('plane', 'car', 'bird', 'cat',
           'deer', 'dog', 'frog', 'horse', 'ship', 'truck')


# network
class Net(nn.Module):
  def __init__(self):
    super(Net, self).__init__()
    self.conv1 = nn.Conv2d(1, 10, kernel_size=5)
    self.conv2 = nn.Conv2d(10, 20, kernel_size=5)
    self.conv2_drop = nn.Dropout2d()
    self.fc1 = nn.Linear(320, 50)
    self.fc2 = nn.Linear(50, 10)

  def forward(self, x):
    x = F.relu(F.max_pool2d(self.conv1(x), 2))
    x = F.relu(F.max_pool2d(self.conv2_drop(self.conv2(x)), 2))
    x = x.view(-1, 320)
    x = F.relu(self.fc1(x))
    x = self.fc2(x)
    return F.log_softmax(x, dim=1)

model = Net().to(device)
optimiser = optim.SGD(model.parameters(), lr=args.lr, momentum=args.momentum)

if args.resume:
  model.load_state_dict(torch.load('model.pth'))
  optimiser.load_state_dict(torch.load('optimiser.pth'))

# training

model.train()
train_losses = []

for i, (data, target) in enumerate(train_loader):
  data = data.to(device=device, non_blocking=True)
  target = target.to(device=device, non_blocking=True)
  optimiser.zero_grad()
  output = model(data)
  loss = F.nll_loss(output, target)
  loss.backward()
  train_losses.append(loss.item())
  optimiser.step()

  if i % 10 == 0:
    print(i, loss.item())
    torch.save(model.state_dict(), 'model.pth')
    torch.save(optimiser.state_dict(), 'optimiser.pth')
    torch.save(train_losses, 'train_losses.pth')


# testing
model.eval()
test_loss, correct = 0, 0

with torch.no_grad():
  for data, target in test_loader:
    data = data.to(device=device, non_blocking=True)
    target = target.to(device=device, non_blocking=True)
    output = model(data)
    test_loss += F.nll_loss(output, target, reduction='sum').item()
    pred = output.argmax(1, keepdim=True)
    correct += pred.eq(target.view_as(pred)).sum().item()

test_loss /= len(test_data)
acc = correct / len(test_data)
print(acc, test_loss)
