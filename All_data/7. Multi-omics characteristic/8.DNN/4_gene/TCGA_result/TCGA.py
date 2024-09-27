import torch
import pandas as pd
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np
import random
def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
set_seed(1)

alltrain = pd.read_csv('alltrain.csv',index_col=0)
alltrain = np.asarray(alltrain)
train_pd = pd.read_csv('train_pd.csv')
train_pd = train_pd['RS']
train_pd = np.array(train_pd)

alltest = pd.read_csv('alltest.csv',index_col=0)
alltest = np.asarray(alltest)
test_pd = pd.read_csv('test_pd.csv')
test_pd = test_pd['RS']
test_pd = np.array(test_pd)

allentire = pd.read_csv('allentire.csv',index_col=0)
allentire = np.asarray(allentire)
entire_pd = pd.read_csv('entire_pd.csv')
entire_pd = entire_pd['RS']
entire_pd = np.array(entire_pd)


alltrain_tensor = torch.tensor(alltrain, dtype=torch.float32)
alltest_tensor = torch.tensor(alltest, dtype=torch.float32)
entire_tensor = torch.tensor(allentire, dtype=torch.float32)

labels_tensor1 = torch.tensor(train_pd, dtype=torch.float32).view(-1, 1)
labels_tensor2 = torch.tensor(test_pd, dtype=torch.float32).view(-1, 1)
labels_tensor3 = torch.tensor(entire_pd, dtype=torch.float32).view(-1, 1)

# 划分数据集
set_seed(42)

X_train1, X_test1, y_train1, y_test1 = train_test_split(alltrain_tensor, labels_tensor1, test_size=0.2, random_state=42)
X_train2, X_test2, y_train2, y_test2 = train_test_split(alltest_tensor, labels_tensor2, test_size=0.2, random_state=42)
X_train3, X_test3, y_train3, y_test3 = train_test_split(entire_tensor, labels_tensor3, test_size=0.2, random_state=42)

# 定义深度神经网络模型
class SimpleNN(nn.Module):
    def __init__(self, input_size):
        super(SimpleNN, self).__init__()
        self.fc1 = nn.Linear(input_size, 256)
        self.relu1 = nn.ReLU()
        self.fc2 = nn.Linear(256, 128)
        self.relu2 = nn.ReLU()
        self.fc3 = nn.Linear(128, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.relu1(self.fc1(x))
        x = self.relu2(self.fc2(x))
        x = self.sigmoid(self.fc3(x))
        return x

set_seed(2)
# 初始化模型、损失函数和优化器
input_size = alltrain.shape[1]
model = SimpleNN(input_size)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)


train_dataloader1 = DataLoader(TensorDataset(X_train1, y_train1), batch_size=32, shuffle=True)
test_dataloader1 = DataLoader(TensorDataset(X_test1, y_test1), batch_size=32, shuffle=False)
all_dataloader1 = DataLoader(TensorDataset(alltrain_tensor, labels_tensor1), batch_size=32, shuffle=False)



train_dataloader2 = DataLoader(TensorDataset(X_train2, y_train2), batch_size=32, shuffle=True)
test_dataloader2 = DataLoader(TensorDataset(X_test2, y_test2), batch_size=32, shuffle=False)
all_dataloader2 = DataLoader(TensorDataset(alltest_tensor, labels_tensor2), batch_size=32, shuffle=False)

train_dataloader3 = DataLoader(TensorDataset(X_train3, y_train3), batch_size=32, shuffle=True)
test_dataloader3 = DataLoader(TensorDataset(X_test3, y_test3), batch_size=32, shuffle=False)
all_dataloader3 = DataLoader(TensorDataset(entire_tensor, labels_tensor3), batch_size=32, shuffle=False)

set_seed(520123)
num_epochs=1000
for epoch in range(num_epochs):
    model.train()
    for inputs, labels in train_dataloader1:
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

# 在测试集上进行评估
model.eval()
all_preds1, all_labels1 = [], []

with torch.no_grad():
    for inputs, labels in test_dataloader1:
        preds = model(inputs)
        all_preds1.extend(preds.numpy())
        all_labels1.extend(labels.numpy())
# 在全数据集上为得到完整preds 单独跑
model.eval()
all_preds1, all_labels1 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader1:
        preds = model(inputs)
        all_preds1.extend(preds.numpy())
        all_labels1.extend(labels.numpy())
x = np.vstack(all_preds1)
x = pd.DataFrame(x)
x.to_csv('Trainpreds.csv')  

# 计算 ROC 曲线
fpr1, tpr1, _ = roc_curve(all_labels1, all_preds1)
roc_auc1 = auc(fpr1, tpr1)

for epoch in range(num_epochs):
    model.train()
    for inputs, labels in train_dataloader2:
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

# 在测试集上进行评估
model.eval()
all_preds2, all_labels2 = [], []

with torch.no_grad():
    for inputs, labels in test_dataloader2:
        preds = model(inputs)
        all_preds2.extend(preds.numpy())
        all_labels2.extend(labels.numpy())
# 在全数据集上为得到完整preds 单独跑
model.eval()
all_preds2, all_labels2 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader2:
        preds = model(inputs)
        all_preds2.extend(preds.numpy())
        all_labels2.extend(labels.numpy())
x = np.vstack(all_preds2)
x = pd.DataFrame(x)
x.to_csv('Testpreds.csv') 
# 计算 ROC 曲线
fpr2, tpr2, _ = roc_curve(all_labels2, all_preds2)
roc_auc2 = auc(fpr2, tpr2)

for epoch in range(num_epochs):
    model.train()
    for inputs, labels in train_dataloader3:
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

# 在测试集上进行评估
model.eval()
all_preds3, all_labels3 = [], []

with torch.no_grad():
    for inputs, labels in test_dataloader3:
        preds = model(inputs)
        all_preds3.extend(preds.numpy())
        all_labels3.extend(labels.numpy())
# 在全数据集上为得到完整preds 单独跑
model.eval()
all_preds3, all_labels3 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader3:
        preds = model(inputs)
        all_preds3.extend(preds.numpy())
        all_labels3.extend(labels.numpy())
x = np.vstack(all_preds3)
x = pd.DataFrame(x)
x.to_csv('Entirepreds.csv') 
# 计算 ROC 曲线
fpr3, tpr3, _ = roc_curve(all_labels3, all_preds3)
roc_auc3 = auc(fpr3, tpr3)


# 画出 ROC 曲线
plt.figure(figsize=(8,8))

plt.plot(fpr1, tpr1, color='r', label=f'Train ROC Curve (AUC = {roc_auc1:.3f})')
plt.plot(fpr2, tpr2, color='y', label=f'Test ROC Curve (AUC = {roc_auc2:.3f})')
plt.plot(fpr3, tpr3, color='b', label=f'Entire ROC Curve (AUC = {roc_auc3:.3f})')
plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend()

plt.savefig('TCGA.pdf')