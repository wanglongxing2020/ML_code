import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import pandas as pd
import random
def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
set_seed(1)
labels = pd.read_csv('./Entire_pd.csv')
labels = labels['RS']
labels = np.array(labels)
np.random.seed(42)
gene_expression_matrix = pd.read_csv('Entire_data.csv',index_col=0)
gene_expression_matrix = np.asarray(gene_expression_matrix).astype(np.float32)
gene_expression_tensor = torch.tensor(gene_expression_matrix)

class StackedAutoencoder(nn.Module):
    def __init__(self, input_dim, encoding_dims):
        super(StackedAutoencoder, self).__init__()

        # 构建编码器层
        encoder_layers = []
        in_dim = input_dim
        for encoding_dim in encoding_dims:
            encoder_layers.append(nn.Linear(in_dim, encoding_dim))
            encoder_layers.append(nn.ReLU(inplace=True))
            in_dim = encoding_dim

        self.encoder = nn.Sequential(*encoder_layers)

        # 构建解码器层
        decoder_layers = []
        for encoding_dim in reversed(encoding_dims):
            decoder_layers.append(nn.Linear(in_dim, encoding_dim))
            decoder_layers.append(nn.ReLU(inplace=True))
            in_dim = encoding_dim

        decoder_layers.append(nn.Linear(in_dim, input_dim))
        self.decoder = nn.Sequential(*decoder_layers)

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded
    
set_seed(22)
input_dim = gene_expression_matrix.shape[1]  # 输入维度，即基因表达矩阵的特征数
encoding_dims = [80, 60, 40, 20]
model = StackedAutoencoder(input_dim, encoding_dims)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 转换为PyTorch数据集和加载器
dataset = TensorDataset(gene_expression_tensor, gene_expression_tensor)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

num_epochs = 1000
for epoch in range(num_epochs):
    for data in dataloader:
        inputs, _ = data
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, inputs)
        loss.backward()
        optimizer.step()

    print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
encoder_model = model.encoder
encoded_features = encoder_model(gene_expression_tensor).detach().numpy()

set_seed(33)
lnc_expression_matrix = pd.read_csv('Entire_lncdata.csv',index_col=0)
lnc_expression_matrix = np.asarray(lnc_expression_matrix).astype(np.float32)
lnc_expression_tensor = torch.tensor(lnc_expression_matrix)
model = StackedAutoencoder(input_dim, encoding_dims)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)
dataset = TensorDataset(lnc_expression_tensor, lnc_expression_tensor)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
num_epochs = 1000
for epoch in range(num_epochs):
    for data in dataloader:
        inputs, _ = data
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, inputs)
        loss.backward()
        optimizer.step()

    print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
encoder_model = model.encoder
encoded_lnc = encoder_model(lnc_expression_tensor).detach().numpy()

set_seed(44)
meth_expression_matrix = pd.read_csv('Entire_methydata.csv',index_col=0)
meth_expression_matrix = np.asarray(meth_expression_matrix).astype(np.float32)
meth_expression_tensor = torch.tensor(meth_expression_matrix)
meth_expression_tensor.shape
model = StackedAutoencoder(input_dim, encoding_dims)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)
dataset = TensorDataset(meth_expression_tensor, meth_expression_tensor)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
num_epochs = 1000
for epoch in range(num_epochs):
    for data in dataloader:
        inputs, _ = data
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, inputs)
        loss.backward()
        optimizer.step()

    print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
encoder_model = model.encoder
encoded_meth = encoder_model(meth_expression_tensor).detach().numpy()

set_seed(55)
mi_expression_matrix = pd.read_csv('Entire_midata.csv',index_col=0)
mi_expression_matrix = np.asarray(mi_expression_matrix).astype(np.float32)
mi_expression_tensor = torch.tensor(mi_expression_matrix)
mi_expression_tensor.shape
model = StackedAutoencoder(input_dim, encoding_dims)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)
dataset = TensorDataset(mi_expression_tensor, mi_expression_tensor)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
num_epochs = 1000
for epoch in range(num_epochs):
    for data in dataloader:
        inputs, _ = data
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, inputs)
        loss.backward()
        optimizer.step()

    print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
encoder_model = model.encoder
encoded_mi = encoder_model(mi_expression_tensor).detach().numpy()


set_seed(6)
allre = np.column_stack((encoded_features,encoded_lnc,encoded_meth,encoded_mi))
allre = pd.DataFrame(allre)
allre.to_csv('allentire.csv')


set_seed(42)
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np

gene_expression_matrix = encoded_features
lnc_expression_matrix = encoded_lnc
meth_expression_matrix = encoded_meth
mi_expression_matrix = encoded_mi
#all_expression_matrix = allre
labels = pd.read_csv('Entire_pd.csv')
labels = labels['RS']
labels = np.array(labels)



# 转换为 PyTorch 张量
gene_expression_tensor = torch.tensor(gene_expression_matrix, dtype=torch.float32)
lnc_expression_tensor = torch.tensor(lnc_expression_matrix, dtype=torch.float32)
meth_expression_tensor = torch.tensor(meth_expression_matrix, dtype=torch.float32)
mi_expression_tensor = torch.tensor(mi_expression_matrix, dtype=torch.float32)
#all_expression_tensor =torch.tensor(all_expression_matrix, dtype=torch.float32)
labels_tensor = torch.tensor(labels, dtype=torch.float32).view(-1, 1)

# 划分数据集
X_train, X_test, y_train, y_test = train_test_split(gene_expression_tensor, labels_tensor, test_size=0.2, random_state=42)
X_train1, X_test1, y_train1, y_test1 = train_test_split(lnc_expression_tensor, labels_tensor, test_size=0.2, random_state=42)
X_train2, X_test2, y_train2, y_test2 = train_test_split(meth_expression_tensor, labels_tensor, test_size=0.2, random_state=42)
X_train3, X_test3, y_train3, y_test3 = train_test_split(mi_expression_tensor, labels_tensor, test_size=0.2, random_state=42)

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

set_seed(7)
# 初始化模型、损失函数和优化器
input_size = gene_expression_matrix.shape[1]
model = SimpleNN(input_size)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 转换为 PyTorch 数据加载器
train_dataloader = DataLoader(TensorDataset(X_train, y_train), batch_size=32, shuffle=True)
test_dataloader = DataLoader(TensorDataset(X_test, y_test), batch_size=32, shuffle=False)
all_dataloader = DataLoader(TensorDataset(gene_expression_tensor, labels_tensor), batch_size=32, shuffle=False)

train_dataloader1 = DataLoader(TensorDataset(X_train1, y_train1), batch_size=32, shuffle=True)
test_dataloader1 = DataLoader(TensorDataset(X_test1, y_test1), batch_size=32, shuffle=False)
all_dataloader1 = DataLoader(TensorDataset(lnc_expression_tensor, labels_tensor), batch_size=32, shuffle=False)

train_dataloader2 = DataLoader(TensorDataset(X_train2, y_train2), batch_size=32, shuffle=True)
test_dataloader2 = DataLoader(TensorDataset(X_test2, y_test2), batch_size=32, shuffle=False)
all_dataloader2 = DataLoader(TensorDataset(meth_expression_tensor, labels_tensor), batch_size=32, shuffle=False)

train_dataloader3 = DataLoader(TensorDataset(X_train3, y_train3), batch_size=32, shuffle=True)
test_dataloader3 = DataLoader(TensorDataset(X_test3, y_test3), batch_size=32, shuffle=False)
all_dataloader3 = DataLoader(TensorDataset(mi_expression_tensor, labels_tensor), batch_size=32, shuffle=False)

# 训练模型
set_seed(51234)
num_epochs = 1000
#mRNA
for epoch in range(num_epochs):
    model.train()
    for inputs, labels in train_dataloader:
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

# 在测试集上进行评估
model.eval()
all_preds, all_labels = [], []

with torch.no_grad():
    for inputs, labels in test_dataloader:
        preds = model(inputs)
        all_preds.extend(preds.numpy())
        all_labels.extend(labels.numpy())
#全数据集得到完整的preds
model.eval()
all_preds, all_labels = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader:
        preds = model(inputs)
        all_preds.extend(preds.numpy())
        all_labels.extend(labels.numpy())
# 计算 ROC 曲线
fpr, tpr, _ = roc_curve(all_labels, all_preds)
roc_auc = auc(fpr, tpr)

#lncRNA
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
#全数据集得到完整的preds
model.eval()
all_preds1, all_labels1 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader1:
        preds = model(inputs)
        all_preds1.extend(preds.numpy())
        all_labels1.extend(labels.numpy())
# 计算 ROC 曲线
fpr1, tpr1, _ = roc_curve(all_labels1, all_preds1)
roc_auc1 = auc(fpr1, tpr1)

#Methylation
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
#全数据集得到完整的preds
model.eval()
all_preds2, all_labels2 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader2:
        preds = model(inputs)
        all_preds2.extend(preds.numpy())
        all_labels2.extend(labels.numpy())
# 计算 ROC 曲线
fpr2, tpr2, _ = roc_curve(all_labels2, all_preds2)
roc_auc2 = auc(fpr2, tpr2)

#miRNA
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
#全数据集得到完整的preds
model.eval()
all_preds3, all_labels3 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader3:
        preds = model(inputs)
        all_preds3.extend(preds.numpy())
        all_labels3.extend(labels.numpy())
# 计算 ROC 曲线
fpr3, tpr3, _ = roc_curve(all_labels3, all_preds3)
roc_auc3 = auc(fpr3, tpr3)


# 画出 ROC 曲线
plt.figure(figsize=(8,8))
plt.plot(fpr, tpr, color='r', label=f'mRNA ROC Curve (AUC = {roc_auc:.3f})')
plt.plot(fpr3, tpr3, color='y', label=f'miRNA ROC Curve (AUC = {roc_auc3:.3f})')
plt.plot(fpr1, tpr1, color='b', label=f'lncRNA ROC Curve (AUC = {roc_auc1:.3f})')
plt.plot(fpr2, tpr2, color='g', label=f'Methylation ROC Curve (AUC = {roc_auc2:.3f})')
plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend()

plt.savefig('Entire.pdf')
