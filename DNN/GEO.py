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
#GSE14520
labels = pd.read_csv('./GSE14520_pd.csv')
labels = labels['RS']
labels = np.array(labels)

gene_expression_matrix = pd.read_csv('GSE14520_data.csv',index_col=0)
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
    
set_seed(2)
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

set_seed(3)
#GSE116174
labels = pd.read_csv('./GSE116174_pd.csv')
labels = labels['RS']
labels = np.array(labels)

gene_expression_matrix = pd.read_csv('GSE116174_data.csv',index_col=0)
gene_expression_matrix = np.asarray(gene_expression_matrix).astype(np.float32)
gene_expression_tensor = torch.tensor(gene_expression_matrix)
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
encoded_features2 = encoder_model(gene_expression_tensor).detach().numpy()

set_seed(4)
#ICGC
labels = pd.read_csv('./ICGC_pd.csv')
labels = labels['RS']
labels = np.array(labels)

gene_expression_matrix = pd.read_csv('ICGC_data.csv',index_col=0)
gene_expression_matrix = np.asarray(gene_expression_matrix).astype(np.float32)
gene_expression_tensor = torch.tensor(gene_expression_matrix)
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
encoded_features3 = encoder_model(gene_expression_tensor).detach().numpy()


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
#省去改名，alltrain=GSE14520;alltest=GSE116174;allentire=ICGC
alltrain = encoded_features
train_pd = pd.read_csv('GSE14520_pd.csv')
train_pd = train_pd['RS']
train_pd = np.array(train_pd)

alltest = encoded_features2
test_pd = pd.read_csv('GSE116174_pd.csv')
test_pd = test_pd['RS']
test_pd = np.array(test_pd)

allentire = encoded_features3
entire_pd = pd.read_csv('ICGC_pd.csv')
entire_pd = entire_pd['RS']
entire_pd = np.array(entire_pd)


alltrain_tensor = torch.tensor(alltrain, dtype=torch.float32)
alltest_tensor = torch.tensor(alltest, dtype=torch.float32)
entire_tensor = torch.tensor(allentire, dtype=torch.float32)

labels_tensor1 = torch.tensor(train_pd, dtype=torch.float32).view(-1, 1)
labels_tensor2 = torch.tensor(test_pd, dtype=torch.float32).view(-1, 1)
labels_tensor3 = torch.tensor(entire_pd, dtype=torch.float32).view(-1, 1)

# 划分数据集

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


set_seed(5)
# 初始化模型、损失函数和优化器
input_size = alltrain.shape[1]
model = SimpleNN(input_size)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)


train_dataloader1 = DataLoader(TensorDataset(X_train1, y_train1), batch_size=32, shuffle=True)
test_dataloader1 = DataLoader(TensorDataset(X_test1, y_test1), batch_size=32, shuffle=False)

#全部的
all_dataloader1 = DataLoader(TensorDataset(alltrain_tensor, labels_tensor1), batch_size=32, shuffle=False)

train_dataloader2 = DataLoader(TensorDataset(X_train2, y_train2), batch_size=32, shuffle=True)
test_dataloader2 = DataLoader(TensorDataset(X_test2, y_test2), batch_size=32, shuffle=False)
all_dataloader2 = DataLoader(TensorDataset(alltest_tensor, labels_tensor2), batch_size=32, shuffle=False)

train_dataloader3 = DataLoader(TensorDataset(X_train3, y_train3), batch_size=32, shuffle=True)
test_dataloader3 = DataLoader(TensorDataset(X_test3, y_test3), batch_size=32, shuffle=False)
all_dataloader3 = DataLoader(TensorDataset(entire_tensor, labels_tensor3), batch_size=32, shuffle=False)


set_seed(123)
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

# 在全数据集上为得到完整preds
model.eval()
all_preds1, all_labels1 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader1:
        preds = model(inputs)
        all_preds1.extend(preds.numpy())
        all_labels1.extend(labels.numpy())
x = np.vstack(all_preds1)
x = pd.DataFrame(x)
x.to_csv('GSE14520preds.csv')        
        
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
# 在全数据集上为得到完整preds
model.eval()
all_preds2, all_labels2 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader2:
        preds = model(inputs)
        all_preds2.extend(preds.numpy())
        all_labels2.extend(labels.numpy())
x = np.vstack(all_preds2)
x = pd.DataFrame(x)
x.to_csv('GSE116174preds.csv') 

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
# 在全数据集上为得到完整preds
model.eval()
all_preds3, all_labels3 = [], []

with torch.no_grad():
    for inputs, labels in all_dataloader3:
        preds = model(inputs)
        all_preds3.extend(preds.numpy())
        all_labels3.extend(labels.numpy())
x = np.vstack(all_preds3)
x = pd.DataFrame(x)
x.to_csv('ICGCpreds.csv') 

# 计算 ROC 曲线
fpr3, tpr3, _ = roc_curve(all_labels3, all_preds3)
roc_auc3 = auc(fpr3, tpr3)



# 画出 ROC 曲线
plt.figure(figsize=(8,8))

plt.plot(fpr1, tpr1, color='r', label=f'GSE14520 (AUC = {roc_auc1:.3f})')
plt.plot(fpr2, tpr2, color='y', label=f'GSE116174 (AUC = {roc_auc2:.3f})')
plt.plot(fpr3, tpr3, color='b', label=f'ICGC (AUC = {roc_auc3:.3f})')
plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend()

plt.savefig('GEO.pdf')





