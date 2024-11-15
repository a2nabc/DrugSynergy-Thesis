import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split

# Load dataset
file_path = 'Data/OmicsExpressionGenesExpectedCountProfile.csv'
data = pd.read_csv(file_path)

# Inspect the data
print(data.head())
 
# Drop 1st column with unnamed header 
data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
print(data.head())

# Normalize gene expression data
data = data.apply(lambda x: np.log1p(x))
scaler = MinMaxScaler()
data_normalized = scaler.fit_transform(data)

# TO DO: DEFINE RIGHT TARGET LABELS!!!!!
# Create target labels (binary classification) (SIMULATES SYNERGY DATA)
targets = np.random.randint(0, 2, size=data.shape[0])

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(data_normalized, targets, test_size=0.2, random_state=42)

################################################### MODEL ###################################################

# Define a simple feedforward neural network
class SimpleNN(nn.Module):
    def __init__(self, input_size):
        super(SimpleNN, self).__init__()
        self.fc1 = nn.Linear(input_size, 128)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(128, 1)
        self.sigmoid = nn.Sigmoid()
    
    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.sigmoid(x)
        return x

# Instantiate model
input_size = X_train.shape[1]
model = SimpleNN(input_size)

# Define loss and optimizer
criterion = nn.BCELoss()  # Binary Cross-Entropy Loss
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Convert data to PyTorch tensors
X_train_tensor = torch.FloatTensor(X_train)
y_train_tensor = torch.FloatTensor(y_train).view(-1, 1)
X_test_tensor = torch.FloatTensor(X_test)
y_test_tensor = torch.FloatTensor(y_test).view(-1, 1)

# Training loop
epochs = 20
for epoch in range(epochs):
    model.train()
    optimizer.zero_grad()
    
    # Forward pass
    outputs = model(X_train_tensor)
    loss = criterion(outputs, y_train_tensor)
    
    # Backward pass
    loss.backward()
    optimizer.step()
    
    print(f'Epoch {epoch+1}/{epochs}, Loss: {loss.item()}')

# Evaluate on test set
model.eval()
with torch.no_grad():
    test_outputs = model(X_test_tensor)
    predictions = (test_outputs > 0.5).float()
