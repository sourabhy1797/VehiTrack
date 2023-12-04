import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
import scipy.io
import mat73


# Reading Node details
filename = 'new_nodes.mat'
path = '../Posteriors Dataset/'+str(filename)
mat_file = scipy.io.loadmat(path)
cell_array = mat_file['new_node_ids']
nodes_list = []
for i in range(len(cell_array)):
    nodes = cell_array[i][0]
    nodes_list.append(nodes)

num_nodes = len(nodes_list)

# Reading the trajectories
mat_file = scipy.io.loadmat('../Posteriors Dataset/new_trajectories.mat')

cell_array = mat_file['new_traj']
trajectories = []
for i in range(10000):
    traj = cell_array[i][0][0].tolist()
    trajectories.append(traj)
    
print("Trajectory File Readed")


node_to_index = {node: index for index, node in enumerate(nodes_list)}
new_trajectories = []
for trajectory in trajectories:
    new_traj = [node_to_index[node] for node in trajectory]
    new_trajectories.append(new_traj)

num_trajectories = len(new_trajectories)
max_trajectory_length = max(max(new_trajectories, key=len))

# One-hot encode the synthetic trajectories
one_hot_trajectories = []

for trajectory in new_trajectories:
    one_hot_trajectory = []
    for node_id in trajectory:
        one_hot_vector = [0] * num_nodes
        one_hot_vector[node_id - 1] = 1  # Convert node_id to 0-based indexing
        one_hot_trajectory.append(one_hot_vector)
    one_hot_trajectories.append(one_hot_trajectory)


k=1
# Reading Posterior Vectors
filename = 'pv'+str(k)+'.mat'
path = '../Posteriors Dataset/'+str(filename)
mat_file = mat73.loadmat(path)
array_name = 'px'+str(k)
cell_array = mat_file[array_name]
pv_data = []
post_vectors = []
for m in range(10000):
    element = cell_array[m]
    pv = []
    for j in range(len(element[0])):
        vec = element[0][j]
        pv.append((vec).tolist())
    post_vectors.append(pv)
print('File '+filename+' is readed')

posterior_data = post_vectors


# Pad the sequences to the maximum trajectory length

one_hot_trajectories = []

for trajectory in new_trajectories:
    one_hot_trajectory = []
    for node_id in trajectory:
        one_hot_vector = np.zeros(num_nodes, dtype=int)
        one_hot_vector[node_id] = 1
        one_hot_trajectory.append(one_hot_vector)
    one_hot_trajectories.append(one_hot_trajectory)

print("On hot encoding done")

trajectory_lengths = [len(trajectory) for trajectory in post_vectors]
max_length = max(trajectory_lengths)
padded_trajectories = torch.zeros((num_trajectories, max_length, num_nodes))
for i, trajectory in enumerate(post_vectors):
    padded_trajectories[i, :len(trajectory), :] = torch.Tensor(trajectory)

padded_trajectories_ids = torch.zeros((num_trajectories, max_length, num_nodes))
for i, trajectory in enumerate(one_hot_trajectories):
    padded_trajectories_ids[i, :len(trajectory), :] = torch.Tensor(trajectory)

print("Padding Done")

X = padded_trajectories
Y = padded_trajectories_ids
print("X and Y declared")


# padded_posterior_data = [pad_sequence(seq, max_trajectory_length, [0.0] * num_nodes) for seq in posterior_data]
# padded_one_hot_trajectories = [pad_sequence(seq, max_trajectory_length, [0] * num_nodes) for seq in one_hot_trajectories]



# Convert data to PyTorch tensors
# X = torch.tensor(padded_posterior_data, dtype=torch.float32)
# Y = torch.tensor(padded_one_hot_trajectories, dtype=torch.float32)

# Define a more complex Bidirectional LSTM model with 5 layers
class BiLSTM:
    def __init__(self, input_size, hidden_size, output_size, num_layers, dropout):
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size
        self.num_layers = num_layers
        self.dropout = dropout
        self.bilstm = self.build_bilstm()
        self.linear = torch.nn.Parameter(torch.Tensor(hidden_size * 2, output_size))  # Multiply by 2 for bidirectional
        self.init_weights()

    def build_bilstm(self):
        return torch.nn.LSTM(
            input_size=self.input_size,
            hidden_size=self.hidden_size,
            num_layers=self.num_layers,
            batch_first=True,
            dropout=self.dropout,
            bidirectional=True  # Use bidirectional
        )

    def init_weights(self):
        torch.nn.init.xavier_normal_(self.linear)

    def forward(self, input_data):
        output, _ = self.bilstm(input_data)
        output = torch.matmul(output, self.linear)
        return output

print("BiLSTM class declared")

# Define hyperparameters
input_size = num_nodes
hidden_size = 128  # Reduced hidden size (because it's doubled in bidirectional)
output_size = num_nodes
num_layers = 5  # Number of LSTM layers
dropout = 0.2  # Dropout rate
num_epochs = 50
learning_rate = 0.001



print("Hyper-parameters declared")

# Initialize the model and optimizer
model = BiLSTM(input_size, hidden_size, output_size, num_layers, dropout)
optimizer = optim.Adam(model.bilstm.parameters(), lr=learning_rate)

print("Initiating the training loop")

# Training loop
for epoch in range(num_epochs):
    optimizer.zero_grad()
    for i in range(num_trajectories):
        input_data = X[i].unsqueeze(0)
        target = Y[i].unsqueeze(0)

        output = model.forward(input_data)
        output = output.view(-1, num_nodes)
        target = Y[i].view(-1, num_nodes)
        loss = torch.nn.CrossEntropyLoss()(output, target)

        loss.backward()
        optimizer.step()
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
    # print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

# Save the trained model (optional)
torch.save(model.state_dict(), './Models_BiLSTM/bilstm_model1.pth')

print("Intial Model prepared")



k = 2
while k <= 21 :
    # Reading Posterior Vectors
    filename = 'pv'+str(k)+'.mat'
    path = '../Posteriors Dataset/'+str(filename)
    mat_file = mat73.loadmat(path)
    array_name = 'px'+str(k)
    cell_array = mat_file[array_name]
    pv_data = []
    post_vectors = []
    for m in range(10000):
        element = cell_array[m]
        pv = []
        for j in range(len(element[0])):
            vec = element[0][j]
            pv.append((vec).tolist())
        post_vectors.append(pv)
    print('File '+filename+' is readed')

    posterior_data = post_vectors
    max_trajectory_length = max(max(new_trajectories, key=len))


    trajectory_lengths = [len(trajectory) for trajectory in post_vectors]
    max_length = max(trajectory_lengths)
    padded_trajectories = torch.zeros((num_trajectories, max_length, num_nodes))
    for i, trajectory in enumerate(post_vectors):
        padded_trajectories[i, :len(trajectory), :] = torch.Tensor(trajectory)

    X = padded_trajectories
    Y = padded_trajectories_ids

    model = BiLSTM(input_size, hidden_size, output_size, num_layers, dropout)

    # Load the saved state_dict into the model.
    filename = 'bilstm_model'+str(k-1)+'.pth'
    path = './Models_BiLSTM/'+str(filename)
    saved_state_dict_path = path  # Provide the path to your saved state_dict.
    saved_state_dict = torch.load(saved_state_dict_path)
    model.load_state_dict(saved_state_dict)
    optimizer = optim.Adam(model.bilstm.parameters(), lr=learning_rate)

    print("Initiating the Loop")

    for epoch in range(num_epochs):
        optimizer.zero_grad()
        for i in range(num_trajectories):
            input_data = X[i].unsqueeze(0)
            target = Y[i].unsqueeze(0)
            output = model.forward(input_data)
            output = output.view(-1, num_nodes)
            target = Y[i].view(-1, num_nodes)
            loss = torch.nn.CrossEntropyLoss()(output, target)
            loss.backward()
            optimizer.step()
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
    
    print("{k} Batch training done")
    
    # Save only the model's state dictionary
    filename_model_state_dict = 'model_state_dict'+str(k)+'.pth'
    model_state_dict = './Models_BiLSTM/'+str(filename_model_state_dict)

    torch.save(model.state_dict(), model_state_dict)
    print(str(k)+"th batch trained")
    k = k+1
