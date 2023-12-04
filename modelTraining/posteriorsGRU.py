
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
import scipy.io
import mat73

# from google.colab import drive
# drive.mount('/content/drive')
# %cd /content/drive/My Drive/dataset_for_posteriors



mat_file = mat73.loadmat('post_vectors.mat')
cell_array = mat_file['post_vectors']

post_vectors = []
for i in range(5000):
  element = cell_array[i]
  pv = []
  for j in range(len(element[0])):
    vec = element[0][j]
    pv.append((vec).tolist())
  post_vectors.append(pv)



mat_file = scipy.io.loadmat('new_trajectories.mat')

cell_array = mat_file['new_traj']
trajectories = []
for i in range(5000):
  traj = cell_array[i][0][0].tolist()
  trajectories.append(traj)


mat_file = scipy.io.loadmat('new_nodes.mat')

cell_array = mat_file['new_node_ids']
nodes_list = []
for i in range(len(cell_array)):
  nodes = cell_array[i][0]
  nodes_list.append(nodes)



node_to_index = {node: index for index, node in enumerate(nodes_list)}
new_trajectories = []
for trajectory in trajectories:
    new_traj = [node_to_index[node] for node in trajectory]
    new_trajectories.append(new_traj)


num_nodes = 1739
feature_length = 1739
num_trajectories = 5000  # Test for first 50 trajectories

# new_trajectories = new_trajectories[:50]    # Test for first 50 trajectories


one_hot_trajectories = []

for trajectory in new_trajectories:
    one_hot_trajectory = []
    for node_id in trajectory:
        one_hot_vector = np.zeros(num_nodes, dtype=int)
        one_hot_vector[node_id] = 1
        one_hot_trajectory.append(one_hot_vector)
    one_hot_trajectories.append(one_hot_trajectory)



trajectory_lengths = [len(trajectory) for trajectory in post_vectors]
max_length = max(trajectory_lengths)
padded_trajectories = torch.zeros((num_trajectories, max_length, feature_length))
for i, trajectory in enumerate(post_vectors):
    padded_trajectories[i, :len(trajectory), :] = torch.Tensor(trajectory)

padded_trajectories_ids = torch.zeros((num_trajectories, max_length, feature_length))
for i, trajectory in enumerate(one_hot_trajectories):
    padded_trajectories_ids[i, :len(trajectory), :] = torch.Tensor(trajectory)


class Encoder(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(Encoder, self).__init__()
        self.gru = nn.GRU(input_size, hidden_size, batch_first=True)

    def forward(self, input_trajectory):
        _, hidden = self.gru(input_trajectory)
        return hidden

class Decoder(nn.Module):
    def __init__(self, output_size, hidden_size):
        super(Decoder, self).__init__()
        self.gru = nn.GRU(output_size, hidden_size, batch_first=True)
        self.out = nn.Linear(hidden_size, output_size)

    def forward(self, input_node, hidden, context):
        output, hidden = self.gru(input_node, hidden)
        predicted_node = self.out(output)
        return predicted_node, hidden

class Seq2Seq(nn.Module):
    def __init__(self, encoder, decoder):
        super(Seq2Seq, self).__init__()
        self.encoder = encoder
        self.decoder = decoder

    def forward(self, input_trajectory, target_trajectory):
        hidden = self.encoder(input_trajectory)

        decoder_hidden = hidden
        context = None  # We don't use context in this version

        outputs = []
        for input_node in input_trajectory.permute(1, 0, 2):
            predicted_node, decoder_hidden = self.decoder(input_node.view(1, 1, -1), decoder_hidden, context)
            outputs.append(predicted_node)

        return torch.cat(outputs, dim=0)


input_size = feature_length
hidden_size = 256
output_size = num_nodes

encoder = Encoder(input_size, hidden_size)
decoder = Decoder(output_size, hidden_size)

model = Seq2Seq(encoder, decoder)

criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)


print("Training Initiated")

num_epochs = 5

for epoch in range(num_epochs):
    for i in range(len(padded_trajectories)):
        target_trajectory = padded_trajectories_ids[i]  # For this simple example, the target is the same as input
        optimizer.zero_grad()
        output = model(padded_trajectories[i].unsqueeze(0), target_trajectory.unsqueeze(0))
        # output = output.squeeze(1)
        output = output.view(-1, num_nodes)
        target_trajectory = target_trajectory.view(-1, num_nodes)
        loss = criterion(output, target_trajectory)
        loss.backward()
        optimizer.step()
        print(f"Epoch {epoch + 1}, Loss: {loss.item()}")

# Save the entire model
torch.save(model, 'model.pth')
# Save only the model's state dictionary
torch.save(model.state_dict(), 'model_state_dict.pth')

