import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt

# Assuming F, G, and D are available as PyTorch Tensors

matrix_dims = (200, 200)
rotation_matrix_to_find = torch.rand(matrix_dims)

# Define the Model
class SimplifiedFluidModel(nn.Module):
    def __init__(self):
        super(SimplifiedFluidModel, self).__init__()
        # A is a learnable matrix. Its dimension depends on your specific problem.
        # Assuming A is a 3x3 matrix for this example
        self.A = nn.Parameter(torch.randn(matrix_dims))

    def forward(self, G, F, D):
        # G.F (matrix multiplication)
        # GF = torch.matmul(G, F)

        # A.GF (matrix multiplication)
        AGF = torch.matmul(G, self.A)

        # V(p) = AGF * D
        Vp = AGF # * D
        return Vp

# Initialize the model
model = SimplifiedFluidModel()

# Define a Loss function and optimizer
criterion = nn.MSELoss()
optimizer = optim.SGD(model.parameters(), lr=0.01)

num_epochs = 10000

def createData():
    Gs = torch.rand(10, matrix_dims[0])
    Fs = torch.tensor([1] * matrix_dims[0] * len(Gs)).resize(10, matrix_dims[0])
    Ds = torch.tensor([0 for _ in range(len(Gs))])
    # Define the 90 degree rotation matrix
    # rotation_matrix = torch.tensor([[1, -1, 0], [1, 0, 1], [8, 0, 1]], dtype=torch.float32)
    # Apply the rotation to all vectors
    targets = torch.matmul(Gs, rotation_matrix_to_find) #rotation_matrix)
    return Gs, Fs, Ds, targets

G, F, D, target_velocities = createData()

losses = []
# Training loop
for epoch in range(num_epochs):
    # Zero the parameter gradients
    optimizer.zero_grad()

    # Forward pass
    outputs = model(G, F, D)

    # Compute the loss
    loss = criterion(outputs, target_velocities) # target_velocities should be the ground truth

    # Backward pass and optimize
    loss.backward()
    optimizer.step()

    losses.append(loss.item())

    # Print statistics
    if (epoch+1) % 100 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
print(model.A - rotation_matrix_to_find)

plt.plot(np.array(losses))
plt.show()