
class LegacyGIN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, num_layers, out_channels, dropout):
        super().__init__()
        self.convs = ModuleList()
        # The Checkpoint has `GIN_p.convs.0.nn.0.weight`. So `convs` list exists.
        
        for i in range(num_layers):
            start_dim = in_channels if i == 0 else hidden_channels
            # Match strict sequential structure: 0:Lin, 1:BN, 2:ReLU, 3:Lin, 4:ReLU
            nn = Sequential(
                Linear(start_dim, hidden_channels),
                torch.nn.BatchNorm1d(hidden_channels),
                ReLU(),
                Linear(hidden_channels, hidden_channels),
                ReLU(),
            )
            self.convs.append(GINConv(nn, train_eps=True))
        self.dropout = dropout

    def forward(self, x, edge_index, *args, **kwargs):
        for i, conv in enumerate(self.convs):
            x = conv(x, edge_index)
            if i < len(self.convs) - 1:
                x = F.relu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
        return x

class GINpKa(LegacyGIN):
