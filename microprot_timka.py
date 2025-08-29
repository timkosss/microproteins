import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# List amino acids, map to indices (random, may be improved)
AMINO_ACIDS = 'IVLMAGPFWYCTSQNKRHED' ##sorted by functional groups and polarity
# (Nonpolar aliphatic (hydrophobicity high→low):IVLMAGP, Aromatic (hydrophobicity high→low):FWY, Polar uncharged (hydrophobicity high→low):
# CTSQN Positively charged (hydrophobicity high→low):KRH Negatively charged (hydrophobicity high→low):ED

AA_TO_IDX = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}

## some params
SEQ_LEN = int(sys.argv[5])
VOCAB_SIZE = len(AMINO_ACIDS)  # 20
batch_size = 64
epochs = int(sys.argv[1])
lr = float(sys.argv[2])



# ------------ functions and classes ------------

def filter_valid_proteins(sequences):
    """
    Filters a list of protein sequences, keeping only those composed of valid amino acids.

    Parameters
    ----------
    sequences : list of str
        List of protein sequences (strings of amino acid letters).

    Returns
    -------
    list of str
        List containing only sequences that have valid amino acid letters.
        Valid amino acids:IVLMAGPFWYCTSQNKRHED.
    """
    valid_aas = set(AMINO_ACIDS)
    filtered = []
    for seq in sequences:
        seq_set = set(seq.upper())
        if seq_set.issubset(valid_aas):
            filtered.append(seq)
    return filtered

def read_fasta_first_k(fasta_path,k=10):
    """
    Reads sequences from a FASTA file and returns only the first `k` residues of each sequence.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.
    k : int, optional (default=10)
        Number of amino acids to keep from the start of each sequence.

    Returns
    -------
    list of str
        List of truncated sequences (length <= k).
    """
    sequences = []
    with open(fasta_path, 'r') as f:
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences.append(seq[:k+1])
                    seq = ''
            else:
                seq += line
        # Add last
        if seq:
            sequences.append(seq[:k+1])
    return sequences

class ProteinDataset(Dataset):
    """
    PyTorch Dataset for protein classification tasks.

    Each item returns:
      - Encoded amino acid indices (tensor)
      - Label (1 <-->, 0 <--> negative)

    Parameters
    ----------
    positive_seqs : list of str
        List of positive protein sequences (of len SEQ_LEN)
    negative_seqs : list of str
        List of negative protein sequences (of len SEQ_LEN).

    Attributes
    ----------
    sequences : list of str
        All protein sequences (positive + negative)
    labels : list of int
        Corresponding binary labels for each sequence
    """
    def __init__(self, positive_seqs, negative_seqs):
        self.sequences = positive_seqs + negative_seqs
        self.labels = [1] * len(positive_seqs) + [0] * len(negative_seqs)
        assert all(len(seq) == SEQ_LEN for seq in self.sequences), "All sequences must be of length "+str(SEQ_LEN)

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq = self.sequences[idx]
        seq_idx = torch.tensor([AA_TO_IDX[aa] for aa in seq], dtype=torch.long)
        label = torch.tensor(self.labels[idx], dtype=torch.float)
        return seq_idx, label

dropout = float(sys.argv[6])

class ProteinModel(nn.Module):
    """
       LSTM-based neural net for protein sequence classification with dropout

       Architecture:
         - Embedding layer: Converts AA to dense vectors...  with dropout
         - LSTM layer: models sequence of AA and long-term dependencies
         - Classifier: Fully connected layers with ReLU activation
         - Sigmoid: Outputs probability of positive class (can replace with softmax)

       Parameters
       ----------
       vocab_size : int
           Number of unique amino acids in the vocabulary.
       embed_dim : int, optional (default=64)
           Dimension of amino acid embeddings
       hidden_dim : int, optional (default=128) --> higher gets very slow and worse
           Hidden dimension of the LSTM
    """

    def __init__(self, vocab_size, embed_dim=64, hidden_dim=128, dropout_prob=dropout):
        super().__init__()
        # Embedding layer
        self.embedding = nn.Embedding(vocab_size, embed_dim)
        self.embedding_dropout = nn.Dropout(p=dropout_prob)  # dropout after embedding, better regularization

        # LSTM
        self.lstm = nn.LSTM(
            embed_dim,
            hidden_dim,
            batch_first=True,
            num_layers=1,   # do not stack LSTMs!!
            dropout=0.0     # dropout inside LSTM if num_layers > 1
        )

        # Classifier with dropout
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim, 64),
            nn.ReLU(),
            nn.Dropout(p=dropout_prob),   # dropout classifier
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        embedded = self.embedding(x)               # [batch, seq_len, embed_dim]
        embedded = self.embedding_dropout(embedded)
        _, (hn, _) = self.lstm(embedded)
        hn = hn.squeeze(0)                         #  hidden state
        out = self.classifier(hn)                  #  classifier
        return out.squeeze(1)



def train(model, dataloader, criterion, optimizer, device):
    """
      Trains the protein model per epoch.

      Parameters
      ----------
      model : ProteinModel
          Model to train
      dataloader : DataLoader
          PyTorch DataLoader for the training data.
      criterion : loss function
          Loss function (e.g., BCELoss or BCEWithLogitsLoss).
      optimizer : torch.optim.Optimizer
          Optimizer (e.g., Adam).
      device : torch.device
          Device to run training on ('cpu' or 'cuda').

      Returns
      -------
      float
          Average training loss over the dataset.
      """
    model.train()
    total_loss = 0
    for seqs, labels in dataloader:
        seqs = seqs.to(device)
        labels = labels.to(device)
        optimizer.zero_grad()
        outputs = model(seqs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * seqs.size(0)
    return total_loss / len(dataloader.dataset)



class ProteinDataset_testing(Dataset):
    """
    PyTorch Dataset for inference on protein sequences (without labels). Good to keep train,test seperate.

    Parameters
    ----------
    sequences : list of str
        Protein sequences (of len SEQ_LEN).
    """
    def __init__(self, sequences):
        self.sequences = sequences
        assert all(len(seq) == SEQ_LEN for seq in self.sequences), "All sequences must be if length "+SEQ_LEN
    def __len__(self):
        return len(self.sequences)
    def __getitem__(self, idx):
        seq = self.sequences[idx]
        seq_idx = torch.tensor([AA_TO_IDX[aa] for aa in seq], dtype=torch.long)
        return seq_idx

def predict(model, sequences, device):
    """
    Applies model to a list of protein sequences for prediction of short/long protein

    Parameters
    ----------
    model : ProteinModel
        Trained model used for prediction
    sequences : list of str
        Protein sequences to predict using model
    device : torch.device
        Device for running

    Returns
    -------
    probs : list of float
        Predicted probabilities of positive class.
    preds : list of int
        Binary predictions (0/1) using threshold of 0.5.
    """
    dataset = ProteinDataset_testing(sequences)
    loader = DataLoader(dataset, batch_size, shuffle=False)
    model.eval()
    preds = []
    probs = []
    with torch.no_grad():
        for seqs in loader:
            seqs = seqs.to(device)
            outputs = model(seqs)  # outputs probabilities (sigmoid)
            probs.extend(outputs.cpu().numpy())
            preds.extend((outputs >= 0.5).cpu().numpy().astype(int))
    return probs, preds


# ------------ training model ------------

pos_s = pd.read_csv('/wistar/auslander/microproteins/TableS2c.csv')
pos_s['len'] = [len(i) for i in pos_s['Sequence']]
pos_s=pos_s[pos_s['len']>=SEQ_LEN]

df_sampled = pos_s.groupby("Cluster_ID", group_keys=False).apply(lambda x: x.sample(1, random_state=42))


ps=list(df_sampled['Sequence']) #read the sequences of short proteins from: https://pubmed.ncbi.nlm.nih.gov/39978337/

positive_seqs=filter_valid_proteins([i[1:1+SEQ_LEN] for i in ps if len(i)>=SEQ_LEN+1 ]) #take the first K (10) AA of short proteins with at least K AA

negative_database = str(sys.argv[3])

ns=read_fasta_first_k(f'/wistar/auslander/Timothy/microproteins_timka/{negative_database}', k=SEQ_LEN)[:len(positive_seqs)] #read all sequences of long proteins (>=KX10) from the same bacteria
negative_seqs = filter_valid_proteins([i[1:1+SEQ_LEN] for i in ns if len(i)>=SEQ_LEN+1])

num_training = int(sys.argv[4])

train_pos=positive_seqs[:num_training] #training, for now the first 1K (should cluster and split) use Cluster_ID for microprots, and need to cluster the larger db
train_neg=negative_seqs[:num_training] #training, for now the first 1K (should cluster and split)

print(f"Sequence Length: {SEQ_LEN}\nEpochs: {epochs}\nLearning Rate: {lr}\nAmount of Training Data: {num_training}\nNumber of Sequences per Dataset: {len(positive_seqs)}\nNegative Database Used: {negative_database}")


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

train_dataset = ProteinDataset(train_pos, train_neg)
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

model = ProteinModel(VOCAB_SIZE).to(device)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=lr)

for epoch in range(1, epochs + 1):
    train_loss = train(model, train_loader, criterion, optimizer, device)
    print(f"Epoch {epoch}: Train loss={train_loss:.4f}")



# ------------testing model ------------

positive_test_seqs = positive_seqs[num_training:]
negative_test_seqs = negative_seqs[num_training:]

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Predict on positive and negative test sequences separately
pos_probs, pos_preds = predict(model, positive_test_seqs, device)
neg_probs, neg_preds = predict(model, negative_test_seqs, device)

# calc AUROC, ToDo: consider other metrics for evaluation
labels = [1 for i in range(len(pos_probs))]+[0 for i in range(len(neg_probs))]
scores = pos_probs+neg_probs
auroc = roc_auc_score(labels, scores)

print(f"AUROC: {auroc:.4f}")

with open('/wistar/auslander/Timothy/microproteins_timka/hyperparameter_results_negative_num_training.csv', 'a') as f:
    f.write(f"{auroc:.4f},{train_loss:.4f},{SEQ_LEN},{epoch},{lr},{negative_database},{num_training},{dropout}\n")