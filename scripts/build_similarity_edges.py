import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize

EMBED_FILE = "/home/funda/Projects/kg/fungal-secretome-kg/embeddings/yeast_secreted_esm2.npz"
OUTFILE = "/home/funda/Projects/kg/fungal-secretome-kg/data/processed/protein_similarity_edges.tsv"

TOP_K = 5
MIN_SIM = 0.80

# Load embeddings
data = np.load(EMBED_FILE)
protein_ids = list(data.keys())

X = np.vstack([data[p] for p in protein_ids])

# Normalize vectors (important)
X = normalize(X)

# Cosine similarity matrix
S = cosine_similarity(X)

edges = []

for i, pid in enumerate(protein_ids):
    sims = S[i]

    # Sort neighbors (excluding self)
    idxs = np.argsort(sims)[::-1]

    count = 0
    for j in idxs:
        if i == j:
            continue
        if sims[j] < MIN_SIM:
            break

        edges.append({
            "protein_a": pid,
            "protein_b": protein_ids[j],
            "cosine_similarity": float(sims[j])
        })

        count += 1
        if count >= TOP_K:
            break

df = pd.DataFrame(edges)
df.to_csv(OUTFILE, sep="\t", index=False)

print(df.head())
print(f"✅ {len(df)} similarity edges written")

