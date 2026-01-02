import os
import sys
import torch
import esm
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# ----------------------------
# Config
# ----------------------------
FASTA_FILE = "/home/funda/Projects/kg/fungal-secretome-kg/data/raw/y_secreted.fasta"
SIGNALP_TSV = "/home/funda/Projects/kg/fungal-secretome-kg/data/processed/yeast_signalp.tsv"
OUT_EMB = "/home/funda/Projects/kg/fungal-secretome-kg/embeddings/yeast_secreted_esm2.npz"
OUT_META = "/home/funda/Projects/kg/fungal-secretome-kg/embeddings/metadata.tsv"

MIN_SIGNALP_PROB = 0.8
BATCH_SIZE = 4            # Safe for CPU
MODEL_NAME = "esm2_t33_650M_UR50D"

# ----------------------------
# Load SignalP filtering
# ----------------------------
signalp = pd.read_csv(SIGNALP_TSV, sep="\t")

allowed_ids = set(
    signalp[
        (signalp["is_secreted"] == True) &
        (signalp["signalp_sp_prob"] >= MIN_SIGNALP_PROB)
    ]["protein_id"]
)

if not allowed_ids:
    print("❌ No proteins passed SignalP filtering.")
    sys.exit(1)

print(f"✅ {len(allowed_ids)} proteins pass SignalP filtering")

# ----------------------------
# Load sequences
# ----------------------------
records = [
    r for r in SeqIO.parse(FASTA_FILE, "fasta")
    if r.id in allowed_ids
]

print(f"✅ {len(records)} sequences loaded")

# ----------------------------
# Resume support
# ----------------------------
existing = {}
if os.path.exists(OUT_EMB):
    npz = np.load(OUT_EMB, allow_pickle=True)
    existing = dict(npz)
    print(f"🔁 Resuming: {len(existing)} embeddings already exist")

# ----------------------------
# Load model
# ----------------------------
print("🚀 Loading ESM model...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()

# ----------------------------
# Embedding loop
# ----------------------------
embeddings = existing.copy()
metadata_rows = []

for i in tqdm(range(0, len(records), BATCH_SIZE)):
    batch = records[i:i+BATCH_SIZE]
    batch = [r for r in batch if r.id not in embeddings]

    if not batch:
        continue

    labels, strs, tokens = batch_converter(
        [(r.id, str(r.seq)) for r in batch]
    )

    with torch.no_grad():
        results = model(tokens, repr_layers=[33])

    reps = results["representations"][33]

    for j, record in enumerate(batch):
        emb = reps[j, 1:len(record.seq)+1].mean(0).cpu().numpy()
        embeddings[record.id] = emb

        metadata_rows.append({
            "protein_id": record.id,
            "length": len(record.seq),
            "model": MODEL_NAME
        })

# ----------------------------
# Save outputs
# ----------------------------
np.savez_compressed(OUT_EMB, **embeddings)
pd.DataFrame(metadata_rows).to_csv(OUT_META, sep="\t", index=False)

print(f"✅ Saved {len(embeddings)} embeddings")
print("🎉 ESM embedding generation complete")

