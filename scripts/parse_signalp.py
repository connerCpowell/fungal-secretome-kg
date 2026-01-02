import pandas as pd

infile = "/home/funda/Projects/kg/fungal-secretome-kg/data/processed/y_secreted_summary.signalp5"
outfile = "/home/funda/Projects/kg/fungal-secretome-kg/data/processed/yeast_signalp.tsv"

rows = []

with open(infile) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        parts = line.strip().split()

        protein_id = parts[0]
        prediction = parts[1]
        sp_prob = float(parts[2])
        other_prob = float(parts[3])

        # Cleavage site probability (optional)
        cleavage_prob = None
        for token in parts:
            if token.startswith(" Pr:"):
                cleavage_prob = float(token.replace(" Pr:", ""))

        rows.append({
            "protein_id": protein_id,
            "prediction": prediction,
            "signalp_sp_prob": sp_prob,
            "signalp_other_prob": other_prob,
            "cleavage_prob": cleavage_prob,
            "is_secreted": prediction.startswith("SP")
        })

df = pd.DataFrame(rows)

# Confidence classes (for graph reasoning later)
df["secreted_confidence_class"] = pd.cut(
    df["signalp_sp_prob"],
    bins=[0, 0.5, 0.8, 1.0],
    labels=["low", "medium", "high"]
)

df.to_csv(outfile, sep="\t", index=False)

print(df.head())
print("\nConfidence class counts:")
print(df["secreted_confidence_class"].value_counts())
