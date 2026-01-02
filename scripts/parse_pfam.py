import pandas as pd

infile = "/home/funda/Projects/kg/fungal-secretome-kg/domain_annotation/pfam.domtblout"
outfile = "/home/funda/Projects/kg/fungal-secretome-kg/data/processed/yeast_pfam_domains.tsv"

rows = []

with open(infile) as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split()
        rows.append({
            "protein_id": parts[0],
            "domain_id": parts[1],
            "domain_name": parts[2],
            "evalue": float(parts[6]),
            "score": float(parts[7]),
            "hmm_from": int(parts[15]),
            "hmm_to": int(parts[16]),
            "ali_from": int(parts[17]),
            "ali_to": int(parts[18])
        })

df = pd.DataFrame(rows)

# Filter high-confidence domains
df = df[
    (df["evalue"] < 1e-5) &
    ((df["ali_to"] - df["ali_from"]) > 30)
]

df.to_csv(outfile, sep="\t", index=False)

print(df.head())
print(f"✅ {df['protein_id'].nunique()} proteins with Pfam domains")

