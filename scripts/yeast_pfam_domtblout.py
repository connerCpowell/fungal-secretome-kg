import pandas as pd

rows = []

with open("domain_annotation/pfam.domtblout") as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split()

        target_name = fields[0]     # Pfam model
        target_acc  = fields[1]     # PFxxxxx.xx
        query_name  = fields[3]     # Protein ID

        evalue = float(fields[6])
        score = float(fields[7])

        hmm_from = int(fields[15])
        hmm_to   = int(fields[16])
        ali_from = int(fields[17])
        ali_to   = int(fields[18])

        rows.append([
            query_name,
            target_acc,
            evalue,
            score,
            hmm_from,
            hmm_to,
            ali_from,
            ali_to
        ])

df = pd.DataFrame(rows, columns=[
    "protein_id",
    "domain_id",
    "evalue",
    "score",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to"
])

df.to_csv(
    "data/processed/yeast_pfam_domains2.tsv",
    sep="\t",
    index=False
)

