import pandas as pd

go = pd.read_csv("data/processed/yeast_protein_go.tsv", sep="\t")

keywords = [
    "hydrolase",
    "carbohydrate",
    "cell wall",
    "extracellular",
    "polysaccharide"
]

mask = go["go_name"].str.lower().apply(
    lambda x: any(k in x for k in keywords)
)

filtered = go[mask]

filtered.to_csv(
    "data/processed/yeast_protein_go_secretome.tsv",
    sep="\t",
    index=False
)

