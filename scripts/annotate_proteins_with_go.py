import pandas as pd

pfam = pd.read_csv(
    "data/processed/yeast_pfam_domains2.tsv", 
    sep="\t"
)

pfam["pfam_id"] = pfam["domain_id"].str.split(".").str[0]

pfam2go = pd.read_csv(
    "data/processed/pfam_to_go.tsv", 
    sep="\t"
)

annot = pfam.merge(
    pfam2go,
    on="pfam_id",
    how="inner"
)

annot = annot[["protein_id", 
               "domain_id", 
               "go_id", 
               "go_name"]].drop_duplicates()

annot.to_csv(
    "data/processed/yeast_protein_go.tsv",
    sep="\t",
    index=False
)

