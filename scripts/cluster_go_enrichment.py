"""
Cluster ESM embeddings of secreted yeast proteins
and perform GO enrichment per cluster.

Pipeline:
1. Load ESM embeddings
2. Cluster proteins by sequence similarity
3. Load GO annotations
4. Join clusters + GO terms
5. Perform GO enrichment (Fisher exact test)

Inputs:
  - data/processed/yeast_secreted_esm2.npz
  - data/processed/yeast_protein_go_secretome.tsv

Outputs:
  - data/processed/yeast_protein_clusters.tsv
  - data/processed/yeast_cluster_go_enrichment.tsv
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.stats import fisher_exact

# -----------------------------
# PARAMETERS
# -----------------------------

K = 5                 # number of clusters
MIN_GO_COUNT = 2      # minimum cluster hits for GO term


# -----------------------------
# STEP C1 — Load embeddings
# -----------------------------
def load_esm_embeddings(path):
    """
    Loads ESM embeddings from .npz file.

    Returns:
      protein_ids : list[str]
      X           : np.ndarray (N x 1280)
    """
    emb = np.load(path)
    protein_ids = list(emb.keys())
    X = np.vstack([emb[p] for p in protein_ids])
    return protein_ids, X


# -----------------------------
# STEP C1 — Cluster embeddings
# -----------------------------
def cluster_embeddings(protein_ids, X, k):
    """
    Clusters proteins using KMeans.

    Returns:
      cluster_df : DataFrame with columns
                   [protein_id, cluster_id]
    """
    model = KMeans(n_clusters=k, random_state=42)
    labels = model.fit_predict(X)

    cluster_df = pd.DataFrame({
        "protein_id": protein_ids,
        "cluster_id": [f"C{c}" for c in labels]
    })

    return cluster_df


# -----------------------------
# STEP C2 — Load GO annotations
# -----------------------------
def load_go_annotations(path):
    """
    Loads protein → Pfam → GO mappings.

    Returns:
      go_df : DataFrame
    """
    return pd.read_csv(path, sep="\t")


# -----------------------------
# STEP C2 — Merge GO + clusters
# -----------------------------
def merge_clusters_and_go(cluster_df, go_df):
    """
    Joins cluster labels to GO annotations.

    Returns:
      merged_df : DataFrame
      background_proteins : set[str]
    """
    merged_df = go_df.merge(cluster_df, on="protein_id")
    background_proteins = set(go_df["protein_id"])
    return merged_df, background_proteins


# -----------------------------
# STEP C3/C4 — GO enrichment
# -----------------------------
def compute_go_enrichment(merged_df, background_proteins):
    """
    Performs GO enrichment per cluster using Fisher's exact test.

    Returns:
      enrichment_df : DataFrame
    """
    results = []

    for cluster_id in merged_df["cluster_id"].unique():
        cluster_proteins = set(
            merged_df.loc[
                merged_df["cluster_id"] == cluster_id,
                "protein_id"
            ]
        )
        cluster_size = len(cluster_proteins)

        for (go_id, go_name), sub in merged_df.groupby(
            ["go_id", "go_name"]
        ):
            go_proteins = set(sub["protein_id"])

            # contingency table
            a = len(cluster_proteins & go_proteins)
            if a < MIN_GO_COUNT:
                continue

            b = cluster_size - a
            c = len(go_proteins - cluster_proteins)
            d = len(background_proteins) - (a + b + c)

            _, p_value = fisher_exact(
                [[a, b], [c, d]],
                alternative="greater"
            )

            results.append({
                "cluster_id": cluster_id,
                "go_id": go_id,
                "go_name": go_name,
                "cluster_size": cluster_size,
                "go_count_in_cluster": a,
                "p_value": p_value
            })

    enrichment_df = (
        pd.DataFrame(results)
        .sort_values(["cluster_id", "p_value"])
    )

    return enrichment_df


# -----------------------------
# MAIN PIPELINE
# -----------------------------
def main():
    print("Loading ESM embeddings...")
    protein_ids, X = load_esm_embeddings(
        "embeddings/yeast_secreted_esm2.npz"
    )

    print(f"Clustering {len(protein_ids)} proteins into {K} clusters...")
    cluster_df = cluster_embeddings(protein_ids, X, K)

    cluster_df.to_csv(
        "data/processed/yeast_protein_clusters.tsv",
        sep="\t",
        index=False
    )

    print("Loading GO annotations...")
    go_df = load_go_annotations(
        "data/processed/yeast_protein_go_secretome.tsv"
    )

    merged_df, background = merge_clusters_and_go(
        cluster_df, go_df
    )

    print("Computing GO enrichment...")
    enrichment_df = compute_go_enrichment(
        merged_df, background
    )

    enrichment_df.to_csv(
        "data/processed/yeast_cluster_go_enrichment.tsv",
        sep="\t",
        index=False
    )

    print("Pipeline complete.")


if __name__ == "__main__":
    main()

