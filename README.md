# Graph-based Functional Inference of Secreted Fungal Enzymes Using Protein Language Model Embeddings


![Fungi](https://upload.wikimedia.org/wikipedia/commons/7/73/Saccharomyces_cerevisiae_cells.jpg)

**A knowledge graph of secreted proteins in *Saccharomyces cerevisiae* integrating sequence embeddings, Pfam domains, and functional annotations.**

---

## 🚀 Project Overview

This project demonstrates the construction of a **knowledge graph (KG)** for the **fungal secretome**, with the goal of linking **protein sequences → embeddings → clusters → functional annotations (GO terms)**. It showcases capabilities in:

- **Bioinformatics preprocessing**: secretome prediction, domain annotation  
- **Protein embeddings**: ESM embeddings for sequence similarity  
- **Functional enrichment**: GO term enrichment per protein cluster  
- **Knowledge graph construction**: ingesting into Neo4j for visualization and exploration  

The final KG allows **exploration of functional patterns, cluster relationships, and protein similarity networks**.

---

## 🧬 Dataset

- Protein sequences for *S. cerevisiae* (secretome subset)  
- Pfam domain annotations  
- GO term mappings (via `pfam2go`)  
- SignalP-5.0 predictions for secretion  

**Processed files** (stored in `/data/processed`):

- `yeast_protein_go.tsv` → protein → GO mapping  
- `yeast_protein_clusters.tsv` → protein → cluster mapping  
- `yeast_cluster_go_enrichment.tsv` → cluster → enriched GO terms  
- `yeast_secreted_esm2.npz` → ESM protein embeddings  

---

## ⚙️ Workflow

1. **Secreted protein prediction**  
   - Run SignalP on raw sequences  
   - Filter proteins predicted as secreted

2. **Protein embeddings**  
   - Compute ESM embeddings  
   - Save as `.npz` for clustering

3. **Pfam domain annotation**  
   - Map proteins to Pfam domains  
   - Map Pfam domains to GO terms (`pfam2go`)  
   - Create `yeast_protein_go_secretome.tsv`

4. **Clustering**  
   - Cluster ESM embeddings (KMeans, k=5)  
   - Save protein → cluster mapping

5. **GO enrichment analysis**  
   - Compute per-cluster GO enrichment  
   - Save results (`yeast_cluster_go_enrichment.tsv`)

6. **Knowledge graph construction in Neo4j**  
   - Nodes: `Protein`, `Cluster`, `GO`  
   - Relationships:  
     - `:IN_CLUSTER` → Protein → Cluster  
     - `:HAS_GO` → Protein → GO  
     - `:ENRICHED_FOR` → Cluster → GO  
   - Optional: future cross-species expansion

---

## 🗂 Project Structure

ungal-secretome-kg/
│
├─ data/
│ ├─ raw/ # raw sequences and pfam2go
│ └─ processed/ # TSVs for ingestion
│
├─ scripts/
│ ├─ compute_esm_embeddings.py
│ ├─ annotate_proteins_with_go.py
│ ├─ cluster_go_enrichment.py
│ └─ build_similarity_edges.py
│
├─ tools/ # SignalP binaries, ESM models
│
├─ Dockerfile / docker-compose # Neo4j setup
│
└─ README.md


---

## 📊 Graph Schema

**Nodes**:

| Label   | Properties                     |
|---------|--------------------------------|
| Protein | protein_id                     |
| Cluster | cluster_id                     |
| GO      | go_id, go_name                 |

**Relationships**:

| Type           | From     | To      | Properties                        |
|----------------|---------|---------|-----------------------------------|
| IN_CLUSTER      | Protein | Cluster | —                                 |
| HAS_GO          | Protein | GO      | —                                 |
| ENRICHED_FOR    | Cluster | GO      | p_value, cluster_size, go_count   |

---

## 🔍 Neo4j Visualization

1. Start Neo4j (Docker recommended).  
2. Import processed TSVs using `LOAD CSV WITH HEADERS`.  
3. Explore graph using queries like:

```cypher
MATCH (p:Protein)-[:IN_CLUSTER]->(c:Cluster)
RETURN p, c
LIMIT 25;


Future Work

Expand to other fungal species

Add protein similarity edges based on embeddings

Cross-species cluster comparison and GO enrichment analysis

Interactive web-based visualization of KG

References

SignalP 5.0: https://services.healthtech.dtu.dk/service.php?SignalP-5.0

Pfam Database: https://pfam.xfam.org/

ESM Protein Language Models: https://github.com/facebookresearch/esm

Neo4j Graph Database: https://neo4j.com/
