from collections import defaultdict

pfam2go = defaultdict(list)

with open("data/raw/pfam2go") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("!"):
            continue

        # Split left/right of >
        try:
            left, right = line.split(">", 1)
        except ValueError:
            continue

        # Extract Pfam ID
        # Example left: Pfam:PF00001 7tm_1
        pfam_id = left.split()[0].replace("Pfam:", "")

        # Extract GO name and GO ID
        # Example right: GO:G protein-coupled receptor activity ; GO:0004930
        if ";" not in right:
            continue

        go_name_part, go_id_part = right.split(";", 1)

        go_name = go_name_part.replace("GO:", "").strip()
        go_id = go_id_part.strip()

        if go_id.startswith("GO:"):
            pfam2go[pfam_id].append((go_id, go_name))

# Write output
with open("data/processed/pfam_to_go.tsv", "w") as out:
    out.write("pfam_id\tgo_id\tgo_name\n")
    for pfam, gos in pfam2go.items():
        for go_id, go_name in gos:
            out.write(f"{pfam}\t{go_id}\t{go_name}\n")

