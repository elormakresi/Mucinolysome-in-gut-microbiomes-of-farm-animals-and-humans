import pandas as pd
from pathlib import Path

#by jerry after Siva explained and used to make the final data file  that can be used to generate heatmaps with cper_combined_images_all_6.py

conf_tsv = "other_ppi_scores_mutput.tsv"   # source_file, iptm, ptm, etc. the one that get_other_metrics_from_json.py will generate
contacts_tsv = "ppi_scores.tsv"  # final_source, num_contacts, avg_if_plddt, pdockq

# 1) Load the confidences table
df1 = pd.read_csv(conf_tsv, sep="\t")

# Derive final_source from source_file path:
# e.g., /array3/.../cpegh20_coh-cpegh31_doc_summary_confidences.json -> cpegh20_coh-cpegh31_doc
def to_final_source(p):
    name = Path(str(p)).name
    return name.replace("_summary_confidences.json", "")

df1["final_source"] = df1["source_file"].map(to_final_source)

# Keep the needed columns (and coerce numeric just in case)
df1 = df1[["final_source", "fraction_disordered", "has_clash", "iptm", "ptm", "ranking_score"]].copy()
for c in ["fraction_disordered", "has_clash", "iptm", "ptm", "ranking_score"]:
    df1[c] = pd.to_numeric(df1[c], errors="coerce")

# Optional: force has_clash to 0/1 integers if it's 0.0/1.0
df1["has_clash"] = df1["has_clash"].fillna(0).round().astype(int)

# 2) Load the contacts table (no header)
df2 = pd.read_csv(
    contacts_tsv, sep="\t", header=None,
    names=["final_source", "num_contacts", "avg_if_plddt", "pdockq"]
)
for c in ["num_contacts", "avg_if_plddt", "pdockq"]:
    df2[c] = pd.to_numeric(df2[c], errors="coerce")

# 3) Merge on final_source
out = pd.merge(df1, df2, on="final_source", how="inner")

# 4) Reorder columns to your desired output
out = out[
    ["final_source", "fraction_disordered", "has_clash",
     "iptm", "ptm", "ranking_score",
     "num_contacts", "avg_if_plddt", "pdockq"]
]

# 5) Save as TSV
out.to_csv("final_metrics.tsv", sep="\t", index=False)

# Preview
print(out.head())
