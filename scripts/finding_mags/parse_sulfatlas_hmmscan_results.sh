#!/bin/bash

#######################################
# to parse sulfatlas results. They do not provide the hmm DB so i could not use hmmscan. I manually run each MAG on their webserver, downaloaded each result and parsed it myself
# Input is .txt file from their webpage
#######################################



set -euo pipefail

indir="${1:-./faa_files_putative_mucins/peptidase_sulfatase_output_sep_16_25/sulfatlas_webpage_hmmscan/sulfatlas_hmm_results}"
outdir="${2:-$indir/sulfatlas_parsed}"
EVALUE_THRESH="${EVALUE_THRESH:-1e-3}"

mkdir -p "$outdir"
shopt -s nullglob

for f in "$indir"/*.txt; do
  awk -v FS='[[:space:]]+' -v OFS='\t' -v thr="$EVALUE_THRESH" '
    # domtblout fields used:
    # 1=target name, 3=tlen, 4=query name, 7=E-value, 16=hmm_from, 17=hmm_to
    { sub(/\r$/, "") }            # strip CR if present
    /^#/ { next }                 # skip comments (if any)

    {
      t=$1; tlen=$3+0; q=$4; e=$7+0; hf=$16+0; ht=$17+0;
      if (e <= thr && tlen > 0) {
        cov = (ht - hf) / tlen;   # fraction (no +1)
        # Keep the BEST (LOWEST) E-value per query; tie-break by higher coverage
        if (!(q in bestE) || e < bestE[q] || (e == bestE[q] && cov > bestC[q])) {
          bestE[q]=e; bestT[q]=t; bestC[q]=cov;
          if (!(q in seen)) { order[++m]=q; seen[q]=1 }
        }
      }
    }

    END {
      for (i=1; i<=m; i++) {
        q = order[i];
        printf "%s\t%s\t%.6f\n", bestT[q], q, bestC[q];
      }
    }
  ' "$f" > "$outdir/$(basename "$f")"
done

echo "Done. Wrote outputs to: $outdir"
