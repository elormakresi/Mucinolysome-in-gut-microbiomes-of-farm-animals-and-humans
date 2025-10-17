#!/usr/bin/env bash
#SBATCH --time=1:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=2
#SBATCH --job-name=tpm2bw
#SBATCH --partition=yinlab,batch,guest
set -euo pipefail

#Written to make tpm bigwigs instead of the CPM or RPKM ones.

GTF="../kol109.gtf"                       # has transcript rows with gene_id
FASTA="../kol109.fna"                     # genome FASTA
TPM_CSV="./TPM_matrix_average.csv"        # first column = gene_id, columns include G_avg,M_avg, must be obtained from the DESEQ2 calc.
COLS="G_avg,M_avg"                        # which TPM columns to export

# --- outputs ---
OUTDIR="./tpm_bigwigs"
BED_GENE="${OUTDIR}/genes.transcript.bed"
CHROMSIZES="${OUTDIR}/chrom.sizes"

mkdir -p "$OUTDIR"

# tools
ml samtools/1.9 || true
ml ucsc-bedgraphtobigwig/455 || true
command -v bedSort >/dev/null 2>&1 || echo "[info] bedSort not found; relying on UNIX sort."

# 1) transcript BED keyed by gene_id (6 columns)
awk -F'\t' 'BEGIN{OFS="\t"}
$0!~/^#/ && $3=="transcript"{
  gid="";
  if (match($9,/gene_id "[^"]+"/)) gid=substr($9,RSTART+9,RLENGTH-10);
  if(gid!=""){
    strand=$7; if (strand!="+" && strand!="-") strand=".";
    print $1, $4-1, $5, gid, 0, strand
  }
}' "$GTF" | LC_ALL=C sort -k1,1 -k2,2n > "$BED_GENE"

# 2) chrom sizes
samtools faidx "$FASTA"
cut -f1,2 "${FASTA}.fai" > "$CHROMSIZES"

# map header -> column index (1-based for awk)
header=$(head -n1 "$TPM_CSV")
IFS=',' read -r -a hdr <<< "$header"
declare -A col2idx
for i in "${!hdr[@]}"; do
  key=$(echo "${hdr[$i]}" | sed 's/"//g' | xargs)
  col2idx["$key"]=$((i+1))
done
IFS=',' read -r -a wantCols <<< "$COLS"

for cname in "${wantCols[@]}"; do
  idx=${col2idx[$cname]:-}
  [[ -n "${idx:-}" ]] || { echo "[ERR] Column '$cname' not found in header: $header" >&2; exit 1; }

  # gene_id -> TPM map  (FORCE TAB OUTPUT)
  awk -F',' -v I="$idx" 'BEGIN{OFS="\t"} NR>1{
    g=$1; sub(/^"*/,"",g); sub(/"*$/,"",g);
    v=$I; sub(/^"*/,"",v); sub(/"*$/,"",v);
    if(v==""||v=="NA"||v=="NaN"||v=="nan") v=0;
    print g, v
  }' "$TPM_CSV" > "${OUTDIR}/tpm_${cname}.map"

  # join TPM onto transcript BED -> raw (possibly overlapping) bedGraph
  # (accept map split on tabs or spaces to be robust)
  awk 'BEGIN{FS=OFS="\t"}
       NR==FNR {
         n=split($0,a,/[\t ]+/); if(n>=2) val[a[1]]=a[2]; next
       }
       { g=$4; v=(g in val ? val[g] : 0); print $1,$2,$3,v }
  ' "${OUTDIR}/tpm_${cname}.map" "$BED_GENE" \
  | LC_ALL=C sort -k1,1 -k2,2n > "${OUTDIR}/${cname}.raw.bedGraph"

  # --- resolve overlaps by disjoining; value = max across overlaps ---
  awk 'BEGIN{OFS="\t"}
    function flush(chr,   i,nBp,prev,segS,segE,maxv,k,covered){
      if (chr=="") return
      nBp=bpN[chr]
      if (nBp<2) { cleanup(chr); return }
      prev=bp[chr,1]
      for (i=2; i<=nBp; i++) {
        segS=prev; segE=bp[chr,i]
        if (segE>segS) {
          maxv=0; covered=0
          nI=idxN[chr]
          for (k=1; k<=nI; k++) {
            if (S[chr,k] < segE && E[chr,k] > segS) {
              covered=1
              if (V[chr,k] > maxv) maxv=V[chr,k]
            }
          }
          if (covered) print chr, segS, segE, maxv
        }
        prev=segE
      }
      cleanup(chr)
    }
    function cleanup(chr, i){
      delete bpN[chr]; delete idxN[chr]
      for (i in bp) if (i ~ "^"chr SUBSEP) delete bp[i]
      for (i in S)  if (i ~ "^"chr SUBSEP) { delete S[i]; delete E[i]; delete V[i] }
    }
    {
      c=$1; s=$2; e=$3; v=$4+0
      if (c!=curChr && curChr!="") { flush(curChr) }
      curChr=c
      ++bpN[c]; bp[c,bpN[c]]=s
      ++bpN[c]; bp[c,bpN[c]]=e
      if (!(c in idxN)) idxN[c]=0
      idxN[c]++; i=idxN[c]; S[c,i]=s; E[c,i]=e; V[c,i]=v
    }
    END{ flush(curChr) }' "${OUTDIR}/${cname}.raw.bedGraph" \
  | LC_ALL=C sort -k1,1 -k2,2n > "${OUTDIR}/${cname}.bedGraph"

  # quick sanity stats
  rawN=$(wc -l < "${OUTDIR}/${cname}.raw.bedGraph" || echo 0)
  disN=$(wc -l < "${OUTDIR}/${cname}.bedGraph" || echo 0)
  echo "[info] ${cname}: raw rows=${rawN}, disjoint rows=${disN}"

  # enforce non-overlap & merge adjacent same-valued tiles
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "${OUTDIR}/${cname}.bedGraph" \
  | awk 'BEGIN{OFS="\t"}
    NR==1 { pc=$1; ps=$2; pe=$3; pv=$4; next }
    {
      c=$1; s=$2; e=$3; v=$4
      if (c!=pc) {
        if (pe>ps) print pc,ps,pe,pv
        pc=c; ps=s; pe=e; pv=v; next
      }
      if (s<pe) s=pe
      if (s>=e) next
      if (s==pe && v==pv) { pe=e } else { if (pe>ps) print pc,ps,pe,pv; ps=s; pe=e; pv=v }
    }
    END { if (pe>ps) print pc,ps,pe,pv }' \
  > "${OUTDIR}/${cname}.bedGraph.fixed"

  mv "${OUTDIR}/${cname}.bedGraph.fixed" "${OUTDIR}/${cname}.bedGraph"
  disN=$(wc -l < "${OUTDIR}/${cname}.bedGraph" || echo 0)
  echo "[info] ${cname}: rows after fix=${disN}"

  command -v bedSort >/dev/null 2>&1 && bedSort "${OUTDIR}/${cname}.bedGraph" "${OUTDIR}/${cname}.bedGraph"

  [[ "$disN" -gt 0 ]] || { echo "[ERR] ${cname}.bedGraph is empty after fixes." >&2; exit 1; }

  bedGraphToBigWig "${OUTDIR}/${cname}.bedGraph" "$CHROMSIZES" "${OUTDIR}/${cname}.tpm.bw"
  echo "[ok] ${OUTDIR}/${cname}.tpm.bw"
done

echo "[done] TPM bigWigs in ${OUTDIR}/"
