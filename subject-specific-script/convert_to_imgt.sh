#!/usr/bin/env bash
# Convert VDJbase haplotype allele CSVs into IMGT-formatted reference FASTAs
# for nf-core/airrflow (--reference_fasta)
#
# Usage: bash convert_to_imgt.sh <TRA.csv> <TRB.csv> <TRD.csv> <TRG.csv> [output_dir]
#
# output_dir defaults to "custom_imgtdb".  Pass "." to write directly into the
# current working directory (the script will zip up human/ + IMGT.yaml only).
#
# Example (writing into the current directory):
#   bash convert_to_imgt.sh TRA.csv TRB.csv TRD.csv TRG.csv .

set -euo pipefail

TRA_CSV="${1}"
TRB_CSV="${2}"
TRD_CSV="${3}"
TRG_CSV="${4}"
OUTDIR="${5:-custom_imgtdb}"
SPECIES="Homo_sapiens"
TODAY=$(date +%Y-%m-%d)

# Resolve OUTDIR to an absolute path so ZIP logic is unambiguous
OUTDIR="$(cd "$OUTDIR" 2>/dev/null && pwd || { mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd; })"

# Output directories
VDJ="${OUTDIR}/human/vdj"
LEAD="${OUTDIR}/human/leader"
CONST="${OUTDIR}/human/constant"
mkdir -p "$VDJ" "$LEAD" "$CONST"

# --- helper: write one IMGT FASTA record ---
# Usage: write_record <allele> <region_type> <sequence> <outfile>
write_record() {
    local allele="$1"
    local region="$2"
    local seq="$3"
    local outfile="$4"

    [ -z "$seq" ] && return   # skip empty sequences

    # count ungapped length (remove dots)
    local ungapped="${seq//./}"
    local len=${#ungapped}
    local gapped_len=${#seq}
    local extra=$(( gapped_len - len ))

    # IMGT header
    printf ">CUSTOM|%s|%s|F|%s|.|%d nt|1| | | | |%d+%d=%d| | |\n" \
        "$allele" "$SPECIES" "$region" "$len" "$len" "$extra" "$gapped_len" >> "$outfile"

    # sequence wrapped at 60 chars (fold)
    echo "$seq" | tr '[:upper:]' '[:lower:]' | fold -w 60 >> "$outfile"
}

# --- helper: process one CSV ---
# Reads CSV, deduplicates by allele, routes each row to the right FASTA
process_csv() {
    local csv="$1"
    local chain="$2"   # TRA, TRB, TRD, or TRG

    # Parse CSV with awk:
    # - Row 1 is the header; build a column-index map
    # - Deduplicate on vdjbase_allele (col index stored in idx["vdjbase_allele"])
    # - For each unique allele, extract the right columns and call write_record
    #
    # Because sequences can contain commas inside quotes, we use a simple approach:
    # the sequences in these CSVs do NOT have commas, so plain FS="," is safe.

    awk -F',' -v chain="$chain" -v vdj="$VDJ" -v lead="$LEAD" -v const_dir="$CONST" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            # strip quotes/spaces from header names
            gsub(/^[ "]+|[ "]+$/, "", $i)
            col[$i] = i
        }
        next
    }
    {
        allele = $col["vdjbase_allele"]
        gsub(/^[ "]+|[ "]+$/, "", allele)
        if (allele == "" || seen[allele]++) next

        gene = $col["gene"]
        gsub(/^[ "]+|[ "]+$/, "", gene)

        vrg    = $col["V-REGION-GAPPED"];   gsub(/^[ "]+|[ "]+$/, "", vrg)
        lp1    = $col["L-PART1"];            gsub(/^[ "]+|[ "]+$/, "", lp1)
        lp2    = $col["L-PART2"];            gsub(/^[ "]+|[ "]+$/, "", lp2)
        jreg   = $col["J-REGION"];           gsub(/^[ "]+|[ "]+$/, "", jreg)
        dreg   = $col["D-REGION"];           gsub(/^[ "]+|[ "]+$/, "", dreg)
        creg   = $col["C-REGION"];           gsub(/^[ "]+|[ "]+$/, "", creg)

        # Determine gene family from prefix
        if (gene ~ /^TRA?V/ || gene ~ /^TRBV/ || gene ~ /^TRDV/ || gene ~ /^TRGV/) {
            # V gene → vdj (gapped) + leader
            if (vrg != "") {
                ungapped = vrg; gsub(/\./, "", ungapped)
                len = length(ungapped); glen = length(vrg); extra = glen - len
                printf ">CUSTOM|%s|%s|F|V-REGION|.|%d nt|1| | | | |%d+%d=%d| | |\n",
                    allele, "Homo_sapiens", len, len, extra, glen \
                    >> (vdj "/imgt_human_" chain "V.fasta")
                seq = tolower(vrg)
                while (length(seq) > 0) {
                    print substr(seq, 1, 60) >> (vdj "/imgt_human_" chain "V.fasta")
                    seq = substr(seq, 61)
                }
            }
            lseq = lp1 lp2
            if (lseq != "") {
                llen = length(lseq)
                printf ">CUSTOM|%s|%s|F|L-PART1+L-PART2|.|%d nt|1| | | | |%d+0=%d| | |\n",
                    allele, "Homo_sapiens", llen, llen, llen \
                    >> (lead "/imgt_human_" chain "L.fasta")
                seq = tolower(lseq)
                while (length(seq) > 0) {
                    print substr(seq, 1, 60) >> (lead "/imgt_human_" chain "L.fasta")
                    seq = substr(seq, 61)
                }
            }
        }
        else if (gene ~ /^TRA?J/ || gene ~ /^TRBJ/ || gene ~ /^TRDJ/ || gene ~ /^TRGJ/) {
            if (jreg != "") {
                jlen = length(jreg)
                printf ">CUSTOM|%s|%s|F|J-REGION|.|%d nt|1| | | | |%d+0=%d| | |\n",
                    allele, "Homo_sapiens", jlen, jlen, jlen \
                    >> (vdj "/imgt_human_" chain "J.fasta")
                seq = tolower(jreg)
                while (length(seq) > 0) {
                    print substr(seq, 1, 60) >> (vdj "/imgt_human_" chain "J.fasta")
                    seq = substr(seq, 61)
                }
            }
        }
        else if (gene ~ /^TRBD/ || gene ~ /^TRDD/) {
            if (dreg != "") {
                dlen = length(dreg)
                # Use chain-specific D file (TRBD or TRDD)
                dfile = vdj "/imgt_human_" chain "D.fasta"
                printf ">CUSTOM|%s|%s|F|D-REGION|.|%d nt|1| | | | |%d+0=%d| | |\n",
                    allele, "Homo_sapiens", dlen, dlen, dlen >> dfile
                seq = tolower(dreg)
                while (length(seq) > 0) {
                    print substr(seq, 1, 60) >> dfile
                    seq = substr(seq, 61)
                }
            }
        }
        else if (gene ~ /^TRA?C/ || gene ~ /^TRBC/ || gene ~ /^TRDC/ || gene ~ /^TRGC/) {
            if (creg != "") {
                clen = length(creg)
                printf ">CUSTOM|%s|%s|F|EX1+EX2+EX3|.|%d nt|1| | | | |%d+0=%d| | |\n",
                    allele, "Homo_sapiens", clen, clen, clen \
                    >> (const_dir "/imgt_human_" chain "C.fasta")
                seq = tolower(creg)
                while (length(seq) > 0) {
                    print substr(seq, 1, 60) >> (const_dir "/imgt_human_" chain "C.fasta")
                    seq = substr(seq, 61)
                }
            }
        }
    }
    ' "$csv"
}

# --- run ---
echo "Processing TRA..."
process_csv "$TRA_CSV" "TRA"

echo "Processing TRB..."
process_csv "$TRB_CSV" "TRB"

echo "Processing TRD..."
process_csv "$TRD_CSV" "TRD"

echo "Processing TRG..."
process_csv "$TRG_CSV" "TRG"

# --- IMGT.yaml ---
cat > "${OUTDIR}/IMGT.yaml" << YAML
source:  custom haplotype-resolved alleles (VDJbase)
date:    ${TODAY}
species:
    - human:Homo+sapiens
YAML

# --- report ---
echo ""
echo "Output files:"
for f in "$VDJ"/*.fasta "$LEAD"/*.fasta "$CONST"/*.fasta; do
    count=$(grep -c "^>" "$f" 2>/dev/null || echo 0)
    printf "  %3d records  %s\n" "$count" "$f"
done

# --- zip ---
# nf-core/airrflow UNZIP_REFERENCE_FASTA runs "unzip custom_imgtdb.zip" and then
# expects to find a directory called exactly "custom_imgtdb" in its working dir.
# So the archive must have a single top-level folder named "custom_imgtdb/".
#
# Strategy: build a temporary staging dir next to OUTDIR, symlink (or copy) the
# content into custom_imgtdb/, zip from the staging dir, then clean up.

ZIP_OUT="${OUTDIR}/custom_imgtdb.zip"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

mkdir -p "${STAGE}/custom_imgtdb"
cp -r "${OUTDIR}/human"      "${STAGE}/custom_imgtdb/"
cp    "${OUTDIR}/IMGT.yaml"  "${STAGE}/custom_imgtdb/"

(cd "$STAGE" && zip -qr "$ZIP_OUT" custom_imgtdb/)

echo ""
echo "✓ ZIP → ${ZIP_OUT}"
echo "  Archive layout: custom_imgtdb/human/... (as required by nf-core/airrflow)"
echo "  Use with: --reference_fasta ${ZIP_OUT}"