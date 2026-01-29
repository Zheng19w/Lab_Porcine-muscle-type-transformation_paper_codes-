#!/usr/bin/env bash
set -euo pipefail

# =======================
# Inputs (EDIT THESE)
# =======================
BEDTOOLS=bedtools         # v2.30.0
BCFTOOLS=bcftools

# consensus peaks (BED; 0-based, half-open)
LEAN_H3K27AC="peaks/LEAN.H3K27ac.consensus.bed"
LEAN_H3K4ME3="peaks/LEAN.H3K4me3.consensus.bed"
CHINESE_H3K27AC="peaks/CHINESE_INDIGENOUS.H3K27ac.consensus.bed"
CHINESE_H3K4ME3="peaks/CHINESE_INDIGENOUS.H3K4me3.consensus.bed"

# cohort PASS VCFs (biallelic recommended)
VCF_SNPS="joint_calling/pass/cohort.snps.PASS.vcf.gz"
VCF_INDELS="joint_calling/pass/cohort.indels.PASS.vcf.gz"

# sample lists (one sample ID per line; must match VCF column names)
LEANTYPE_SAMPLES="lists/lean_type.samples"                 # LargeWhite + Duroc
CHINESE_INDIGENOUS_SAMPLES="lists/chinese_indigenous.samples"  # JianLi + LaiWu

# AF threshold to define "present in group" for intersection (tune as needed)
AF_PRESENT="0.05"

# output
OUT="regvar"
mkdir -p "${OUT}"/{reg,vars,overlap,summary,tmp}

# =======================
# Helpers
# =======================
sort_merge() {
  # $1 in.bed  $2 out.bed
  "${BEDTOOLS}" sort -i "$1" | "${BEDTOOLS}" merge -i - > "$2"
}

bed_n() { [[ -s "$1" ]] && wc -l < "$1" || echo 0; }

vcf_to_bed_af() {
  # $1 vcf.gz  $2 sample_list  $3 af_threshold  $4 out.bed
  "${BCFTOOLS}" view -S "$2" -Ou "$1" \
    | "${BCFTOOLS}" +fill-tags -Ou -- -t AF \
    | "${BCFTOOLS}" view -i "INFO/AF>=${3}" -Ou \
    | "${BCFTOOLS}" query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' \
    | LC_ALL=C sort -k1,1 -k2,2n > "$4"
}

# =======================
# 1) Define active regulatory regions
# Promoters: H3K4me3 peaks overlapping H3K27ac by >=50%
# Enhancers: H3K27ac peaks not overlapping H3K4me3
# =======================
# sort/merge inputs
sort_merge "${LEAN_H3K27AC}" "${OUT}/tmp/lean_type.H3K27ac.merged.bed"
sort_merge "${LEAN_H3K4ME3}" "${OUT}/tmp/lean_type.H3K4me3.merged.bed"
sort_merge "${CHINESE_H3K27AC}" "${OUT}/tmp/chinese_indigenous.H3K27ac.merged.bed"
sort_merge "${CHINESE_H3K4ME3}" "${OUT}/tmp/chinese_indigenous.H3K4me3.merged.bed"

# promoters
"${BEDTOOLS}" intersect -a "${OUT}/tmp/lean_type.H3K4me3.merged.bed" -b "${OUT}/tmp/lean_type.H3K27ac.merged.bed" -f 0.5 -u \
  | "${BEDTOOLS}" sort -i - | "${BEDTOOLS}" merge -i - > "${OUT}/reg/lean_type.promoter.bed"

"${BEDTOOLS}" intersect -a "${OUT}/tmp/chinese_indigenous.H3K4me3.merged.bed" -b "${OUT}/tmp/chinese_indigenous.H3K27ac.merged.bed" -f 0.5 -u \
  | "${BEDTOOLS}" sort -i - | "${BEDTOOLS}" merge -i - > "${OUT}/reg/chinese_indigenous.promoter.bed"

# enhancers
"${BEDTOOLS}" intersect -a "${OUT}/tmp/lean_type.H3K27ac.merged.bed" -b "${OUT}/tmp/lean_type.H3K4me3.merged.bed" -v \
  | "${BEDTOOLS}" sort -i - | "${BEDTOOLS}" merge -i - > "${OUT}/reg/lean_type.enhancer.bed"

"${BEDTOOLS}" intersect -a "${OUT}/tmp/chinese_indigenous.H3K27ac.merged.bed" -b "${OUT}/tmp/chinese_indigenous.H3K4me3.merged.bed" -v \
  | "${BEDTOOLS}" sort -i - | "${BEDTOOLS}" merge -i - > "${OUT}/reg/chinese_indigenous.enhancer.bed"

# common vs specific (reciprocal >=50% overlap)
"${BEDTOOLS}" intersect -a "${OUT}/reg/lean_type.promoter.bed" -b "${OUT}/reg/chinese_indigenous.promoter.bed" -f 0.5 -r -u \
  | "${BEDTOOLS}" sort -i - | "${BEDTOOLS}" merge -i - > "${OUT}/reg/ACP.bed"   # common promoters
"${BEDTOOLS}" intersect -a "${OUT}/reg/lean_type.enhancer.bed" -b "${OUT}/reg/chinese_indigenous.enhancer.bed" -f 0.5 -r -u \
  | "${BEDTOOLS}" sort -i - | "${BEDTOOLS}" merge -i - > "${OUT}/reg/ACE.bed"   # common enhancers

"${BEDTOOLS}" intersect -a "${OUT}/reg/lean_type.promoter.bed" -b "${OUT}/reg/ACP.bed" -v > "${OUT}/reg/LSP.bed"
"${BEDTOOLS}" intersect -a "${OUT}/reg/lean_type.enhancer.bed" -b "${OUT}/reg/ACE.bed" -v > "${OUT}/reg/LSE.bed"
"${BEDTOOLS}" intersect -a "${OUT}/reg/chinese_indigenous.promoter.bed" -b "${OUT}/reg/ACP.bed" -v > "${OUT}/reg/CSP.bed"
"${BEDTOOLS}" intersect -a "${OUT}/reg/chinese_indigenous.enhancer.bed" -b "${OUT}/reg/ACE.bed" -v > "${OUT}/reg/CSE.bed"

# =======================
# 2) SNP/Indel intersection between breed groups
# =======================
"${BCFTOOLS}" query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "${VCF_SNPS}"   | LC_ALL=C sort -k1,1 -k2,2n > "${OUT}/vars/ALL.SNP.bed"
"${BCFTOOLS}" query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "${VCF_INDELS}" | LC_ALL=C sort -k1,1 -k2,2n > "${OUT}/vars/ALL.INDEL.bed"

vcf_to_bed_af "${VCF_SNPS}"   "${LEANTYPE_SAMPLES}"           "${AF_PRESENT}" "${OUT}/vars/lean_type.SNP.present.bed"
vcf_to_bed_af "${VCF_SNPS}"   "${CHINESE_INDIGENOUS_SAMPLES}" "${AF_PRESENT}" "${OUT}/vars/chinese_indigenous.SNP.present.bed"
vcf_to_bed_af "${VCF_INDELS}" "${LEANTYPE_SAMPLES}"           "${AF_PRESENT}" "${OUT}/vars/lean_type.INDEL.present.bed"
vcf_to_bed_af "${VCF_INDELS}" "${CHINESE_INDIGENOUS_SAMPLES}" "${AF_PRESENT}" "${OUT}/vars/chinese_indigenous.INDEL.present.bed"

# SNPs
"${BEDTOOLS}" intersect -a "${OUT}/vars/lean_type.SNP.present.bed" -b "${OUT}/vars/chinese_indigenous.SNP.present.bed" -v > "${OUT}/vars/LSS.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/chinese_indigenous.SNP.present.bed" -b "${OUT}/vars/lean_type.SNP.present.bed" -v > "${OUT}/vars/CSS.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/lean_type.SNP.present.bed" -b "${OUT}/vars/chinese_indigenous.SNP.present.bed" -u > "${OUT}/vars/ACS.bed"

# INDELs
"${BEDTOOLS}" intersect -a "${OUT}/vars/lean_type.INDEL.present.bed" -b "${OUT}/vars/chinese_indigenous.INDEL.present.bed" -v > "${OUT}/vars/LSI.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/chinese_indigenous.INDEL.present.bed" -b "${OUT}/vars/lean_type.INDEL.present.bed" -v > "${OUT}/vars/CSI.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/lean_type.INDEL.present.bed" -b "${OUT}/vars/chinese_indigenous.INDEL.present.bed" -u > "${OUT}/vars/ACI.bed"

# =======================
# 3) Variants within active regulatory elements
# =======================
"${BEDTOOLS}" intersect -a "${OUT}/vars/LSS.bed" -b "${OUT}/reg/LSP.bed" -u > "${OUT}/overlap/LSS_in_LSP.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/LSS.bed" -b "${OUT}/reg/LSE.bed" -u > "${OUT}/overlap/LSS_in_LSE.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/LSI.bed" -b "${OUT}/reg/LSP.bed" -u > "${OUT}/overlap/LSI_in_LSP.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/LSI.bed" -b "${OUT}/reg/LSE.bed" -u > "${OUT}/overlap/LSI_in_LSE.bed"

"${BEDTOOLS}" intersect -a "${OUT}/vars/CSS.bed" -b "${OUT}/reg/CSP.bed" -u > "${OUT}/overlap/CSS_in_CSP.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/CSS.bed" -b "${OUT}/reg/CSE.bed" -u > "${OUT}/overlap/CSS_in_CSE.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/CSI.bed" -b "${OUT}/reg/CSP.bed" -u > "${OUT}/overlap/CSI_in_CSP.bed"
"${BEDTOOLS}" intersect -a "${OUT}/vars/CSI.bed" -b "${OUT}/reg/CSE.bed" -u > "${OUT}/overlap/CSI_in_CSE.bed"

# =======================
# 4) Minimal enrichment summary
# =======================
summary="${OUT}/summary/enrichment.tsv"
echo -e "set\ttype\treg\tset_total\tset_in_reg\tbg_total\tbg_in_reg\tfold_enrichment" > "${summary}"

calc_enrich() {
  local set_name="$1" typ="$2" set_bed="$3" reg_name="$4" reg_bed="$5" bg_bed="$6"
  local set_total bg_total set_in bg_in fold
  set_total=$(bed_n "${set_bed}")
  bg_total=$(bed_n "${bg_bed}")
  set_in=$("${BEDTOOLS}" intersect -a "${set_bed}" -b "${reg_bed}" -u | wc -l | tr -d ' ')
  bg_in=$("${BEDTOOLS}" intersect -a "${bg_bed}"  -b "${reg_bed}" -u | wc -l | tr -d ' ')
  fold=$(awk -v a="${set_in}" -v b="${set_total}" -v c="${bg_in}" -v d="${bg_total}" 'BEGIN{
      if(b==0||d==0||c==0) {print "NA"; exit}
      printf "%.6f", (a/b)/(c/d)
    }')
  echo -e "${set_name}\t${typ}\t${reg_name}\t${set_total}\t${set_in}\t${bg_total}\t${bg_in}\t${fold}" >> "${summary}"
}

calc_enrich "LSS" "SNP"   "${OUT}/vars/LSS.bed" "LSP" "${OUT}/reg/LSP.bed" "${OUT}/vars/ALL.SNP.bed"
calc_enrich "LSS" "SNP"   "${OUT}/vars/LSS.bed" "LSE" "${OUT}/reg/LSE.bed" "${OUT}/vars/ALL.SNP.bed"
calc_enrich "LSI" "INDEL" "${OUT}/vars/LSI.bed" "LSP" "${OUT}/reg/LSP.bed" "${OUT}/vars/ALL.INDEL.bed"
calc_enrich "LSI" "INDEL" "${OUT}/vars/LSI.bed" "LSE" "${OUT}/reg/LSE.bed" "${OUT}/vars/ALL.INDEL.bed"

calc_enrich "CSS" "SNP"   "${OUT}/vars/CSS.bed" "CSP" "${OUT}/reg/CSP.bed" "${OUT}/vars/ALL.SNP.bed"
calc_enrich "CSS" "SNP"   "${OUT}/vars/CSS.bed" "CSE" "${OUT}/reg/CSE.bed" "${OUT}/vars/ALL.SNP.bed"
calc_enrich "CSI" "INDEL" "${OUT}/vars/CSI.bed" "CSP" "${OUT}/reg/CSP.bed" "${OUT}/vars/ALL.INDEL.bed"
calc_enrich "CSI" "INDEL" "${OUT}/vars/CSI.bed" "CSE" "${OUT}/reg/CSE.bed" "${OUT}/vars/ALL.INDEL.bed"

echo "[OK] Outputs:"
echo "  Regulatory: ${OUT}/reg/{LSP,LSE,CSP,CSE,ACP,ACE}.bed"
echo "  Variants:   ${OUT}/vars/{LSS,LSI,CSS,CSI,ACS,ACI}.bed"
echo "  Overlaps:   ${OUT}/overlap/*.bed"
echo "  Summary:    ${OUT}/summary/enrichment.tsv"