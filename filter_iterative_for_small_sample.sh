#!/bin/bash
set -euo pipefail

VCF="Joint_Medicago.vcf.gz"
OUT_PREFIX="Joint_Medicago.filtered"
MIN_QUAL1=10
MIN_DP1=1
# require present in at least 2 samples: 2/6 = 0.333 -> max-missing = 0.34
MAX_MISSING_2SAMPLES=0.34
MAC1=1    # keep minor allele count >=1 (we will filter by presence via max-missing)
MIN_ALLELES=2
MAX_ALLELES=2

# Relaxed settings to allow singletons (present in 1/6 samples): 1/6 = 0.1667
MAX_MISSING_SINGLETON=0.17

# helper to print counts
site_count() {
  local f=$1
  if [[ ! -f "$f" ]]; then
    echo "0"
    return
  fi
  zgrep -v '^#' "$f" | wc -l
}

echo "1) Quick VCF stats (bcftools/vcftools checks)"
echo "Total sites in VCF:"
bcftools view -h "$VCF" >/dev/null 2>&1 || { echo "ERROR: cannot read $VCF"; exit 1; }
bcftools query -f '%CHROM\t%POS\n' "$VCF" | wc -l || true
echo

echo "Per-sample missingness (vcftools):"
vcftools --gzvcf "$VCF" --missing-indv --out prefilter_missing || true
cat prefilter_missing.imiss || true
echo

echo "Site mean depth (vcftools):"
vcftools --gzvcf "$VCF" --site-mean-depth --out prefilter_depth || true
echo "Preview site mean depth (first 10 lines):"
head -n 12 prefilter_depth | sed -n '1,12p'
echo

# ---------- First filter: require presence in >=2 samples ----------
echo "2) Attempting first (recommended) filtering: require site present in >=2 samples (max-missing ${MAX_MISSING_2SAMPLES})"
vcftools --gzvcf "$VCF" \
  --minQ ${MIN_QUAL1} \
  --min-meanDP ${MIN_DP1} \
  --max-missing ${MAX_MISSING_2SAMPLES} \
  --mac ${MAC1} \
  --remove-indels \
  --min-alleles ${MIN_ALLELES} --max-alleles ${MAX_ALLELES} \
  --recode --recode-INFO-all \
  --out "${OUT_PREFIX}.2plus"

echo "Sites retained (>=2 samples):"
SITES_2PLUS=$(site_count "${OUT_PREFIX}.2plus.recode.vcf")
echo "${SITES_2PLUS}"

if [[ "${SITES_2PLUS}" -gt 0 ]]; then
  echo "Using ${OUT_PREFIX}.2plus.recode.vcf as final filtered VCF."
  FINAL_VCF="${OUT_PREFIX}.2plus.recode.vcf"
else
  echo "No sites retained with >=2-sample presence filter. Relaxing filters to keep singletons."
  # ---------- Relaxed filter: allow singletons (present in >=1 sample) ----------
  vcftools --gzvcf "$VCF" \
    --minQ ${MIN_QUAL1} \
    --min-meanDP ${MIN_DP1} \
    --max-missing ${MAX_MISSING_SINGLETON} \
    --mac ${MAC1} \
    --remove-indels \
    --min-alleles ${MIN_ALLELES} --max-alleles ${MAX_ALLELES} \
    --recode --recode-INFO-all \
    --out "${OUT_PREFIX}.singleton"

  echo "Sites retained (singletons allowed):"
  SITES_SINGLE=$(site_count "${OUT_PREFIX}.singleton.recode.vcf")
  echo "${SITES_SINGLE}"

  if [[ "${SITES_SINGLE}" -gt 0 ]]; then
    echo "Using ${OUT_PREFIX}.singleton.recode.vcf as final filtered VCF."
    FINAL_VCF="${OUT_PREFIX}.singleton.recode.vcf"
  else
    echo "Still zero sites after relaxed filters. Next step: allow multiallelic sites and lower QUAL to 5."
    # ---------- More relaxed: allow multiallelic and very low QUAL ----------
    vcftools --gzvcf "$VCF" \
      --minQ 5 \
      --min-meanDP 1 \
      --max-missing ${MAX_MISSING_SINGLETON} \
      --mac 1 \
      --recode --recode-INFO-all \
      --out "${OUT_PREFIX}.very_relaxed"

    echo "Sites retained (very relaxed):"
    SITES_VREL=$(site_count "${OUT_PREFIX}.very_relaxed.recode.vcf")
    echo "${SITES_VREL}"

    if [[ "${SITES_VREL}" -gt 0 ]]; then
      FINAL_VCF="${OUT_PREFIX}.very_relaxed.recode.vcf"
      echo "Using ${FINAL_VCF} as final VCF (very relaxed)."
    else
      echo "No variants retained under multiple relaxed settings. There may be a problem upstream (variant calling). Exiting."
      exit 2
    fi
  fi
fi

# ---------- Final diagnostics ----------
echo; echo "FINAL VCF: ${FINAL_VCF}"
echo "Variant count:"
zcat "${FINAL_VCF}" | grep -v '^#' | wc -l

echo "Per-sample missingness on final VCF:"
vcftools --vcf "${FINAL_VCF}" --missing-indv --out final_missing || true
cat final_missing.imiss || true

echo "Site mean depth on final VCF (sample):"
vcftools --vcf "${FINAL_VCF}" --site-mean-depth --out final_depth || true
head -n 12 final_depth || true

echo "Generate freq file (SFS) for final VCF:"
vcftools --vcf "${FINAL_VCF}" --freq --out final_freq || true
head -n 12 final_freq || true

# Optional: convert to PLINK (if sites > 0)
SITES_FINAL=$(zcat "${FINAL_VCF}" | grep -v '^#' | wc -l)
if [[ "${SITES_FINAL}" -gt 0 ]]; then
  echo "Converting final VCF to PLINK (for PCA/fastSTRUCTURE)..."
  plink --vcf "${FINAL_VCF}" --double-id --allow-extra-chr --make-bed --out final_plink || true
  echo "PLINK conversion done: final_plink.*"
else
  echo "No sites to convert to PLINK. Stop here."
fi

echo "Filtering script completed."
