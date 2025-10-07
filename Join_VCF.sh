#!/bin/bash
# Joint variant calling pipeline for Medicago genomes
# Aligns each genome to Reference genome.fna, calls SNPs,
# and merges them into a single joint VCF.

# Requirements: minimap2, samtools, bcftools

# Reference genome
REF="Reference genome.fna"

# Index the reference (only once)
samtools faidx "$REF"

# List of query genomes (exclude Reference)
GENOMES=(
  "medtr.HM020.gnm1.PNGN.genome_main.fna"
  "medtr.HM010.gnm1.YDPY.genome_main.fna"
  "medtr.HM006.gnm1.XKYH.genome_main.fna"
  "medtr.HM005.gnm1.BQQH.genome_main.fna"
  "medtr.HM004.gnm1.L7PM.genome_main.fna"
  "medtr.HM002.gnm1.TWKJ.genome_main.fna"
)

VCF_LIST=()

# Step 1: Per-genome variant calling
for GENOME in "${GENOMES[@]}"; do
    PREFIX=$(basename "$GENOME" .fna)
    echo ">>> Processing $GENOME"

    # Align genome to reference
    minimap2 -ax asm5 "$REF" "$GENOME" > "${PREFIX}.sam"

    # Convert to BAM, sort and index
    samtools view -bS "${PREFIX}.sam" | samtools sort -o "${PREFIX}.sorted.bam"
    samtools index "${PREFIX}.sorted.bam"
    rm "${PREFIX}.sam"  # save space

    # Variant calling
    bcftools mpileup -Ou -f "$REF" "${PREFIX}.sorted.bam" | \
      bcftools call -mv -Oz -o "${PREFIX}.vcf.gz"

    # Index VCF
    bcftools index "${PREFIX}.vcf.gz"

    VCF_LIST+=("${PREFIX}.vcf.gz")

    echo ">>> Finished $GENOME → ${PREFIX}.vcf.gz"
done

# Step 2: Merge all VCFs into a joint multi-sample VCF
echo ">>> Merging all VCFs into joint file..."
bcftools merge -Oz -o Joint_Medicago.vcf.gz "${VCF_LIST[@]}"

# Index the joint VCF
bcftools index Joint_Medicago.vcf.gz

echo ">>> Joint VCF created: Joint_Medicago.vcf.gz"
