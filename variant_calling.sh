#!/bin/bash
# Script to align Medicago genomes against Reference genome
# and call variants (SNPs/indels).
# Requirements: minimap2, samtools, bcftools installed

# Reference genome
REF="Reference genome.fna"

# Index the reference (only once)
samtools faidx "$REF"

# List of query genomes
GENOMES=(
  "medtr.HM020.gnm1.PNGN.genome_main.fna"
  "medtr.HM010.gnm1.YDPY.genome_main.fna"
  "medtr.HM006.gnm1.XKYH.genome_main.fna"
  "medtr.HM005.gnm1.BQQH.genome_main.fna"
  "medtr.HM004.gnm1.L7PM.genome_main.fna"
  "medtr.HM002.gnm1.TWKJ.genome_main.fna"
)

# Loop through genomes
for GENOME in "${GENOMES[@]}"; do
    # Get prefix (file basename without extension)
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

    echo ">>> Finished $GENOME → ${PREFIX}.vcf.gz"
done

echo "All genomes processed successfully."
