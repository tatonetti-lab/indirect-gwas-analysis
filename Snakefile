rule all:
  input:
    # Inputs
    multiext("data/geno/hapmap3_variants_white_british_100k", ".pgen", ".psam", ".pvar"),
    # Outputs
    multiext("data/geno/geno_500k", ".pgen", ".psam", ".pvar"),

rule filter_genotypes:
  input:
    geno = multiext("data/geno/hapmap3_variants_white_british_100k", ".pgen", ".psam", ".pvar"),
  output:
    multiext("data/geno/geno_500k", ".pgen", ".psam", ".pvar"),
  params:
    n_variants = 500_000,
    input_prefix = "data/geno/hapmap3_variants_white_british_100k",
    output_prefix = "data/geno/geno_500k"
  shell:
    """
    plink2 \
      --pfile {params.input_prefix} \
      --thin-count {params.n_variants} \
      --make-pgen \
      --out {params.output_prefix}
    """
