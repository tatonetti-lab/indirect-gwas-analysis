import polars as pl


phenotype_names = (
    pl.read_csv("data/pheno/pheno_jan2024.tsv", separator="\t", n_rows=0)
    .drop(["#FID", "IID"])
    .columns
)


rule all:
  input:
    # Inputs
    multiext("data/geno/hapmap3_variants_white_british_100k", ".pgen", ".psam", ".pvar"),
    "data/pheno/pheno_jan2024.tsv",
    "data/pheno/covar.tsv",
    # Outputs
    multiext("data/geno/geno_500k", ".pgen", ".psam", ".pvar"),
    multiext("data/grm/gcta.grm", ".bin", ".id", ".N.bin", ".sp"),
    expand(
        "data/gwas/plink.{phenotype}.glm.linear.zst",
        phenotype=phenotype_names
    ),



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


rule filter_binary_pheno:
  input:
    raw_pheno = "data/pheno/pheno_jan2024.tsv",
    psam = "data/geno/geno_500k.psam",
  output:
    "data/pheno/binary_pheno.tsv",
  run:
    psam_df = pl.scan_csv(input.psam, separator="\t").select("IID")
    (
        pl.scan_csv(input.raw_pheno, separator="\t")
        .select(pl.col("#FID").alias("FID"), "IID", "^b_.*$")
        .join(psam_df, on="IID")
        .sink_csv(output[0], separator="\t")
    )


rule gcta_grm_part:
  input:
    geno = multiext("data/geno/geno_500k", ".pgen", ".psam", ".pvar"),
  output:
    multiext("data/grm/gcta.part_10_{i}.grm", ".bin", ".id", ".N.bin"),
  params:
    geno_prefix = "data/geno/geno_500k",
    output_prefix = "data/grm/gcta",
  threads: 35
  benchmark: "benchmarks/gcta_grm_part_{i}.txt"
  shell:
    """
    gcta \
      --pfile {params.geno_prefix} \
      --make-grm-part 10 {wildcards.i} \
      --thread-num {threads} \
      --out {params.output_prefix}
    """


rule gcta_collect_grm:
  input:
    collect(
      "data/grm/gcta.part_10_{i}.grm.{ext}",
      i=[f"{i:02}" for i in range(1, 11)],
      ext=["bin", "id", "N.bin"]
    ),
  output:
    bin = "data/grm/gcta.grm.bin",
    id = "data/grm/gcta.grm.id",
    N = "data/grm/gcta.grm.N.bin",
  benchmark: "benchmarks/gcta_collect_grm.txt"
  params:
    bin_files = " ".join(f"data/grm/gcta.part_10_{i:02}.grm.bin" for i in range(1, 11)),
    id_files = " ".join(f"data/grm/gcta.part_10_{i:02}.grm.id" for i in range(1, 11)),
    N_files = " ".join(f"data/grm/gcta.part_10_{i:02}.grm.N.bin" for i in range(1, 11))
  shell:
    """
    cat {params.bin_files} > {output.bin}
    cat {params.id_files} > {output.id}
    cat {params.N_files} > {output.N}
    """


rule gcta_sparsify_grm:
  input:
    multiext("data/grm/gcta.grm", ".bin", ".id", ".N.bin"),
  output:
    "data/grm/gcta.grm.sp",
  params:
    prefix = "data/grm/gcta",
  threads: 35
  benchmark: "benchmarks/gcta_sparsify_grm.txt"
  shell:
    """
    gcta \
      --grm {params.prefix} \
      --make-bK-sparse 0.05 \
      --thread-num {threads} \
      --out {params.prefix}
    """


rule gwas:
  input:
    geno = multiext("data/geno/geno_500k", ".pgen", ".psam", ".pvar"),
    pheno = "data/pheno/pheno_jan2024.tsv",
    covar = "data/pheno/covar.tsv",
  output:
    expand(
        "data/gwas/plink.{phenotype}.glm.linear.zst",
        phenotype=phenotype_names
    ),
  params:
    geno_prefix = "data/geno/geno_500k",
    output_prefix = "data/gwas/plink",
  threads: 35
  benchmark: "benchmarks/gwas.txt"
  shell:
    """
    plink2 \
      --pfile {params.geno_prefix} \
      --glm zs hide-covar \
      --pheno {input.pheno} \
      --covar {input.covar} \
      --threads {threads} \
      --out {params.output_prefix}
    """
