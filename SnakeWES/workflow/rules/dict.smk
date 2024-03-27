rule create_dict:
    input:
        "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        "GRCh38_full_analysis_set_plus_decoy_hla.dict"
    log:
        "logs/create_dict.log",
    params:
        extra="",  # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v3.5.3/bio/picard/createsequencedictionary"