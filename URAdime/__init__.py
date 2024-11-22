"""URAdime - A package for analyzing primers in BAM files."""

from .URAdime import (
    load_primers,
    is_match,
    find_primers_in_region,
    process_read_chunk,
    create_analysis_summary,
    downsample_reads,
    bam_to_fasta_parallel,
    parallel_analysis_pipeline,
    parse_arguments,
    validate_inputs,
    create_primer_statistics,
    save_results,
    main
)

__version__ = "0.1.6"
