def format_wildcard_metrics(
    wc_name: str,
    total_reads: int,
    adapter_reads: int,
    unique_count: int,
    num_to_ninety: int,
    expected_count: int,
    whitelist_count: int,
    total_read_from_expected: float,
    pct_for_50: str,
    pct_for_96: str,
) -> str:
    return (
        f"\n######### Info for {wc_name}\n"
        f"Total Reads: {total_reads}  Reads with Adapter: {adapter_reads}\n"
        f"The number of unique strings in the 8mer column is {unique_count}.\n"
        f"Ninety percent (90%) of the reads come from total of {num_to_ninety} 8mers.\n"
        f"Total of {expected_count} out of {whitelist_count} expected 8-mers "
        f"accounted for {total_read_from_expected:.1%} of the reads\n"
        f"Top 50 8mers represent {pct_for_50} fraction of reads\n"
        f"Top 96 8mers represent {pct_for_96} fraction of reads"
    )


def format_hilo_metrics(
    label: str,
    total_hi_warn: int,
    total_lo_warn: int,
    total_mers: int,
) -> str:
    out = (
        f"\n######### Info for {label}\n"
        f"Total hiWarn: {total_hi_warn}\tTotal loWarn: {total_lo_warn}"
    )
    if total_mers > 0:
        out += (
            f"\nPct hiWarn: {(total_hi_warn / total_mers):.3f}\t"
            f"Pct loWarn: {(total_lo_warn / total_mers):.3f}"
        )
    else:
        out += "\nPct hiWarn: N/A\tPct loWarn: N/A"
    return out
