MISSING = {None, "", ".", "nan", "None"}


def parse_format(format_value, sample_value):
    """Map FORMAT keys to the corresponding values for one sample."""
    if format_value in MISSING or sample_value in MISSING:
        return {}
    return dict(zip(str(format_value).split(":"), str(sample_value).split(":")))


def require_format_fields(format_value, required_fields, caller):
    """Fail clearly if an upstream caller changes its expected FORMAT schema."""
    format_fields = set() if format_value in MISSING else set(str(format_value).split(":"))
    missing_fields = sorted(set(required_fields) - format_fields)
    if missing_fields:
        missing = ", ".join(missing_fields)
        raise ValueError(f"{caller} FORMAT is missing required field(s): {missing}")


def get_format_value(format_value, sample_value, key, default="."):
    value = parse_format(format_value, sample_value).get(key, default)
    return default if value in MISSING else value


def get_alt_depth(format_value, sample_value):
    """Return Sniffles2 variant-supporting read depth (DV)."""
    return get_format_value(format_value, sample_value, "DV", default="0")


def get_depth(format_value, sample_value):
    """Calculate Sniffles2 total depth as reference reads (DR) plus variant reads (DV)."""
    values = parse_format(format_value, sample_value)
    reference_depth = values.get("DR")
    variant_depth = values.get("DV")
    if reference_depth not in MISSING and variant_depth not in MISSING:
        try:
            return str(int(reference_depth) + int(variant_depth))
        except (TypeError, ValueError):
            pass
    return "."


def convert_ont_sample_value(format_value, sample_value, variant_type):
    """Convert one ONT sample field to the correct report layout."""
    genotype = get_format_value(format_value, sample_value, "GT")
    if variant_type == "CNV":
        # Spectre GT:HO:GQ:CN to the GT:CN layout expected by get_CN().
        require_format_fields(format_value, {"GT", "CN"}, "Spectre CNV")
        copy_number = get_format_value(format_value, sample_value, "CN")
        return f"{genotype}:{copy_number}"

    # Sniffles GT:GQ:DR:DV:PS:ID to GT:AD:DP:PS. The report only reads the
    # alternate part of AD, so 0 is used as a reference-depth placeholder.
    require_format_fields(format_value, {"GT", "DR", "DV"}, "Sniffles2 SV")
    alternate_depth = get_alt_depth(format_value, sample_value)
    depth = get_depth(format_value, sample_value)
    phase_set = get_format_value(format_value, sample_value, "PS")
    return f"{genotype}:0,{alternate_depth}:{depth}:{phase_set}"
