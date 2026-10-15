#!/usr/bin/env python3
"""Write the FragPipe 24.0 SDRF template (``sdrf.tsv``) for a run.

Faithful port of the FragPipe writer that runs after a FragPipe job when
``workflow.misc.save-sdrf=true``:

* FragpipeRun.java          (save trigger, TMT label, FASTA protein headers)
* tabs/TabWorkflow.java     sdrfSave(), manifestLoad(), on(MessageSDRFsave)
* util/SDRFtable.java       header layout, enzyme map, organism guess, Unimod match, printTable()
* tabs/TabMsfragger.java    getSDRFenzymes(), getSDRFmods(), getOffsetsAsMods(),
                            getPrecTolString(), getProdTolString(), getInstrumentMap()
* batmass-io UnimodOboReader (unimodDB entries "mass;id;name;sites")
* tools/tmtintegrator/QuantLabel.java (label types and reagent names)

The MSFragger UI state is rebuilt from the ``msfragger.*`` keys of FragPipe's
saved UI state (``fragpipe.workflow``, every tab key) exactly as FragPipe loads
it (ui.translations, UiCombo, document filters).  Java number formatting (HALF_UP ``%.Nf``,
``Double.toString``) and ``float`` arithmetic are emulated so the output is
byte-identical.
"""

import argparse
import os
import re
import struct
import sys
from decimal import ROUND_HALF_UP, Decimal
from fractions import Fraction
from math import isnan

__version__ = "1.0.0"

# --------------------------------------------------------------------------------------
# Java emulation helpers
# --------------------------------------------------------------------------------------


def _f32_bits(x):
    return struct.unpack("<I", struct.pack("<f", x))[0]


def _bits_f32(b):
    return struct.unpack("<f", struct.pack("<I", b))[0]


def java_parse_float(s):
    """Java Float.parseFloat: correctly rounded decimal -> float32 (no double rounding)."""
    s = s.strip()
    if s[-1:] in ("f", "F", "d", "D"):
        s = s[:-1]
    d = float(s)
    if d != d or d in (float("inf"), float("-inf")) or d == 0.0:
        return struct.unpack("<f", struct.pack("<f", d))[0]
    exact = Fraction(s)
    f = struct.unpack("<f", struct.pack("<f", d))[0]
    bits = _f32_bits(f)
    best = f
    best_err = abs(Fraction(f) - exact)
    for nb in (bits - 1, bits + 1):
        if nb < 0 or nb > 0xFFFFFFFF:
            continue
        g = _bits_f32(nb)
        if g != g:
            continue
        err = abs(Fraction(g) - exact)
        if err < best_err or (err == best_err and (nb & 1) == 0):
            best, best_err = g, err
    return best


def java_parse_int(s):
    """Java Integer.parseInt: no surrounding whitespace allowed."""
    if not re.fullmatch(r"[+-]?\d+", s):
        raise ValueError('For input string: "%s"' % s)
    return int(s)


def to_f32(x):
    """Java (float) cast of a double."""
    return struct.unpack("<f", struct.pack("<f", x))[0]


def java_fixed(x, digits):
    """Java String.format("%.<digits>f", double): HALF_UP on the shortest decimal repr."""
    if isnan(x):
        return "NaN"
    q = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP)
    return format(q, "f")


def java_double_to_string(x):
    """Java Double.toString(double)."""
    if isnan(x):
        return "NaN"
    if x in (float("inf"), float("-inf")):
        return "Infinity" if x > 0 else "-Infinity"
    if x == 0.0:
        return "-0.0" if str(x).startswith("-") else "0.0"
    sign = "-" if x < 0 else ""
    d = Decimal(repr(abs(x)))
    t = d.as_tuple()
    digits = "".join(str(i) for i in t.digits).rstrip("0") or "0"
    # exponent of the leading digit
    exp10 = len(t.digits) + t.exponent - 1
    if 1e-3 <= abs(x) < 1e7:
        if exp10 >= 0:
            int_part = digits[: exp10 + 1].ljust(exp10 + 1, "0")
            frac = digits[exp10 + 1 :] or "0"
        else:
            int_part = "0"
            frac = "0" * (-exp10 - 1) + digits
        return f"{sign}{int_part}.{frac}"
    mant = digits[0] + "." + (digits[1:] or "0")
    return f"{sign}{mant}E{exp10}"


def java_string_hash(s):
    h = 0
    for unit in struct.unpack("<%dH" % (len(s.encode("utf-16-le")) // 2), s.encode("utf-16-le")):
        h = (31 * h + unit) & 0xFFFFFFFF
    return h


def java_hashmap_order(keys_in_insertion_order):
    """Iteration order of a java.util.HashMap<String, ?> filled in this order (no treeified bins)."""
    n = len(keys_in_insertion_order)
    cap = 16
    while n > cap * 0.75:
        cap *= 2

    def bucket(k):
        h = java_string_hash(k)
        return (h ^ (h >> 16)) & (cap - 1)

    indexed = list(enumerate(keys_in_insertion_order))
    indexed.sort(key=lambda t: (bucket(t[1]), t[0]))
    return [k for _, k in indexed]


def java_split(s, regex):
    """Java String.split(regex): trailing empty strings removed."""
    parts = re.split(regex, s)
    while parts and parts[-1] == "":
        parts.pop()
    return parts if parts else ([s] if s == "" else [])


def java_file_name(path):
    """Paths.get(p).getFileName().toString(); manifests written on Windows use '\\'."""
    p = path.rstrip("/\\")
    return re.split(r"[/\\]", p)[-1]


# Java '.' does not match \n, \r, \u0085,  ,
_JDOT = "[^\n\r\u0085  ]"

# --------------------------------------------------------------------------------------
# Java .properties reader (java.util.Properties.load)
# --------------------------------------------------------------------------------------


def _logical_lines(text):
    lines = re.split(r"\r\n|\r|\n", text)
    i = 0
    while i < len(lines):
        line = lines[i].lstrip(" \t\f")
        i += 1
        if not line or line[0] in "#!":
            continue
        # an odd number of trailing backslashes continues the logical line
        while (len(line) - len(line.rstrip("\\"))) % 2 == 1:
            line = line[:-1]
            if i >= len(lines):
                break
            line += lines[i].lstrip(" \t\f")
            i += 1
        yield line


def _unescape(s):
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "\\" and i + 1 < len(s):
            n = s[i + 1]
            if n == "u":
                out.append(chr(int(s[i + 2 : i + 6], 16)))
                i += 6
                continue
            out.append({"t": "\t", "n": "\n", "r": "\r", "f": "\f"}.get(n, n))
            i += 2
            continue
        if c == "\\":
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def read_properties(path):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    props = {}
    for line in _logical_lines(text):
        key_end = value_start = len(line)
        sep_is_ws = False
        i = 0
        while i < len(line):
            c = line[i]
            if c == "\\":
                i += 2
                continue
            if c in "=:":
                key_end, value_start = i, i + 1
                break
            if c in " \t\f":
                key_end, value_start, sep_is_ws = i, i + 1, True
                break
            i += 1
        j = value_start
        while j < len(line) and line[j] in " \t\f":
            j += 1
        if sep_is_ws and j < len(line) and line[j] in "=:":
            j += 1
            while j < len(line) and line[j] in " \t\f":
                j += 1
        props[_unescape(line[:key_end])] = _unescape(line[j:])
    return props


# --------------------------------------------------------------------------------------
# FragPipe constants
# --------------------------------------------------------------------------------------

# SDRFtable.java column names
COL_source = "source name"
COL_organism = "characteristics[organism]"
COL_strain = "characteristics[strain/breed]"
COL_cultivar = "characteristics[ecotype/cultivar]"
COL_ancestry = "characteristics[ancestry category]"
COL_age = "characteristics[age]"
COL_stage = "characteristics[developmental stage]"
COL_sex = "characteristics[sex]"
COL_disease = "characteristics[disease]"
COL_part = "characteristics[organism part]"
COL_cellType = "characteristics[cell type]"
COL_indvidual = "characteristics[individual]"
COL_cellLine = "characteristics[cell line]"
COL_bioReplicate = "characteristics[biological replicate]"
COL_technology = "technology type"
COL_assay = "assay name"
COL_datafile = "comment[data file]"
COL_replicate = "comment[technical replicate]"
COL_fraction = "comment[fraction identifier]"
COL_label = "comment[label]"
COL_instrument = "comment[instrument]"
COL_enzyme = "comment[cleavage agent details]"
COL_mods = "comment[modification parameters]"
COL_precTol = "comment[precursor mass tolerance]"
COL_prodTol = "comment[fragment mass tolerance]"
EMPTY_STR = "not available"

SDRF_TYPES = ["Default", "Human", "CellLines", "Vertebrates", "NonVertebrates", "Plants"]

# SDRFtable.enzymesMap (keys "custom " and "nocleavage " carry a trailing space in FragPipe)
ENZYMES_MAP = {
    "custom ": "custom",
    "null": "no cleavage",
    "trypsin": "Trypsin",
    "stricttrypsin": "Trypsin/P",
    "argc": "Arg-C",
    "aspn": "Asp-N",
    "chymotrypsin": "Chymotrypsin",
    "cnbr": "CNBr",
    "elastase": "leukocyte elastase",
    "formicacid": "Formic acid",
    "gluc": "glutamyl endopeptidase",
    "gluc_bicarb": "glutamyl endopeptidase",
    "lysc": "Lys-C",
    "lysc-p": "Lys-C/P",
    "lysn": "Lys-N",
    "lysn_promisc": "Lys-N",
    "nonspecific": "unspecific cleavage",
    "trypsin/chymotrypsin": "TrypChymo",
    "trypsin/cnbr": "Trypsin_CNBr",
    "trypsin_gluc": "Trypsin_GluC",
    "trypsin_k": "Lys-C",
    "trypsin_r": "Arg-C",
    "thermolysin": "Thermolysin",
    "nocleavage ": "no cleavage",
}


def _reagents(s):
    return java_split(s, r"[,\s]+")


# QuantLabel.LABELS: name -> (type, reagent names)
QUANT_LABELS = [
    ("TMT", "TMT-0", ["126"]),
    ("TMT", "TMT-2", ["126", "127C"]),
    ("TMT", "TMT-6", _reagents("126, 127N, 128C, 129N, 130C, 131N")),
    ("TMT", "TMT-10", _reagents("126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N")),
    ("TMT", "TMT-11", _reagents("126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C")),
    (
        "TMT",
        "TMT-16",
        _reagents("126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N"),
    ),
    (
        "TMT",
        "TMT-18",
        _reagents(
            "126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N"
        ),
    ),
    (
        "TMT",
        "TMT-35",
        _reagents(
            "126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, "
            "135N, 127D, 128ND, 128CD, 129ND, 129CD, 130ND, 130CD, 131ND, 131CD, 132ND, 132CD, 133ND, 133CD, "
            "134ND, 134CD, 135ND, 135CD"
        ),
    ),
    ("TMT", "iodoTMT-6", _reagents("126, 127, 128, 129, 130, 131")),
    ("iTRAQ", "iTRAQ-4", _reagents("114, 115, 116, 117")),
    ("iTRAQ", "iTRAQ-8", _reagents("113, 114, 115, 116, 117, 118, 119, 121")),
    ("sCLIP", "sCLIP-6", _reagents("300, 301N, 301C, 302N, 302O, 302C")),
    (
        "IBT",
        "IBT-16",
        _reagents("114, 115N, 115C, 116N, 116C, 117N, 117C, 118N, 118C, 119N, 119C, 120N, 120C, 121N, 121C, 122"),
    ),
    ("DiLeu", "DiLeu-12", _reagents("115a, 115b, 116a, 116b, 116c, 117a, 117b, 117c, 118a, 118b, 118c, 118d")),
    ("DiLeu", "DiLeu-1", _reagents("114")),
    (
        "DeAla",
        "DeAla-13",
        ["100", "101a", "101b", "102a", "102b", "102c", "103a", "103b", "103c", "104a", "104b", "105a", "105b"],
    ),
    ("Trp", "Trp-2", _reagents("597, 603")),
    ("xtag", "develop", ["xtag%d" % i for i in range(1, 33)]),
    (
        "xtag2",
        "develop2",
        _reagents(
            "114, 115a, 115b, 115c, 116a, 116b, 116c, 116d, 116e, 117a, 117b, 117c, 117d, 117e, 117f, 118a, 118b, "
            "118c, 118d, 118e, 118f, 118g, 119a, 119b, 119c, 119d, 119e, 119f, 119g"
        ),
    ),
]


def get_quant_label(name):
    """TmtiPanel.getSelectedLabel(): case-insensitive match on QuantLabel name."""
    for qtype, qname, reagents in QUANT_LABELS:
        if qname.lower() == name.lower():
            return qtype, reagents
    raise ValueError("Unknown label name in combo box: " + name)


# ui.translations entries (UI value <=> config value) for the msfragger keys used here
UI_TRANSLATIONS = {
    "msfragger.calibrate_mass": [
        ("None", "0"),
        ("Off", "0"),
        ("Mass calibration", "1"),
        ("Mass calibration, parameter optimization", "2"),
    ],
    "msfragger.fragment_mass_units": [("PPM", "1"), ("Da", "0")],
    "msfragger.precursor_mass_units": [("PPM", "1"), ("Da", "0")],
    "msfragger.labile_search_mode": [("off", "off"), ("off", "standard"), ("nglycan", "nglycan"), ("labile", "labile")],
}

CALIBRATE_LABELS = ["None", "Mass calibration", "Mass calibration, parameter optimization"]
UNITS_MODEL = ["PPM", "Da"]  # PrecursorMassTolUnits / MassTolUnits enum order
SENSE_MODEL = ["N", "C"]
GLYCO_MODEL = ["off", "nglycan", "labile"]
GLYCO_OPTION_nglycan = "nglycan"


def translate_to_ui(key, value):
    """FragpipeCacheUtils.translateValuesToUi for a single entry."""
    ts = UI_TRANSLATIONS.get(key)
    if ts is None:
        return value
    if any(ui == value for ui, _ in ts):
        return value
    for ui, conf in ts:
        if conf.lower() == value.lower():
            return ui
    return value


def combo_from_string(model, value):
    """UiCombo.fromString: unknown values select index 0."""
    return value if value in model else model[0]


# --------------------------------------------------------------------------------------
# MSFragger UI state from the workflow
# --------------------------------------------------------------------------------------


class Mod:
    __slots__ = ("mass_delta", "sites", "is_enabled", "max_occurrences")

    def __init__(self, mass_delta, sites, is_enabled, max_occurrences):
        self.mass_delta = mass_delta  # Java float (stored as python float holding a float32 value)
        self.sites = sites
        self.is_enabled = is_enabled
        self.max_occurrences = max_occurrences


def parse_mods_table(s):
    """ModsTable.fromString + Mod.fromString + ModsTableModel.getModifications."""
    mods = []
    for chunk in s.split("; "):
        split = java_split(chunk, ",")
        if len(split) != 4:
            continue  # ValidationException -> logged and skipped by FragPipe
        mod = Mod(java_parse_float(split[0]), split[1], split[2].lower() == "true", java_parse_int(split[3]))
        if mod.sites.strip() == "":
            continue
        mods.append(mod)
    return mods


class MsfraggerUi:
    # props: FragPipe's saved UI state (fragpipe/fragpipe.workflow, FRAGPIPESAVEWORKFLOW), which holds every
    # MSFragger tab key: FragPipe writes sdrf.tsv from that tab state (FragpipeRun.java:583-587, TabWorkflow.java:1410)
    def __init__(self, props):

        def ui(key):
            return translate_to_ui(key, props[key])

        def cuts(s):
            return re.sub(r"[^A-Z@]", "@", s)

        self.enzyme_name_1 = props["msfragger.search_enzyme_name_1"]
        self.cuts_1 = cuts(props["msfragger.search_enzyme_cut_1"])
        self.nocuts_1 = cuts(props["msfragger.search_enzyme_nocut_1"])
        self.sense_1 = combo_from_string(SENSE_MODEL, props["msfragger.search_enzyme_sense_1"])
        self.enzyme_name_2 = props["msfragger.search_enzyme_name_2"]
        self.cuts_2 = cuts(props["msfragger.search_enzyme_cut_2"])
        self.nocuts_2 = cuts(props["msfragger.search_enzyme_nocut_2"])
        self.sense_2 = combo_from_string(SENSE_MODEL, props["msfragger.search_enzyme_sense_2"])
        self.var_mods = parse_mods_table(props["msfragger.table.var-mods"])
        self.fix_mods = parse_mods_table(props["msfragger.table.fix-mods"])
        self.mass_offsets = props["msfragger.mass_offsets"]
        self.mass_offsets_detailed = props["msfragger.mass_offsets_detailed"]
        self.use_detailed_offsets = props["msfragger.use_detailed_offsets"].lower() == "true"
        self.restrict_deltamass_to = re.sub(r"[^A-Zal\-nc\[\]^*]", "", props["msfragger.restrict_deltamass_to"])
        self.glyco = combo_from_string(GLYCO_MODEL, ui("msfragger.labile_search_mode"))
        self.prec_units = combo_from_string(UNITS_MODEL, ui("msfragger.precursor_mass_units"))
        self.prec_lo = float(props["msfragger.precursor_mass_lower"])
        self.prec_hi = float(props["msfragger.precursor_mass_upper"])
        self.frag_units = combo_from_string(UNITS_MODEL, ui("msfragger.fragment_mass_units"))
        self.frag_tol = float(props["msfragger.fragment_mass_tolerance"])
        self.calibrate = combo_from_string(CALIBRATE_LABELS, ui("msfragger.calibrate_mass"))


# --------------------------------------------------------------------------------------
# Unimod (batmass-io UnimodOboReader) and SDRFtable.matchUnimod
# --------------------------------------------------------------------------------------

_P_ID = re.compile(r"id:\s+UNIMOD:(\d+)")
_P_MASS = re.compile(r'xref:\s+delta_mono_mass\s+"([\d.+-]+)"')
_P_SITE = re.compile(r'xref:\s+spec_\d+_site\s+"(\S+)"')
_ID_UNSET = -2147483648


def read_unimod(path):
    """UnimodOboReader.unimodDB as (mass_float32, name, sites_string) tuples."""
    db = []
    term = None

    def flush(t):
        db.append((t["mass"], t["name"], "".join(t["sites"])))

    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("[Term]"):
                if term is not None:
                    flush(term)
                term = {"id": _ID_UNSET, "mass": float("nan"), "sites": [], "name": ""}
                continue
            if line.startswith("id:"):
                m = _P_ID.fullmatch(line)
                if not m:
                    raise RuntimeError("Could not parse unimod id from line: " + line)
                term["id"] = int(m.group(1))
                continue
            if term is None or term["id"] == _ID_UNSET:
                continue
            if line.startswith("xref: delta_mono_mass"):
                m = _P_MASS.fullmatch(line)
                if not m:
                    raise RuntimeError("Could not parse unimod mass from line: " + line)
                term["mass"] = java_parse_float(m.group(1))
                continue
            if line.startswith("name:"):
                term["name"] = line.split(": ")[1]
                continue
            m = _P_SITE.fullmatch(line)
            if m:
                site = m.group(1)
                ch = "n" if site == "N-term" else "c" if site == "C-term" else site[0]
                if ch not in term["sites"]:
                    term["sites"].append(ch)
    if term is not None:
        flush(term)
    return db


def match_unimod(unimod_db, mass, sites_with_termini):
    """SDRFtable.matchUnimod: the single match, or None for 0 or >1 matches."""
    matches = []
    for unimod_mass, name, unimod_sites in unimod_db:
        if abs(mass - unimod_mass) < 0.0001:
            if unimod_sites == "":
                # Java: "mass;id;name;".split(";") has 3 elements -> splits[3] throws
                raise IndexError("Unimod entry '%s' has no sites (FragPipe throws here)" % name)
            all_matched = True
            for site in sites_with_termini:
                if site == "N-term":
                    site = "n"
                elif site == "C-term":
                    site = "c"
                if site not in unimod_sites:
                    all_matched = False
                    break
            if all_matched:
                matches.append(name)
    return matches[0] if len(matches) == 1 else None


# --------------------------------------------------------------------------------------
# TabMsfragger SDRF helpers
# --------------------------------------------------------------------------------------


def map_enzyme_to_sdrf(name):
    return ENZYMES_MAP.get(name, name)


def get_sdrf_enzymes(ui):
    enzymes = []
    for name, cuts, nocuts, sense in (
        (ui.enzyme_name_1, ui.cuts_1, ui.nocuts_1, ui.sense_1),
        (ui.enzyme_name_2, ui.cuts_2, ui.nocuts_2, ui.sense_2),
    ):
        if cuts:
            regex = "(?%s[%s])" % ("<=" if sense == "C" else "=", cuts)
            if nocuts:
                regex += "(?![%s])" % nocuts
            enzymes.append("NT=%s;CS=%s" % (map_enzyme_to_sdrf(name), regex))
    return enzymes


def _varmod_term(sites):
    if "[" in sites:
        return "Protein N-term"
    if "]" in sites:
        return "Protein C-term"
    if "n" in sites:
        return "Any N-term"
    if "c" in sites:
        return "Any C-term"
    return "Anywhere"


def _fixedmod_term(sites):
    low = sites.lower()
    if "n-term" in low:
        return "Protein N-term" if "protein" in low else "Any N-term"
    if "c-term" in low:
        return "Protein C-term" if "protein" in low else "Any C-term"
    return "Anywhere"


def _update_sites_term(terminal, sites):
    out = list(sites)
    if terminal in ("Protein N-term", "Any N-term"):
        out.append("N-term")
    elif terminal in ("Protein C-term", "Any C-term"):
        out.append("C-term")
    return out


def _mod_name(mass, sites, unimod_name):
    return "%s:%s" % (sites, java_fixed(mass, 0)) if unimod_name is None else unimod_name


_MASS_OFFSET_SITES = re.compile(r"([A-Zalnc\[\]^*-])")
_MASS_OFFSET_MASS = re.compile(r"([\d.-]+)\(")
_MASS_OFFSET_AA = re.compile(r"aa=([A-Z-nc\[\]^*]+)")


def get_offsets_as_mods(ui):
    """TabMsfragger.getOffsetsAsMods (+ getMassOffsetSet / getDetailedOffsets)."""
    out = []
    if ui.use_detailed_offsets:
        if ui.mass_offsets_detailed:
            for split in java_split(ui.mass_offsets_detailed, ";"):
                m = _MASS_OFFSET_MASS.search(split)
                mass = java_parse_float(m.group(1)) if m else 0.0
                a = _MASS_OFFSET_AA.search(split)
                residues = _MASS_OFFSET_SITES.findall(a.group(1)) if a else []
                out.append(Mod(mass, "".join(residues), True, 1))
    else:
        values = set()
        for s in java_split(ui.mass_offsets.strip(), r"[/\s]+"):
            if s == "":
                continue
            v = java_parse_float(s)
            if abs(v) > to_f32(0.01):
                values.add(v)
        residues = ui.restrict_deltamass_to.replace("-", "")
        if residues == "all":
            residues = "ACDEFGHIKLMNPQRSTVWY"
        for v in sorted(values):
            out.append(Mod(v, residues, True, 1))
    return out


def get_sdrf_mods(ui, unimod_db):
    mods = []
    for mod in ui.var_mods:
        if mod.is_enabled:
            sites = re.findall(r"[A-Z]", mod.sites)
            terminal = _varmod_term(mod.sites)
            sites_with_termini = _update_sites_term(terminal, sites)
            sites_temp = "".join(sites)
            unimod_name = match_unimod(unimod_db, mod.mass_delta, sites_with_termini)
            name = _mod_name(mod.mass_delta, mod.sites, unimod_name)
            site_str = "" if sites_temp == "" else ";TA=%s" % sites_temp
            mods.append(
                "NT=%s;MT=Variable;PP=%s%s;MM=%s" % (name, terminal, site_str, java_fixed(mod.mass_delta, 5))
            )
    for mod in ui.fix_mods:
        if mod.is_enabled and mod.mass_delta != 0.0:
            terminal = _fixedmod_term(mod.sites)
            sites = []
            if terminal == "Anywhere":
                sites = re.findall(r"[A-Z]", mod.sites)
                sites_str = ";TA=%s" % "".join(sites)
            else:
                sites_str = ""
            sites_with_termini = _update_sites_term(terminal, sites)
            unimod_name = match_unimod(unimod_db, mod.mass_delta, sites_with_termini)
            name = _mod_name(mod.mass_delta, mod.sites, unimod_name)
            mods.append("NT=%s;MT=Fixed;PP=%s%s;MM=%s" % (name, terminal, sites_str, java_fixed(mod.mass_delta, 5)))
    for offset in get_offsets_as_mods(ui):
        if offset.mass_delta != 0.0:
            sites = list(offset.sites) if offset.sites else [""]  # Java "".split("") == [""]
            unimod_name = match_unimod(unimod_db, offset.mass_delta, sites)
            name = unimod_name if unimod_name is not None else "offset:%s" % java_fixed(offset.mass_delta, 5)
            glyco_site = ""
            offset_residues = offset.sites
            if ui.glyco == GLYCO_OPTION_nglycan:
                glyco_site = "TS=N[^P][ST];"
                offset_residues = "N"
            mods.append(
                "NT=%s;MT=Variable;PP=%s;TA=%s;%sMM=%s"
                % (name, "Anywhere", offset_residues, glyco_site, java_fixed(offset.mass_delta, 5))
            )
    return mods


def get_prec_tol_string(ui):
    tolerance = ui.prec_hi - ui.prec_lo
    if ui.prec_units == "PPM":
        return "%s %s" % (java_fixed(tolerance, 1), ui.prec_units.lower())
    return "%s %s" % (java_fixed(tolerance, 2), ui.prec_units)


_NEW_FRAG_TOL = re.compile(r"New fragment_mass_tolerance = ([0-9\.]+) ([\w]+)", re.ASCII)


def get_prod_tol_string(ui, log_text):
    use_log = False
    if ui.calibrate == "Mass calibration, parameter optimization":
        m = _NEW_FRAG_TOL.search(log_text)
        if m:
            tolerance = float(m.group(1))
            unit = m.group(2)
            use_log = True
    if not use_log:
        unit = ui.frag_units
        tolerance = ui.frag_tol
    if unit == "PPM":
        return "%s %s" % (java_double_to_string(tolerance), unit.lower())
    return "%s %s" % (java_fixed(tolerance, 2), unit)


_INSTRUMENT = re.compile("(%s+): Scans %s+ Instrument = ([^;\\n\\r]+)" % (_JDOT, _JDOT))


def get_instrument_map(log_text):
    out = {}
    for m in _INSTRUMENT.finditer(log_text):
        out[java_file_name(m.group(1)).strip()] = m.group(2).strip()
    return out


# --------------------------------------------------------------------------------------
# SDRFtable
# --------------------------------------------------------------------------------------

_UNIPROT_OS = re.compile(r"OS=([\w\s]+) (OX=|GN=|PE=|SV=)", re.ASCII)
DATABASE_SIZE_LIMIT = 1 << 30  # TabDatabase.databaseSizeLimit


def read_protein_headers(fasta):
    """FragpipeRun: sorted '>' lines of the FASTA if it is smaller than 1 GiB."""
    if os.path.getsize(fasta) >= DATABASE_SIZE_LIMIT:
        return []
    try:
        with open(fasta, encoding="utf-8", errors="strict", newline=None) as fh:
            headers = [line.rstrip("\n") for line in fh if line.startswith(">")]
    except UnicodeDecodeError:
        # On purpose, as FragPipe 24.0: Files.lines(p) (UTF-8) throws on a malformed byte, the catch (Exception e)
        # only logs "Cannot get the sample of ...", and proteinHeaders stays the empty list it was initialised to
        # (FragpipeRun.java:436-454).
        return []
    return sorted(headers)


def guess_organism(protein_headers):
    counts = {}
    order = []
    for header in protein_headers:
        m = _UNIPROT_OS.search(header)
        if m:
            org = m.group(1)
            if org not in counts:
                counts[org] = 0
                order.append(org)
            counts[org] += 1
    highest = 0
    highest_org = EMPTY_STR
    for org in java_hashmap_order(order):
        if counts[org] > highest:
            highest_org = org
            highest = counts[org]
    return highest_org


def header_characteristics(sdrf_type):
    """SDRFtable.initHeaderCharacteristics (NonVertebrates/Plants fall through: no break)."""
    h = [COL_source, COL_organism]
    if sdrf_type == "Human":
        h += [COL_ancestry, COL_age, COL_stage, COL_sex, COL_indvidual]
    elif sdrf_type == "Vertebrates":
        h += [COL_age, COL_stage, COL_sex, COL_indvidual]
    elif sdrf_type in ("NonVertebrates", "Plants", "CellLines"):
        if sdrf_type == "NonVertebrates":
            h += [COL_strain, COL_age, COL_stage, COL_sex, COL_indvidual]
        if sdrf_type in ("NonVertebrates", "Plants"):
            h += [COL_cultivar, COL_age, COL_stage, COL_indvidual]
        h += [COL_cellLine]
    h += [COL_disease, COL_part, COL_cellType, COL_bioReplicate, COL_technology, COL_assay]
    return h


class SdrfTable:
    def __init__(self, sdrf_type, num_enzymes, num_mods, protein_headers):
        self.rows = []
        self.guessed_organism = guess_organism(protein_headers)
        self.header = header_characteristics(sdrf_type)
        self.header += [COL_datafile, COL_replicate, COL_fraction, COL_label, COL_instrument, COL_precTol, COL_prodTol]
        self.first_enzyme_index = len(self.header)
        self.header += [COL_enzyme] * num_enzymes
        self.first_mod_index = len(self.header)
        self.header += [COL_mods] * num_mods

    def _row(self, file_name, replicate, enzymes, mods, prec_tol, prod_tol, instrument):
        row = [None] * len(self.header)
        idx = self.header.index
        row[idx(COL_organism)] = self.guessed_organism
        row[idx(COL_datafile)] = file_name
        row[idx(COL_replicate)] = replicate
        row[idx(COL_precTol)] = prec_tol
        row[idx(COL_prodTol)] = prod_tol
        row[idx(COL_instrument)] = EMPTY_STR if instrument is None else "NT=" + instrument
        for i, e in enumerate(enzymes):
            row[self.first_enzyme_index + i] = e
        for i, m in enumerate(mods):
            row[self.first_mod_index + i] = m
        return row

    def add_sample_lfq(self, file_name, replicate, enzymes, mods, prec_tol, prod_tol, instrument):
        self.rows.append(self._row(file_name, replicate, enzymes, mods, prec_tol, prod_tol, instrument))

    def add_sample_tmt(self, file_name, replicate, enzymes, mods, quant_label, prec_tol, prod_tol, instrument):
        qtype, reagents = quant_label
        for label in reagents:
            row = self._row(file_name, replicate, enzymes, mods, prec_tol, prod_tol, instrument)
            row[self.header.index(COL_label)] = "%s%s" % (qtype, label)
            self.rows.append(row)

    def render(self):
        lines = ["\t".join(self.header)]
        for row in self.rows:
            lines.append("\t".join(EMPTY_STR if v is None else v for v in row))
        return "\n".join(lines)


# --------------------------------------------------------------------------------------
# Manifest (TabWorkflow.manifestLoad) and main
# --------------------------------------------------------------------------------------


def read_manifest(path):
    """Returns [(file_name, replicate_string)] in manifest order."""
    files = []
    with open(path, encoding="utf-8") as fh:
        for line in re.split(r"\r\n|\r|\n", fh.read()):
            if line.strip() == "" or line.startswith("//") or line.startswith("#"):
                continue
            split = java_split(line.strip(), "\t")
            replicate = ""
            if len(split) >= 3 and split[2].strip() != "":
                replicate = str(java_parse_int(split[2]))
            files.append((java_file_name(split[0]), replicate))
    return files


def read_searched_names(path):
    """Manifest file name -> name of the file MSFragger searched (tab-separated, one per line)."""
    names = {}
    if path is None:
        return names
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            manifest_name, searched_name = line.split("\t")
            names[manifest_name] = searched_name
    return names


def build_sdrf(workflow, manifest, fasta, unimod, tmt_label, logs, searched_names=None):
    props = read_properties(workflow)
    sdrf_type = combo_from_string(SDRF_TYPES, props.get("workflow.misc.sdrf-type", SDRF_TYPES[0]))
    ui = MsfraggerUi(props)
    unimod_db = read_unimod(unimod)
    log_text = "\n".join(open(p, encoding="utf-8").read() for p in logs)
    label = get_quant_label(tmt_label) if tmt_label else None

    enzymes = get_sdrf_enzymes(ui)
    mods = get_sdrf_mods(ui, unimod_db)
    prec_tol = get_prec_tol_string(ui)
    prod_tol = get_prod_tol_string(ui, log_text)
    instrument_map = get_instrument_map(log_text)
    searched = read_searched_names(searched_names)

    table = SdrfTable(sdrf_type, len(enzymes), len(mods), read_protein_headers(fasta))
    for file_name, replicate in read_manifest(manifest):
        if label is not None:
            table.add_sample_tmt(
                file_name, replicate, enzymes, mods, label, prec_tol, prod_tol, instrument_map.get(searched.get(file_name, file_name))
            )
        else:
            table.add_sample_lfq(
                file_name, replicate, enzymes, mods, prec_tol, prod_tol, instrument_map.get(searched.get(file_name, file_name).strip())
            )
    return table.render()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Write the FragPipe SDRF template (sdrf.tsv).")
    ap.add_argument("--workflow", required=True, help="FragPipe's saved UI state (fragpipe.workflow) of the run")
    ap.add_argument("--manifest", required=True, help="fragpipe-files.fp-manifest")
    ap.add_argument("--searched-names", default=None, help="TSV: manifest file name, name of the file MSFragger searched (instrument lookup, TabMsfragger.java:2388)")
    ap.add_argument("--fasta", required=True, help="protein database FASTA")
    ap.add_argument("--unimod", required=True, help="unimod.obo shipped in FragPipe tools/")
    ap.add_argument("--tmt-label", default="", help="TMT-Integrator QuantLabel name (e.g. TMT-16); empty for LFQ")
    ap.add_argument("--log", action="append", default=[], help="MSFragger log (repeatable, concatenated in order)")
    ap.add_argument("--output", required=True, help="output sdrf.tsv")
    ap.add_argument("--version", action="version", version="%(prog)s " + __version__)
    a = ap.parse_args(argv)
    text = build_sdrf(a.workflow, a.manifest, a.fasta, a.unimod, a.tmt_label.strip(), a.log, a.searched_names)
    with open(a.output, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    print("Saved SDRF template to: " + a.output, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
