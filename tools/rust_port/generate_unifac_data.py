"""Generate rust/src/unifac_data.rs from src/unifacdata.f90.

The UNIFAC group database is a pure data source in the original fork:
117 sub-groups with their main-group assignment, Rk and Qk parameters,
20 component-to-group mappings and 1603 main-group interaction records.
The generated module keeps the source order so that first-match lookups
behave exactly like the Fortran `getUNIFACujkdbidx` scan.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "src/unifacdata.f90"

text = source.read_text(encoding="utf-8")


def strip_comments(text):
    """Remove Fortran comments while keeping string literals intact."""
    out = []
    for line in text.splitlines():
        cleaned = []
        in_string = False
        index = 0
        while index < len(line):
            char = line[index]
            if char == '"':
                in_string = not in_string
                cleaned.append(char)
            elif char == "!" and not in_string:
                break
            else:
                cleaned.append(char)
            index += 1
        out.append("".join(cleaned))
    return "\n".join(out)


clean = strip_comments(text)


def mask_strings(text):
    """Blank out string literal contents, keeping the original length.

    Some source formula strings contain the literal text of a constructor
    name, so constructor calls are located in the masked copy while the
    argument text is taken from the original.
    """
    masked = list(text)
    in_string = False
    for index, char in enumerate(text):
        if char == '"':
            in_string = not in_string
            masked[index] = " "
        elif in_string:
            masked[index] = " "
    return "".join(masked)


masked = mask_strings(clean)


def calls(name):
    """Yield the argument text of every `name(...)` constructor call."""
    for match in re.finditer(r"\b" + name + r"\s*\(", masked):
        start = match.end()
        depth = 1
        index = start
        while index < len(masked) and depth:
            char = masked[index]
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
            index += 1
        if depth:
            raise ValueError(f"unbalanced {name} call")
        yield clean[start : index - 1]


def split_arguments(body):
    """Split a Fortran argument list on top-level commas."""
    parts = []
    current = []
    depth = 0
    in_string = False
    for char in body:
        if char == '"':
            in_string = not in_string
        if not in_string:
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
            elif char == "," and depth == 0:
                parts.append("".join(current).strip())
                current = []
                continue
        current.append(char)
    parts.append("".join(current).strip())
    return parts


def number(text):
    """Parse a Fortran real literal and format it so that it round-trips.

    The guard rejects a value that the fixed format would not reproduce, so
    that a future data change cannot silently lose significant digits.
    """
    value = float(text.replace("d", "e").replace("D", "e"))
    formatted = f"{value:.10e}"
    if float(formatted) != value:
        raise SystemExit(f"value {text} does not round-trip through {formatted}")
    return value, formatted


def number_value(text):
    return number(text)[0]


subgroups = []
for body in calls("unifacprm"):
    args = split_arguments(body)
    subgroups.append(
        {
            "sub_group": int(args[0]),
            "main_group": int(args[1]),
            "formula": args[2].strip().strip('"'),
            "rk": number_value(args[3]),
            "qk": number_value(args[4]),
            "data_source": args[5].strip().strip('"'),
        }
    )

components = []
for body in calls("unifacComp"):
    args = split_arguments(body)
    groups = re.findall(r"-?\d+", args[2].split("(/")[-1].split("/)")[0])
    components.append(
        {
            "formula": args[0].strip().strip('"'),
            "uid": args[1].strip().strip('"'),
            "v": [int(g) for g in groups],
        }
    )

interactions = []
for body in calls("unifacUij"):
    args = split_arguments(body)
    interactions.append(
        {
            "main_group_i": int(args[0]),
            "main_group_j": int(args[1]),
            "aij": number_value(args[2]),
            "bij": number_value(args[3]),
            "cij": number_value(args[4]),
            "data_source": args[5].strip().strip('"'),
        }
    )

n_sub_groups = int(re.search(r"nSubGroups\s*=\s*(\d+)", clean).group(1))
n_components = int(re.search(r"nUnifacComp\s*=\s*(\d+)", clean).group(1))
n_interactions = int(re.search(r"nUnifacUij\s*=\s*(\d+)", clean).group(1))

if len(subgroups) != n_sub_groups:
    raise SystemExit(f"expected {n_sub_groups} sub-groups, parsed {len(subgroups)}")
if len(components) != n_components:
    raise SystemExit(f"expected {n_components} components, parsed {len(components)}")
if len(interactions) != n_interactions:
    raise SystemExit(f"expected {n_interactions} interactions, parsed {len(interactions)}")

for component in components:
    if len(component["v"]) != n_sub_groups:
        raise SystemExit(
            f"component {component['uid']} has {len(component['v'])} group counts"
        )

# The Fortran lookup returns the first matching record, so duplicated
# (main group i, main group j) pairs must keep the earliest entry.
seen = {}
duplicates = []
for position, record in enumerate(interactions):
    key = (record["main_group_i"], record["main_group_j"])
    if key in seen:
        duplicates.append((key, seen[key], position))
    else:
        seen[key] = position

# Report main-group pairs that the source database cannot resolve.
missing = []
for group_i in sorted({s["main_group"] for s in subgroups}):
    for group_j in sorted({s["main_group"] for s in subgroups}):
        if (group_i, group_j) not in seen:
            missing.append((group_i, group_j))

print(f"Parsed {len(subgroups)} sub-groups, {len(components)} components, "
      f"{len(interactions)} interaction records.")
print(f"Duplicate interaction keys: {len(duplicates)}")
print(f"Unresolvable main-group pairs: {len(missing)}")

out = []
out.append("//! Source-derived UNIFAC group, component and interaction database.")
out.append("//! Automatically generated by tools/rust_port/generate_unifac_data.py from src/unifacdata.f90.")
out.append("//! Copyright thermotools and Contributors. Apache-2.0.\n")
out.append("#![allow(clippy::excessive_precision)]\n")

out.append("/// Number of UNIFAC sub-groups in the source database.")
out.append(f"pub const N_SUB_GROUPS: usize = {n_sub_groups};")
out.append("/// Number of UNIFAC component records in the source database.")
out.append(f"pub const N_COMPONENTS: usize = {n_components};")
out.append("/// Number of main-group interaction records in the source database.")
out.append(f"pub const N_INTERACTIONS: usize = {n_interactions};\n")

out.append("/// A UNIFAC sub-group with its main-group assignment and Rk/Qk parameters.")
out.append("#[derive(Debug, Clone, Copy)]")
out.append("pub struct UnifacSubGroup {")
out.append("    /// One-based sub-group index in the source database.")
out.append("    pub sub_group: u8,")
out.append("    /// One-based main-group index used by the interaction parameters.")
out.append("    pub main_group: u8,")
out.append("    /// Group formula, for example CH3 or CH2=CH.")
out.append("    pub formula: &'static str,")
out.append("    /// Volume parameter Rk in [-].")
out.append("    pub rk: f64,")
out.append("    /// Surface area parameter Qk in [-].")
out.append("    pub qk: f64,")
out.append("    /// Source reference string of the record.")
out.append("    pub data_source: &'static str,")
out.append("}\n")

out.append("/// A component with its number of sub-groups of each kind.")
out.append("#[derive(Debug, Clone, Copy)]")
out.append("pub struct UnifacComponent {")
out.append("    /// Molecular formula recorded in the source.")
out.append("    pub formula: &'static str,")
out.append("    /// Component identifier used by the component list.")
out.append("    pub uid: &'static str,")
out.append("    /// Number of sub-groups of each of the N_SUB_GROUPS kinds.")
out.append(f"    pub v: [u8; N_SUB_GROUPS],")
out.append("}\n")

out.append("/// A main-group interaction record with its three-parameter energy expression.")
out.append("#[derive(Debug, Clone, Copy)]")
out.append("pub struct UnifacInteraction {")
out.append("    /// One-based first main-group index.")
out.append("    pub main_group_i: u8,")
out.append("    /// One-based second main-group index.")
out.append("    pub main_group_j: u8,")
out.append("    /// Attraction parameter aij in Kelvin.")
out.append("    pub aij: f64,")
out.append("    /// Interaction parameter bij in [-].")
out.append("    pub bij: f64,")
out.append("    /// Interaction parameter cij in [1/K].")
out.append("    pub cij: f64,")
out.append("    /// Source reference string of the record.")
out.append("    pub data_source: &'static str,")
out.append("}\n")

out.append(f"/// All {n_sub_groups} source-derived UNIFAC sub-group records.")
out.append(f"pub static UNIFAC_SUBGROUPS: [UnifacSubGroup; {n_sub_groups}] = [")
for record in subgroups:
    out.append("    UnifacSubGroup {")
    out.append(f'        sub_group: {record["sub_group"]},')
    out.append(f'        main_group: {record["main_group"]},')
    out.append(f'        formula: "{record["formula"]}",')
    out.append(f'        rk: {record["rk"]:.10e},')
    out.append(f'        qk: {record["qk"]:.10e},')
    out.append(f'        data_source: "{record["data_source"]}",')
    out.append("    },")
out.append("];\n")

out.append(f"/// All {n_components} source-derived UNIFAC component records.")
out.append(f"pub static UNIFAC_COMPONENTS: [UnifacComponent; {n_components}] = [")
for record in components:
    groups = ", ".join(str(g) for g in record["v"])
    out.append("    UnifacComponent {")
    out.append(f'        formula: "{record["formula"]}",')
    out.append(f'        uid: "{record["uid"]}",')
    out.append(f"        v: [{groups}],")
    out.append("    },")
out.append("];\n")

out.append(f"/// All {n_interactions} source-derived main-group interaction records.")
out.append(
    f"pub static UNIFAC_INTERACTIONS: [UnifacInteraction; {n_interactions}] = ["
)
for record in interactions:
    out.append("    UnifacInteraction {")
    out.append(f'        main_group_i: {record["main_group_i"]},')
    out.append(f'        main_group_j: {record["main_group_j"]},')
    out.append(f'        aij: {record["aij"]:.10e},')
    out.append(f'        bij: {record["bij"]:.10e},')
    out.append(f'        cij: {record["cij"]:.10e},')
    out.append(f'        data_source: "{record["data_source"]}",')
    out.append("    },")
out.append("];\n")

out.append("""/// Look up a component record by its identifier.
pub fn lookup_unifac_component(uid: &str) -> Option<&'static UnifacComponent> {
    UNIFAC_COMPONENTS.iter().find(|c| c.uid == uid)
}

/// Look up a main-group interaction record, mirroring the source scan that
/// returns the first matching (main group i, main group j) pair.
pub fn lookup_unifac_interaction(
    main_group_i: u8,
    main_group_j: u8,
) -> Option<&'static UnifacInteraction> {
    UNIFAC_INTERACTIONS
        .iter()
        .find(|r| r.main_group_i == main_group_i && r.main_group_j == main_group_j)
}
""")

target = ROOT / "rust/src/unifac_data.rs"
target.write_text("\n".join(out), encoding="utf-8")
print(f"Generated {target}.")
