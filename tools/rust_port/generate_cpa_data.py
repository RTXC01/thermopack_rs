"""Generate rust/src/cpa_data.rs from compdatadb.f90 and mixdatadb.f90."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
compdatadb = ROOT / "src/compdatadb.f90"
mixdatadb = ROOT / "src/mixdatadb.f90"

with open(compdatadb, "r", encoding="utf-8") as f:
    lines = f.readlines()

pure_entries = []
i = 0
while i < len(lines):
    line = lines[i]
    if "type(CPAdata), parameter :: cpa" in line:
        cpa_id = line.split("::")[1].split("=")[0].strip()
        body_lines = []
        i += 1
        while i < len(lines):
            body_line = lines[i]
            body_lines.append(body_line)
            if ")" in body_line and ("compName" in "".join(body_lines) or "CPAdata" in "".join(body_lines)):
                if body_line.strip().startswith(")") or body_line.strip().endswith(")"):
                    break
            i += 1
        body = " ".join(body_lines).replace("&", "").replace("\n", " ")
        comp = re.search(r'compName\s*=\s*"([^"]+)"', body)
        eosid = re.search(r'eosid\s*=\s*"([^"]+)"', body)
        ref = re.search(r'ref\s*=\s*"([^"]+)"', body)
        a0 = re.search(r'a0\s*=\s*([\d\.e\+\-]+)', body)
        b = re.search(r'b\s*=\s*([\d\.e\+\-]+)', body)
        eps = re.search(r'eps\s*=\s*([\d\.e\+\-]+)', body)
        beta = re.search(r'beta\s*=\s*([\d\.e\+\-]+)', body)
        alpha = re.search(r'alphaParams\s*=\s*\((.*?)\)', body)
        alphacorr = re.search(r'alphacorridx\s*=\s*([^\s,]+)', body)
        scheme = re.search(r'assoc_scheme\s*=\s*([^\s,]+)', body)
        simp_rdf = re.search(r'simplified_rdf\s*=\s*([^\s,]+)', body)

        pure_entries.append({
            "id": cpa_id,
            "comp": comp.group(1) if comp else "",
            "eosid": eosid.group(1) if eosid else "",
            "ref": ref.group(1) if ref else "",
            "a0": float(a0.group(1)) if a0 else 0.0,
            "b": float(b.group(1)) if b else 0.0,
            "eps": float(eps.group(1)) if eps else 0.0,
            "beta": float(beta.group(1)) if beta else 0.0,
            "c1": float(alpha.group(1).replace("/", "").split(",")[0].strip()) if alpha else 0.0,
            "scheme": scheme.group(1) if scheme else "no_assoc",
            "simplified_rdf": simp_rdf.group(1).lower().startswith(".t") if simp_rdf else True,
        })
    i += 1

# Binary entries from mixdatadb.f90
with open(mixdatadb, "r", encoding="utf-8") as f:
    mix_lines = f.readlines()

binary_entries = []
i = 0
while i < len(mix_lines):
    line = mix_lines[i]
    if "type(CPAkijdata), parameter :: cpa" in line:
        body_lines = []
        i += 1
        while i < len(mix_lines):
            body_line = mix_lines[i]
            body_lines.append(body_line)
            if ")" in body_line and ("uid1" in "".join(body_lines) or "CPAkijdata" in "".join(body_lines)):
                if body_line.strip().startswith(")") or body_line.strip().endswith(")"):
                    break
            i += 1
        body = " ".join(body_lines).replace("&", "").replace("\n", " ")
        eosid = re.search(r'eosid\s*=\s*"([^"]+)"', body)
        ref = re.search(r'ref\s*=\s*"([^"]+)"', body)
        uid1 = re.search(r'uid1\s*=\s*"([^"]+)"', body)
        uid2 = re.search(r'uid2\s*=\s*"([^"]+)"', body)
        kij_a = re.search(r'kij_a\s*=\s*([\d\.e\+\-]+)', body)
        kij_eps = re.search(r'kij_eps\s*=\s*([\d\.e\+\-]+)', body)
        kij_beta = re.search(r'kij_beta\s*=\s*([\d\.e\+\-]+)', body)
        eps_comb = re.search(r'eps_comb_rule\s*=\s*([^\s,]+)', body)
        beta_comb = re.search(r'beta_comb_rule\s*=\s*([^\s,]+)', body)

        binary_entries.append({
            "eosid": eosid.group(1) if eosid else "",
            "ref": ref.group(1) if ref else "",
            "uid1": uid1.group(1) if uid1 else "",
            "uid2": uid2.group(1) if uid2 else "",
            "kij_a": float(kij_a.group(1)) if kij_a else 0.0,
            "kij_eps": float(kij_eps.group(1)) if kij_eps else 0.0,
            "kij_beta": float(kij_beta.group(1)) if kij_beta else 0.0,
            "eps_comb": eps_comb.group(1) if eps_comb else "ariComb",
            "beta_comb": beta_comb.group(1) if beta_comb else "geoComb",
        })
    i += 1

# Generate Rust code
code = []
code.append("//! CPA (Cubic-Plus-Association) parameter database.")
code.append("//! Automatically generated from compdatadb.f90 and mixdatadb.f90.")
code.append("//! Copyright thermotools and Contributors. Apache-2.0.\n")
code.append("#![allow(clippy::excessive_precision)]\n")

code.append("/// Association scheme classification.")
code.append("#[derive(Debug, Clone, Copy, PartialEq, Eq)]")
code.append("pub enum AssocScheme {")
code.append("    NoAssoc,")
code.append("    Scheme1,")
code.append("    Scheme2A,")
code.append("    Scheme2B,")
code.append("    Scheme2C,")
code.append("    Scheme3A,")
code.append("    Scheme3B,")
code.append("    Scheme4A,")
code.append("    Scheme4B,")
code.append("    Scheme4C,")
code.append("    Scheme1ea,")
code.append("}\n")

code.append("impl AssocScheme {")
code.append("    /// Number of association sites.")
code.append("    pub const fn num_sites(self) -> usize {")
code.append("        match self {")
code.append("            AssocScheme::NoAssoc => 0,")
code.append("            AssocScheme::Scheme1 | AssocScheme::Scheme1ea => 1,")
code.append("            AssocScheme::Scheme2A | AssocScheme::Scheme2B | AssocScheme::Scheme2C => 2,")
code.append("            AssocScheme::Scheme3A | AssocScheme::Scheme3B => 3,")
code.append("            AssocScheme::Scheme4A | AssocScheme::Scheme4B | AssocScheme::Scheme4C => 4,")
code.append("        }")
code.append("    }")
code.append("}\n")

code.append("/// Pure component CPA parameters.")
code.append("#[derive(Debug, Clone)]")
code.append("pub struct CpaPureRecord {")
code.append("    pub comp_name: &'static str,")
code.append("    pub eos: &'static str,")
code.append("    pub reference: &'static str,")
code.append("    pub a0: f64,")
code.append("    pub b: f64,")
code.append("    pub eps: f64,")
code.append("    pub beta: f64,")
code.append("    pub c1: f64,")
code.append("    pub scheme: AssocScheme,")
code.append("    pub simplified_rdf: bool,")
code.append("}\n")

code.append("/// Binary interaction CPA parameters.")
code.append("#[derive(Debug, Clone)]")
code.append("pub struct CpaBinaryRecord {")
code.append("    pub comp1: &'static str,")
code.append("    pub comp2: &'static str,")
code.append("    pub eos: &'static str,")
code.append("    pub reference: &'static str,")
code.append("    pub kij_a: f64,")
code.append("    pub kij_eps: f64,")
code.append("    pub kij_beta: f64,")
code.append("    pub eps_comb_arithmetic: bool,")
code.append("    pub beta_comb_arithmetic: bool,")
code.append("}\n")

def map_scheme(s):
    s = s.lower()
    if "no_assoc" in s:
        return "AssocScheme::NoAssoc"
    if "1ea" in s:
        return "AssocScheme::Scheme1ea"
    if "scheme_1" in s:
        return "AssocScheme::Scheme1"
    if "scheme_2a" in s:
        return "AssocScheme::Scheme2A"
    if "scheme_2b" in s:
        return "AssocScheme::Scheme2B"
    if "scheme_2c" in s:
        return "AssocScheme::Scheme2C"
    if "scheme_3a" in s:
        return "AssocScheme::Scheme3A"
    if "scheme_3b" in s:
        return "AssocScheme::Scheme3B"
    if "scheme_4a" in s:
        return "AssocScheme::Scheme4A"
    if "scheme_4b" in s:
        return "AssocScheme::Scheme4B"
    if "scheme_4c" in s:
        return "AssocScheme::Scheme4C"
    return "AssocScheme::NoAssoc"

code.append(f"pub const CPA_PURE_RECORDS: [CpaPureRecord; {len(pure_entries)}] = [")
for e in pure_entries:
    code.append("    CpaPureRecord {")
    code.append(f'        comp_name: "{e["comp"]}",')
    code.append(f'        eos: "{e["eosid"]}",')
    code.append(f'        reference: "{e["ref"]}",')
    code.append(f'        a0: {e["a0"]:.16e},')
    code.append(f'        b: {e["b"]:.16e},')
    code.append(f'        eps: {e["eps"]:.16e},')
    code.append(f'        beta: {e["beta"]:.16e},')
    code.append(f'        c1: {e["c1"]:.16e},')
    code.append(f'        scheme: {map_scheme(e["scheme"])},')
    code.append(f'        simplified_rdf: {"true" if e["simplified_rdf"] else "false"},')
    code.append("    },")
code.append("];\n")

code.append(f"pub const CPA_BINARY_RECORDS: [CpaBinaryRecord; {len(binary_entries)}] = [")
for b in binary_entries:
    code.append("    CpaBinaryRecord {")
    code.append(f'        comp1: "{b["uid1"]}",')
    code.append(f'        comp2: "{b["uid2"]}",')
    code.append(f'        eos: "{b["eosid"]}",')
    code.append(f'        reference: "{b["ref"]}",')
    code.append(f'        kij_a: {b["kij_a"]:.16e},')
    code.append(f'        kij_eps: {b["kij_eps"]:.16e},')
    code.append(f'        kij_beta: {b["kij_beta"]:.16e},')
    code.append(f'        eps_comb_arithmetic: {"true" if "ari" in b["eps_comb"].lower() else "false"},')
    code.append(f'        beta_comb_arithmetic: {"true" if "ari" in b["beta_comb"].lower() else "false"},')
    code.append("    },")
code.append("];\n")

code.append("""/// Lookup pure CPA parameters by component identifier and optional reference.
pub fn lookup_cpa_pure(comp: &str, preferred_ref: Option<&str>) -> Option<&'static CpaPureRecord> {
    let comp_up = comp.to_ascii_uppercase();
    if let Some(pref) = preferred_ref {
        let pref_up = pref.to_ascii_uppercase();
        if let Some(r) = CPA_PURE_RECORDS.iter().find(|r| r.comp_name.eq_ignore_ascii_case(&comp_up) && r.reference.to_ascii_uppercase().contains(&pref_up)) {
            return Some(r);
        }
    }
    // Default fallback
    CPA_PURE_RECORDS.iter().find(|r| r.comp_name.eq_ignore_ascii_case(&comp_up) && r.reference.starts_with("Default"))
        .or_else(|| CPA_PURE_RECORDS.iter().find(|r| r.comp_name.eq_ignore_ascii_case(&comp_up)))
}

/// Lookup binary CPA parameters for a pair of components.
pub fn lookup_cpa_binary(comp1: &str, comp2: &str) -> Option<&'static CpaBinaryRecord> {
    let c1 = comp1.to_ascii_uppercase();
    let c2 = comp2.to_ascii_uppercase();
    CPA_BINARY_RECORDS.iter().find(|r| {
        (r.comp1.eq_ignore_ascii_case(&c1) && r.comp2.eq_ignore_ascii_case(&c2))
            || (r.comp1.eq_ignore_ascii_case(&c2) && r.comp2.eq_ignore_ascii_case(&c1))
    })
}
""")

target = ROOT / "rust/src/cpa_data.rs"
target.write_text("\n".join(code), encoding="utf-8")
print(f"Generated {target} with {len(pure_entries)} pure records and {len(binary_entries)} binary records.")
