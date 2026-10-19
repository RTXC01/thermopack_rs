"""Validate PC-SAFT and sPC-SAFT implementations against exact-source reference DLL."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/rust_port"))
from reference import load

# Load reference library
load(".toolchains/gcc14/mingw64/bin")
from thermopack.pcsaft import pcsaft

LIQPH = 1
VAPPH = 2

def probe_pure(ident, temperatures, densities):
    eos = pcsaft(ident)
    results = []
    for T in temperatures:
        for rho in densities:
            V = 1.0 / rho
            n = [1.0]
            P, dpdt, dpdv, dpdn = eos.pressure_tv(T, V, n, dpdt=True, dpdv=True, dpdn=True)
            u, = eos.internal_energy_tv(T, V, n)
            h, = eos.enthalpy_tv(T, V, n)
            s, = eos.entropy_tv(T, V, n)
            a, = eos.helmholtz_tv(T, V, n)
            mu, = eos.chemical_potential_tv(T, V, n)
            lnphi, = eos.thermo_tvp(T, V, n)
            c = eos.speed_of_sound_tv(T, V, n)
            results.append({
                "fluid": ident,
                "T": T,
                "rho": rho,
                "P": P,
                "dpdt": dpdt,
                "dpdv": dpdv,
                "dpdn": float(dpdn[0]),
                "u": u,
                "h": h,
                "s": s,
                "a": a,
                "mu": float(mu[0]),
                "lnphi": float(lnphi[0]),
                "c": c,
            })
    return results

def probe_mixture(comps, simplified, T, V, n):
    eos = pcsaft(comps, simplified=simplified)
    P, dpdt, dpdv, dpdn = eos.pressure_tv(T, V, n, dpdt=True, dpdv=True, dpdn=True)
    lnphi, = eos.thermo_tvp(T, V, n)
    h, = eos.enthalpy_tv(T, V, n)
    s, = eos.entropy_tv(T, V, n)
    return {
        "comps": comps,
        "simplified": simplified,
        "T": T,
        "V": V,
        "n": n,
        "P": P,
        "dpdt": dpdt,
        "dpdv": dpdv,
        "dpdn": [float(x) for x in dpdn],
        "lnphi": [float(x) for x in lnphi],
        "h": h,
        "s": s,
    }

if __name__ == "__main__":
    print("Testing pure C3...")
    res_c3 = probe_pure("C3", [250.0, 300.0, 350.0], [50.0, 500.0, 10000.0])
    print(f"Generated {len(res_c3)} state points for C3")
    print("Sample C3 point 0:", res_c3[0])

    print("\nTesting mixture CO2,C1...")
    mix1 = probe_mixture("CO2,C1", False, 280.0, 0.002, [0.3, 0.7])
    print("CO2,C1 orig:", mix1)
    mix2 = probe_mixture("CO2,C1", True, 280.0, 0.002, [0.3, 0.7])
    print("CO2,C1 simp:", mix2)
