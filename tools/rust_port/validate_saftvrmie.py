"""Validate SAFT-VR Mie and SAFT-VRQ Mie implementation against exact-source reference DLL."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/rust_port"))
from reference import load

# Load reference library
load(str(ROOT / ".toolchains/gcc14/mingw64/bin"))
from thermopack.saftvrmie import saftvrmie

def probe_pure(ident, temperatures, densities, ref="Default"):
    eos = saftvrmie(ident, parameter_reference=ref)
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
            a_hs = eos.a_hs(T, V, n)
            a1 = eos.a1(T, V, n)
            a2 = eos.a2(T, V, n)
            a3 = eos.a3(T, V, n)
            results.append({
                "fluid": ident,
                "ref": ref,
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
                "a_hs": a_hs,
                "a1": a1,
                "a2": a2,
                "a3": a3,
            })
    return results

def probe_mixture(comps, T, V, n, ref="Default"):
    eos = saftvrmie(comps, parameter_reference=ref)
    P, dpdt, dpdv, dpdn = eos.pressure_tv(T, V, n, dpdt=True, dpdv=True, dpdn=True)
    lnphi, = eos.thermo_tvp(T, V, n)
    h, = eos.enthalpy_tv(T, V, n)
    s, = eos.entropy_tv(T, V, n)
    c = eos.speed_of_sound_tv(T, V, n)
    return {
        "comps": comps,
        "ref": ref,
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
        "c": c,
    }

if __name__ == "__main__":
    print("Probing pure CO2...")
    res_co2 = probe_pure("CO2", [250.0, 300.0, 350.0], [50.0, 500.0, 5000.0, 15000.0])
    print(f"Generated {len(res_co2)} points for CO2")
    print("CO2 point 0:", res_co2[0])
    print("CO2 dense liquid point:", res_co2[-1])

    print("\nProbing pure C1...")
    res_c1 = probe_pure("C1", [150.0, 200.0, 300.0], [50.0, 2000.0, 12000.0])
    print(f"Generated {len(res_c1)} points for C1")
    print("C1 point 0:", res_c1[0])

    print("\nProbing mixture CO2,C1...")
    mix1 = probe_mixture("CO2,C1", 280.0, 0.0002, [0.4, 0.6])
    print("CO2,C1 mix:", mix1)

    print("\nProbing quantum fluid H2 (FH0/Default)...")
    res_h2 = probe_pure("H2", [30.0, 50.0, 100.0], [100.0, 2000.0, 10000.0])
    print(f"Generated {len(res_h2)} points for H2")
    print("H2 point 0:", res_h2[0])
