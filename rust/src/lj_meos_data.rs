//! Generated from src/multiparameter_lj.f90 by generate_lj_meos_data.py.
//! Copyright thermotools and Contributors. Apache-2.0.
#![allow(clippy::excessive_precision)]
use crate::multiparameter_data::{
    AncillarySatDensity, ExpExpTerm, ExpTerm, IdealGasHelmholtz, MultiparameterRecord, PolTerm,
};
#[derive(Debug, Clone, Copy)]
pub(crate) struct Gas {
    pub component: &'static str,
    pub truncated: bool,
    pub epsilon: f64,
    pub sigma: f64,
}
pub(crate) const GASES: [Gas; 8] = [
    Gas {
        component: "NE",
        truncated: false,
        epsilon: 33.921,
        sigma: 2.801e-10,
    },
    Gas {
        component: "AR",
        truncated: false,
        epsilon: 116.79,
        sigma: 3.3952e-10,
    },
    Gas {
        component: "KR",
        truncated: false,
        epsilon: 162.58,
        sigma: 3.6274e-10,
    },
    Gas {
        component: "XE",
        truncated: false,
        epsilon: 226.14,
        sigma: 3.949e-10,
    },
    Gas {
        component: "NE",
        truncated: true,
        epsilon: 39.83,
        sigma: 2.8e-10,
    },
    Gas {
        component: "AR",
        truncated: true,
        epsilon: 137.9,
        sigma: 3.3916e-10,
    },
    Gas {
        component: "KR",
        truncated: true,
        epsilon: 191.52,
        sigma: 3.6233e-10,
    },
    Gas {
        component: "XE",
        truncated: true,
        epsilon: 274.86,
        sigma: 3.946e-10,
    },
];
static LJ_POL: [PolTerm; 6] = [
    PolTerm {
        n: 0.005208073,
        t: 1.0,
        d: 4,
    },
    PolTerm {
        n: 2.186252,
        t: 0.32,
        d: 1,
    },
    PolTerm {
        n: -2.161016,
        t: 0.505,
        d: 1,
    },
    PolTerm {
        n: 1.4527,
        t: 0.672,
        d: 2,
    },
    PolTerm {
        n: -2.041792,
        t: 0.843,
        d: 2,
    },
    PolTerm {
        n: 0.18695286,
        t: 0.898,
        d: 3,
    },
];
static LJ_EXP: [ExpTerm; 6] = [
    ExpTerm {
        n: -0.090988445,
        t: 1.294,
        d: 5,
        l: 1,
    },
    ExpTerm {
        n: -0.4974561,
        t: 2.59,
        d: 2,
        l: 2,
    },
    ExpTerm {
        n: 0.10901431,
        t: 1.786,
        d: 2,
        l: 1,
    },
    ExpTerm {
        n: -0.80055922,
        t: 2.77,
        d: 3,
        l: 2,
    },
    ExpTerm {
        n: -0.568839,
        t: 1.786,
        d: 1,
        l: 2,
    },
    ExpTerm {
        n: -0.6208625,
        t: 1.205,
        d: 1,
        l: 1,
    },
];
static LJ_EXPEXP: [ExpExpTerm; 11] = [
    ExpExpTerm {
        n: -1.4667177,
        t: 2.83,
        d: 1,
        eta: 2.067,
        beta: 0.625,
        gamma: 0.71,
        epsilon: 0.2053,
    },
    ExpExpTerm {
        n: 1.891469,
        t: 2.548,
        d: 1,
        eta: 1.522,
        beta: 0.638,
        gamma: 0.86,
        epsilon: 0.409,
    },
    ExpExpTerm {
        n: -0.1383701,
        t: 4.65,
        d: 2,
        eta: 8.82,
        beta: 3.91,
        gamma: 1.94,
        epsilon: 0.6,
    },
    ExpExpTerm {
        n: -0.3869645,
        t: 1.385,
        d: 3,
        eta: 1.722,
        beta: 0.156,
        gamma: 1.48,
        epsilon: 1.203,
    },
    ExpExpTerm {
        n: 0.1265702,
        t: 1.46,
        d: 3,
        eta: 0.679,
        beta: 0.157,
        gamma: 1.49,
        epsilon: 1.829,
    },
    ExpExpTerm {
        n: 0.605781,
        t: 1.351,
        d: 2,
        eta: 1.883,
        beta: 0.153,
        gamma: 1.945,
        epsilon: 1.397,
    },
    ExpExpTerm {
        n: 1.179189,
        t: 0.66,
        d: 1,
        eta: 3.925,
        beta: 1.16,
        gamma: 3.02,
        epsilon: 1.39,
    },
    ExpExpTerm {
        n: -0.47732679,
        t: 1.496,
        d: 2,
        eta: 2.461,
        beta: 1.73,
        gamma: 1.11,
        epsilon: 0.539,
    },
    ExpExpTerm {
        n: -9.9218575,
        t: 1.83,
        d: 3,
        eta: 28.2,
        beta: 383.0,
        gamma: 1.17,
        epsilon: 0.934,
    },
    ExpExpTerm {
        n: -0.5747932,
        t: 1.616,
        d: 1,
        eta: 0.753,
        beta: 0.112,
        gamma: 1.33,
        epsilon: 2.369,
    },
    ExpExpTerm {
        n: 0.003772923,
        t: 4.97,
        d: 1,
        eta: 0.82,
        beta: 0.119,
        gamma: 0.24,
        epsilon: 2.43,
    },
];
static LJ_N_LIQSAT: [f64; 5] = [1.362, 2.093, -2.11, 0.329, 1.41];
static LJ_EXPO_LIQSAT: [f64; 5] = [0.313, 0.94, 1.63, 17.0, 2.4];
static LJ_N_VAPSAT: [f64; 6] = [-6.9655, -103.31, -2.0325, -44.481, -18.463, -0.0002607];
static LJ_EXPO_VAPSAT: [f64; 6] = [1.32, 19.24, 0.36, 8.78, 4.04, 41.6];
pub(crate) const LJ_LIMITS: (f64, f64) = (9.0, 65.0);
pub(crate) const LJ_TEMPLATE: MultiparameterRecord = MultiparameterRecord {
    ident: "LJ",
    tc: 1.32,
    pc: 0.0,
    rc: 0.31,
    acf: -0.042286,
    molar_mass: 1.0,
    t_triple: 0.661,
    p_triple: 0.0,
    rholiq_triple: 0.0,
    rhovap_triple: 0.0,
    rgas_fit: 1.0,
    ideal: IdealGasHelmholtz::Einstein {
        c_tau: 1.5,
        a: [6.262265814, -1.515151515],
        b: &[],
        v: &[],
    },
    pol_terms: &LJ_POL,
    exp_terms: &LJ_EXP,
    expexp_terms: &LJ_EXPEXP,
    ancillary: AncillarySatDensity::LemmonSpan {
        n_liq: &LJ_N_LIQSAT,
        expo_liq: &LJ_EXPO_LIQSAT,
        n_vap: &LJ_N_VAPSAT,
        expo_vap: &LJ_EXPO_VAPSAT,
    },
};
static LJTS_POL: [PolTerm; 6] = [
    PolTerm {
        n: 0.015606084,
        t: 1.0,
        d: 4,
    },
    PolTerm {
        n: 1.7917527,
        t: 0.304,
        d: 1,
    },
    PolTerm {
        n: -1.9613228,
        t: 0.583,
        d: 1,
    },
    PolTerm {
        n: 1.3045604,
        t: 0.662,
        d: 2,
    },
    PolTerm {
        n: -1.8117673,
        t: 0.87,
        d: 2,
    },
    PolTerm {
        n: 0.15483997,
        t: 0.87,
        d: 3,
    },
];
static LJTS_EXP: [ExpTerm; 6] = [
    ExpTerm {
        n: -0.094885204,
        t: 1.25,
        d: 5,
        l: 1,
    },
    ExpTerm {
        n: -0.20092412,
        t: 3.0,
        d: 2,
        l: 2,
    },
    ExpTerm {
        n: 0.11639644,
        t: 1.7,
        d: 2,
        l: 1,
    },
    ExpTerm {
        n: -0.50607364,
        t: 2.4,
        d: 3,
        l: 2,
    },
    ExpTerm {
        n: -0.58422807,
        t: 1.96,
        d: 1,
        l: 2,
    },
    ExpTerm {
        n: -0.47510982,
        t: 1.286,
        d: 1,
        l: 1,
    },
];
static LJTS_EXPEXP: [ExpExpTerm; 9] = [
    ExpExpTerm {
        n: 0.0094333106,
        t: 3.6,
        d: 1,
        eta: 4.7,
        beta: 20.0,
        gamma: 1.0,
        epsilon: 0.55,
    },
    ExpExpTerm {
        n: 0.30444628,
        t: 2.08,
        d: 1,
        eta: 1.92,
        beta: 0.77,
        gamma: 0.5,
        epsilon: 0.7,
    },
    ExpExpTerm {
        n: -0.0010820946,
        t: 5.24,
        d: 2,
        eta: 2.7,
        beta: 0.5,
        gamma: 0.8,
        epsilon: 2.0,
    },
    ExpExpTerm {
        n: -0.099693391,
        t: 0.96,
        d: 3,
        eta: 1.49,
        beta: 0.8,
        gamma: 1.5,
        epsilon: 1.14,
    },
    ExpExpTerm {
        n: 0.0091193522,
        t: 1.36,
        d: 3,
        eta: 0.65,
        beta: 0.4,
        gamma: 0.7,
        epsilon: 1.2,
    },
    ExpExpTerm {
        n: 0.12970543,
        t: 1.655,
        d: 2,
        eta: 1.73,
        beta: 0.43,
        gamma: 1.6,
        epsilon: 1.31,
    },
    ExpExpTerm {
        n: 0.02303603,
        t: 0.9,
        d: 1,
        eta: 3.7,
        beta: 8.0,
        gamma: 1.3,
        epsilon: 1.14,
    },
    ExpExpTerm {
        n: -0.082671073,
        t: 0.86,
        d: 2,
        eta: 1.9,
        beta: 3.3,
        gamma: 0.6,
        epsilon: 0.53,
    },
    ExpExpTerm {
        n: -2.2497821,
        t: 3.95,
        d: 3,
        eta: 13.2,
        beta: 114.0,
        gamma: 1.3,
        epsilon: 0.96,
    },
];
static LJTS_N_LIQSAT: [f64; 4] = [1.45, -0.172, -0.298, 0.295];
static LJTS_EXPO_LIQSAT: [f64; 4] = [0.334, 0.667, 1.25, 1.92];
static LJTS_N_VAPSAT: [f64; 4] = [1.59809, -0.09975, -0.4774, -2.33736];
static LJTS_EXPO_VAPSAT: [f64; 4] = [1.0, 1.5, 5.94, 0.41452];
pub(crate) const LJTS_LIMITS: (f64, f64) = (11.0, 70.0);
pub(crate) const LJTS_TEMPLATE: MultiparameterRecord = MultiparameterRecord {
    ident: "LJTS",
    tc: 1.086,
    pc: 0.0,
    rc: 0.319,
    acf: 0.0244417,
    molar_mass: 1.0,
    t_triple: 0.56,
    p_triple: 0.0,
    rholiq_triple: 0.0,
    rhovap_triple: 0.0,
    rgas_fit: 1.0,
    ideal: IdealGasHelmholtz::Einstein {
        c_tau: 1.5,
        a: [0., -2.5 * 0.8 / 1.086],
        b: &[],
        v: &[],
    },
    pol_terms: &LJTS_POL,
    exp_terms: &LJTS_EXP,
    expexp_terms: &LJTS_EXPEXP,
    ancillary: AncillarySatDensity::LemmonSpan {
        n_liq: &LJTS_N_LIQSAT,
        expo_liq: &LJTS_EXPO_LIQSAT,
        n_vap: &LJTS_N_VAPSAT,
        expo_vap: &LJTS_EXPO_VAPSAT,
    },
};
