//! Dufal association integral and its TVN derivatives at unit total moles.
use super::{N_A, PreparedSaftVrMie};
use crate::Result;
use crate::association::{AssocRdfFirst, AssocRdfSecond, AssocRdfValue, AssociationRdfProvider};
use crate::rdf::RdfDerivatives;

pub(super) struct MieRdf<'a> {
    state: &'a PreparedSaftVrMie<'a>,
    density: f64,
    reduced: f64,
}

impl<'a> MieRdf<'a> {
    pub(super) fn new(state: &'a PreparedSaftVrMie<'a>, density: f64) -> Self {
        let mut reduced = 0.0;
        for (i, xi) in state.xs.iter().enumerate() {
            for (j, xj) in state.xs.iter().enumerate() {
                reduced += xi * xj * state.model.sigma_ij[i * state.xs.len() + j].powi(3);
            }
        }
        Self {
            state,
            density,
            reduced: reduced * N_A * density * state.m_bar,
        }
    }

    // Horner evaluation of sum a_ij rho^i T^j. Differentiate the
    // recurrences alongside the value only to the requested order.
    fn integral<const ORDER: u8>(&self, i: usize, j: usize) -> [f64; 6] {
        let model = self.state.model;
        let idx = i * model.num_components() + j;
        let eps = model.eps_divk_ij[idx];
        let t = self.state.temperature / eps;
        let r = self.reduced;
        let mut coef = model.integral_coefficients[idx].iter().rev();
        let (mut g, mut gr, mut gt, mut grr, mut grt, mut gtt) = (0., 0., 0., 0., 0., 0.);
        for degree in 0..=10 {
            let (mut b, mut bt, mut btt) = (0., 0., 0.);
            for _ in 0..=degree {
                if ORDER >= 2 {
                    btt = btt * t + 2. * bt;
                }
                if ORDER >= 1 {
                    bt = bt * t + b;
                }
                b = b * t + coef.next().unwrap();
            }
            if ORDER >= 2 {
                grr = grr * r + 2. * gr;
                grt = grt * r + gt;
                gtt = gtt * r + btt;
            }
            if ORDER >= 1 {
                gr = gr * r + g;
                gt = gt * r + bt;
            }
            g = g * r + b;
        }
        [g, gr, gt / eps, grr, grt / eps, gtt / (eps * eps)]
    }

    fn reduced_moles(&self) -> Vec<f64> {
        let p = self.state;
        let nc = p.model.num_components();
        (0..nc)
            .map(|k| {
                let cross: f64 =
                    p.xs.iter()
                        .enumerate()
                        .map(|(l, x)| x * p.model.sigma_ij[k * nc + l].powi(3))
                        .sum();
                p.model.components[k].m * (2. * N_A * self.density * cross - self.reduced / p.m_bar)
            })
            .collect()
    }
}

impl AssociationRdfProvider for MieRdf<'_> {
    fn pair_dependent(&self) -> bool {
        true
    }

    fn rdf_value(&self, i: usize, j: usize) -> Result<AssocRdfValue> {
        Ok(AssocRdfValue {
            g: self.integral::<0>(i, j)[0],
            covol: N_A,
        })
    }

    fn rdf_first(&self, i: usize, j: usize) -> Result<AssocRdfFirst> {
        let [g, gr, gt, ..] = self.integral::<1>(i, j);
        Ok(AssocRdfFirst {
            g,
            g_t: gt,
            g_v: -gr * self.reduced * self.density,
            g_n: self.reduced_moles().into_iter().map(|v| gr * v).collect(),
            covol: N_A,
            covol_t: 0.,
        })
    }

    fn rdf_second(&self, i: usize, j: usize) -> Result<AssocRdfSecond> {
        let p = self.state;
        let nc = p.model.num_components();
        let [g, gr, gt, grr, grt, gtt] = self.integral::<2>(i, j);
        let rn = self.reduced_moles();
        let rv = -self.reduced * self.density;
        let mut nn = Vec::with_capacity(nc * nc);
        for k in 0..nc {
            let mk = p.model.components[k].m;
            for l in 0..nc {
                let ml = p.model.components[l].m;
                let rkl =
                    (2. * N_A * self.density * mk * ml * p.model.sigma_ij[k * nc + l].powi(3)
                        - mk * rn[l]
                        - ml * rn[k])
                        / p.m_bar;
                nn.push(grr * rn[k] * rn[l] + gr * rkl);
            }
        }
        Ok(AssocRdfSecond {
            covol: N_A,
            covol_t: 0.,
            covol_tt: 0.,
            rdf: RdfDerivatives {
                g,
                temperature: gt,
                volume: gr * rv,
                moles: rn.iter().map(|d| gr * d).collect(),
                temperature_second: gtt,
                temperature_volume: grt * rv,
                temperature_moles: rn.iter().map(|d| grt * d).collect(),
                volume_second: grr * rv * rv - 2. * gr * rv * self.density,
                volume_moles: rn
                    .iter()
                    .map(|d| (grr * rv - gr * self.density) * d)
                    .collect(),
                moles_second: nn,
            },
        })
    }
}
