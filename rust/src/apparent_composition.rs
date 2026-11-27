//! Apparent-composition electrolyte transformation utilities.
//!
//! Ported from `src/apparent_composition.f90` (`module apparent_compostion`).
//! Maps between apparent component mole numbers `n` (solvents + neutral salts,
//! length `nc`) and real species mole numbers `ne` (solvents + dissociated ions,
//! length `nce`), transforms real Helmholtz derivatives to apparent composition
//! derivatives, and converts real fugacity coefficients and their `(T, P, n)`
//! derivatives into apparent fugacity coefficients and modified fugacities.

use crate::error::positive;
use crate::{Error, Result};

/// Minimum mole-number floor used in apparent composition derivative denominators,
/// matching `thermopack_constants::min_mol_num`.
pub const MIN_MOL_NUM: f64 = 1.0e-150;

/// Apparent composition stoichiometry container (`apparent_container`).
///
/// The first `ncsym` components are identical in apparent and real bases
/// (neutral solvents/non-electrolytes). Apparent components `ncsym..nc`
/// dissociate into real species `ncsym..nce` according to the stoichiometric
/// matrix `v_stoich[i * nce + j]`, giving the number of moles of real species
/// `j` per mole of apparent component `i`.
#[derive(Debug, Clone, PartialEq)]
pub struct ApparentComposition {
    nc: usize,
    nce: usize,
    ncsym: usize,
    v_stoich: Vec<f64>,
    v_sum: Vec<f64>,
}

/// Real-to-apparent transformed Helmholtz derivatives (`real_to_apparent_differentials`).
#[derive(Debug, Clone, PartialEq)]
pub struct ApparentDifferentials {
    /// Apparent first mole derivatives `F_n` (length `nc`).
    pub f_n: Option<Vec<f64>>,
    /// Apparent mixed `(T, n)` derivatives `F_Tn` (length `nc`).
    pub f_tn: Option<Vec<f64>>,
    /// Apparent mixed `(V, n)` derivatives `F_Vn` (length `nc`).
    pub f_vn: Option<Vec<f64>>,
    /// Apparent second mole derivatives `F_nn` (row-major `nc * nc`).
    pub f_nn: Option<Vec<f64>>,
}

/// Apparent logarithmic fugacity coefficients and optional `(T, P, n)` derivatives
/// (`TP_lnfug_apparent`).
#[derive(Debug, Clone, PartialEq)]
pub struct ApparentFugacityResult {
    /// Apparent `ln(phi_i)` (length `nc`).
    pub ln_fugacity: Vec<f64>,
    /// Apparent `d ln(phi_i) / dT` (length `nc`).
    pub dt: Option<Vec<f64>>,
    /// Apparent `d ln(phi_i) / dP` (length `nc`).
    pub dp: Option<Vec<f64>>,
    /// Apparent `d ln(phi_i) / dn_j` (row-major `nc * nc`, index `i * nc + j`).
    pub dn: Option<Vec<f64>>,
}

impl ApparentComposition {
    /// Construct a new apparent composition container matching
    /// `apparent_constructor(nc, nce)`.
    ///
    /// Initializes `ncsym = nc`, `v_stoich = 0`, and `v_sum = 0`.
    pub fn new(nc: usize, nce: usize) -> Result<Self> {
        if nc == 0 || nce == 0 || nce < nc {
            return Err(Error::InvalidInput(
                "apparent composition requires 0 < nc <= nce",
            ));
        }
        Ok(Self {
            nc,
            nce,
            ncsym: nc,
            v_stoich: vec![0.0; nc * nce],
            v_sum: vec![0.0; nc],
        })
    }

    /// Construct a new apparent composition container with an explicit `ncsym`
    /// (number of non-dissociating solvent/molecular species).
    pub fn with_ncsym(nc: usize, nce: usize, ncsym: usize) -> Result<Self> {
        let mut container = Self::new(nc, nce)?;
        container.set_ncsym(ncsym)?;
        Ok(container)
    }

    /// Number of apparent components `nc`.
    #[inline]
    pub fn nc(&self) -> usize {
        self.nc
    }

    /// Number of real species `nce`.
    #[inline]
    pub fn nce(&self) -> usize {
        self.nce
    }

    /// Number of symmetric (unchanged molecular) species `ncsym`.
    #[inline]
    pub fn ncsym(&self) -> usize {
        self.ncsym
    }

    /// Set the number of symmetric (unchanged molecular) species `ncsym <= nc`.
    pub fn set_ncsym(&mut self, ncsym: usize) -> Result<()> {
        if ncsym > self.nc {
            return Err(Error::InvalidInput("ncsym must not exceed nc"));
        }
        self.ncsym = ncsym;
        Ok(())
    }

    /// Get stoichiometric coefficient `v_stoich(i, j)` (0-based indices).
    pub fn v_stoich(&self, i: usize, j: usize) -> Result<f64> {
        if i >= self.nc || j >= self.nce {
            return Err(Error::InvalidInput("stoichiometry index out of range"));
        }
        Ok(self.v_stoich[i * self.nce + j])
    }

    /// Row sums `v_sum[i] = sum_j v_stoich(i, j)`.
    #[inline]
    pub fn v_sum(&self) -> &[f64] {
        &self.v_sum
    }

    /// Set stoichiometric coefficient `v_stoich(i, j)` (0-based indices) and
    /// update `v_sum[i]`.
    pub fn set_v_stoich(&mut self, i: usize, j: usize, value: f64) -> Result<()> {
        if i >= self.nc || j >= self.nce {
            return Err(Error::InvalidInput("stoichiometry index out of range"));
        }
        if !value.is_finite() || value < 0.0 {
            return Err(Error::InvalidInput(
                "stoichiometric coefficient must be finite and non-negative",
            ));
        }
        self.v_stoich[i * self.nce + j] = value;
        self.update_v_sum();
        Ok(())
    }

    /// Recompute `v_sum[i] = sum_j v_stoich(i, j)` (`update_v_sum`).
    pub fn update_v_sum(&mut self) {
        for i in 0..self.nc {
            let row = &self.v_stoich[i * self.nce..(i + 1) * self.nce];
            self.v_sum[i] = row.iter().copied().sum();
        }
    }

    /// Map apparent mole numbers `n` (length `nc`) to real species mole numbers
    /// `ne` (length `nce`), matching `apparent_to_real_mole_numbers`.
    pub fn apparent_to_real_mole_numbers(&self, n: &[f64]) -> Result<Vec<f64>> {
        if n.len() != self.nc {
            return Err(Error::InvalidInput(
                "apparent mole vector length must equal nc",
            ));
        }
        if n.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("finite apparent mole numbers required"));
        }
        if self.nc == self.nce {
            return Ok(n.to_vec());
        }
        let mut ne = vec![0.0; self.nce];
        ne[..self.ncsym].copy_from_slice(&n[..self.ncsym]);
        for (i, ne_i) in ne.iter_mut().enumerate().skip(self.ncsym) {
            let mut sum = 0.0;
            for (j, &n_j) in n.iter().enumerate().skip(self.ncsym) {
                sum += self.v_stoich[j * self.nce + i] * n_j;
            }
            *ne_i = sum;
        }
        Ok(ne)
    }

    /// Map a real-species first-order derivative vector `fe_n` (length `nce`)
    /// to apparent composition derivatives `f_n` (length `nc`), matching
    /// `real_to_apparent_diff`.
    pub fn real_to_apparent_diff(&self, fe_n: &[f64]) -> Result<Vec<f64>> {
        if fe_n.len() != self.nce {
            return Err(Error::InvalidInput(
                "real derivative vector length must equal nce",
            ));
        }
        if self.nc == self.nce {
            return Ok(fe_n.to_vec());
        }
        let mut f_n = vec![0.0; self.nc];
        f_n[..self.ncsym].copy_from_slice(&fe_n[..self.ncsym]);
        for (i, out) in f_n.iter_mut().enumerate().skip(self.ncsym) {
            let mut sum = 0.0;
            for (j, &val) in fe_n.iter().enumerate().skip(self.ncsym) {
                sum += val * self.v_stoich[i * self.nce + j];
            }
            *out = sum;
        }
        Ok(f_n)
    }

    /// Map real Helmholtz composition derivatives `(Fe_n, Fe_Tn, Fe_Vn, Fe_nn)`
    /// into apparent composition derivatives (`real_to_apparent_differentials`).
    ///
    /// `fe_nn` is stored row-major with shape `nce * nce` (`fe_nn[r * nce + c]`
    /// corresponds to Fortran `Fe_nn(r+1, c+1)`). Output `f_nn` is row-major
    /// `nc * nc` (`f_nn[r * nc + c]` corresponds to Fortran `F_nn(r+1, c+1)`).
    pub fn real_to_apparent_differentials(
        &self,
        fe_n: Option<&[f64]>,
        fe_tn: Option<&[f64]>,
        fe_vn: Option<&[f64]>,
        fe_nn: Option<&[f64]>,
    ) -> Result<ApparentDifferentials> {
        let f_n = fe_n.map(|v| self.real_to_apparent_diff(v)).transpose()?;
        let f_tn = fe_tn.map(|v| self.real_to_apparent_diff(v)).transpose()?;
        let f_vn = fe_vn.map(|v| self.real_to_apparent_diff(v)).transpose()?;
        let f_nn = if let Some(fe_nn) = fe_nn {
            if fe_nn.len() != self.nce * self.nce {
                return Err(Error::InvalidInput(
                    "real Hessian length must equal nce * nce",
                ));
            }
            if self.nc == self.nce {
                Some(fe_nn.to_vec())
            } else {
                let mut out = vec![0.0; self.nc * self.nc];
                for k in 0..self.ncsym {
                    for i in 0..self.ncsym {
                        out[k * self.nc + i] = fe_nn[k * self.nce + i];
                    }
                }
                // Solvent-salt differentials
                for k in 0..self.ncsym {
                    for i in self.ncsym..self.nc {
                        let mut sum_ki = 0.0;
                        let mut sum_ik = 0.0;
                        for j in self.ncsym..self.nce {
                            let v_ij = self.v_stoich[i * self.nce + j];
                            sum_ki += fe_nn[k * self.nce + j] * v_ij;
                            sum_ik += fe_nn[j * self.nce + k] * v_ij;
                        }
                        out[k * self.nc + i] = sum_ki;
                        out[i * self.nc + k] = sum_ik;
                    }
                }
                // Salt-salt differentials:
                // Fortran line 146: F_nn(k,i) += Fe_nn(l,j) * v_stoich(k,j) * v_stoich(i,l)
                for k in self.ncsym..self.nc {
                    for i in self.ncsym..self.nc {
                        let mut sum = 0.0;
                        for j in self.ncsym..self.nce {
                            let v_kj = self.v_stoich[k * self.nce + j];
                            if v_kj == 0.0 {
                                continue;
                            }
                            for l in self.ncsym..self.nce {
                                let v_il = self.v_stoich[i * self.nce + l];
                                if v_il == 0.0 {
                                    continue;
                                }
                                sum += fe_nn[l * self.nce + j] * v_kj * v_il;
                            }
                        }
                        out[k * self.nc + i] = sum;
                    }
                }
                Some(out)
            }
        } else {
            None
        };
        Ok(ApparentDifferentials {
            f_n,
            f_tn,
            f_vn,
            f_nn,
        })
    }

    /// Convert real logarithmic fugacity coefficients and optional `(T, P, n)`
    /// derivatives to apparent composition (`TP_lnfug_apparent`).
    ///
    /// Unlike the Fortran baseline (`src/apparent_composition.f90:227`), this
    /// implementation skips zero stoichiometric entries (`v_stoich[i, j] == 0`)
    /// in the apparent salt fugacity sum so absent non-participating ions do not
    /// trigger `0.0 * (-inf) = NaN` in multi-salt mixtures.
    #[allow(clippy::too_many_arguments)]
    pub fn tp_lnfug_apparent(
        &self,
        n: &[f64],
        ne: &[f64],
        pressure: f64,
        lnfug_real: &[f64],
        dlnfug_dt_real: Option<&[f64]>,
        dlnfug_dp_real: Option<&[f64]>,
        dlnfug_dn_real: Option<&[f64]>,
    ) -> Result<ApparentFugacityResult> {
        positive(pressure, "pressure")?;
        if n.len() != self.nc || ne.len() != self.nce || lnfug_real.len() != self.nce {
            return Err(Error::InvalidInput(
                "input vector lengths must match apparent/real dimensions",
            ));
        }
        if self.nc == self.nce {
            if let Some(v) = dlnfug_dt_real
                && v.len() != self.nce
            {
                return Err(Error::InvalidInput("dlnfug_dt_real length mismatch"));
            }
            if let Some(v) = dlnfug_dp_real
                && v.len() != self.nce
            {
                return Err(Error::InvalidInput("dlnfug_dp_real length mismatch"));
            }
            if let Some(v) = dlnfug_dn_real
                && v.len() != self.nce * self.nce
            {
                return Err(Error::InvalidInput("dlnfug_dn_real length mismatch"));
            }
            return Ok(ApparentFugacityResult {
                ln_fugacity: lnfug_real.to_vec(),
                dt: dlnfug_dt_real.map(|v| v.to_vec()),
                dp: dlnfug_dp_real.map(|v| v.to_vec()),
                dn: dlnfug_dn_real.map(|v| v.to_vec()),
            });
        }

        let sumn: f64 = n.iter().copied().sum();
        let sumne: f64 = ne.iter().copied().sum();
        positive(sumn, "total apparent moles")?;
        positive(sumne, "total real moles")?;
        if n.iter().any(|&v| !v.is_finite() || v < 0.0)
            || ne.iter().any(|&v| !v.is_finite() || v < 0.0)
        {
            return Err(Error::InvalidInput(
                "mole numbers must be finite and non-negative",
            ));
        }

        let mut ln_fugacity = vec![0.0; self.nc];
        let log_ratio = (sumn / sumne).ln();
        for i in 0..self.ncsym {
            ln_fugacity[i] = lnfug_real[i] + log_ratio;
        }
        for i in self.ncsym..self.nc {
            let xi = n[i] / sumn;
            if xi > 0.0 {
                let mut sum = 0.0;
                for j in self.ncsym..self.nce {
                    let v_ij = self.v_stoich[i * self.nce + j];
                    if v_ij == 0.0 {
                        continue;
                    }
                    let xej = ne[j] / sumne;
                    if xej <= 0.0 {
                        return Err(Error::InvalidInput(
                            "participating real ion has zero mole fraction for present salt",
                        ));
                    }
                    sum += v_ij * (xej.ln() + lnfug_real[j]);
                }
                ln_fugacity[i] = sum - xi.ln();
            } else {
                ln_fugacity[i] = 0.0;
            }
        }

        let dt = dlnfug_dt_real
            .map(|v| self.real_to_apparent_diff(v))
            .transpose()?;
        let dp = dlnfug_dp_real
            .map(|v| self.real_to_apparent_diff(v))
            .transpose()?;

        let dn = if let Some(dln_real) = dlnfug_dn_real {
            if dln_real.len() != self.nce * self.nce {
                return Err(Error::InvalidInput(
                    "dlnfug_dn_real length must equal nce * nce",
                ));
            }
            let mut out = vec![0.0; self.nc * self.nc];
            let inv_sumn = 1.0 / sumn;
            let inv_sumne = 1.0 / sumne;

            // Upper-left block (solvent-solvent)
            for i in 0..self.ncsym {
                for j in 0..self.ncsym {
                    out[i * self.nc + j] = dln_real[i * self.nce + j] - inv_sumne + inv_sumn;
                }
            }

            // Lower-right block (salt-salt)
            for k in self.ncsym..self.nc {
                for i in self.ncsym..self.nc {
                    let mut val = 0.0;
                    for j in self.ncsym..self.nce {
                        let v_kj = self.v_stoich[k * self.nce + j];
                        if v_kj == 0.0 {
                            continue;
                        }
                        for l in self.ncsym..self.nce {
                            let v_il = self.v_stoich[i * self.nce + l];
                            if v_il == 0.0 {
                                continue;
                            }
                            val += dln_real[l * self.nce + j] * v_kj * v_il;
                        }
                    }
                    let dik = if i == k && n[k] > 0.0 {
                        1.0 / n[k].max(MIN_MOL_NUM)
                    } else {
                        0.0
                    };
                    val += inv_sumn - dik - self.v_sum[k] * self.v_sum[i] * inv_sumne;
                    for j in self.ncsym..self.nce {
                        if ne[j] > 0.0 {
                            let v_kj = self.v_stoich[k * self.nce + j];
                            let v_ij = self.v_stoich[i * self.nce + j];
                            val += v_kj * v_ij / ne[j].max(MIN_MOL_NUM);
                        }
                    }
                    out[k * self.nc + i] = val;
                }
            }

            // Off-diagonal solvent-salt and salt-solvent blocks
            for k in 0..self.ncsym {
                for i in self.ncsym..self.nc {
                    let mut val_ki = 0.0;
                    let mut val_ik = 0.0;
                    for j in self.ncsym..self.nce {
                        let v_ij = self.v_stoich[i * self.nce + j];
                        val_ki += dln_real[k * self.nce + j] * v_ij;
                        val_ik += dln_real[j * self.nce + k] * v_ij;
                    }
                    let shift = inv_sumn - self.v_sum[i] * inv_sumne;
                    out[k * self.nc + i] = val_ki + shift;
                    out[i * self.nc + k] = val_ik + shift;
                }
            }
            Some(out)
        } else {
            None
        };

        Ok(ApparentFugacityResult {
            ln_fugacity,
            dt,
            dp,
            dn,
        })
    }

    /// Back-calculate modified logarithmic fugacity coefficient by removing
    /// explicit `log(x_i)` and `log(xe_k)` ideal-dissociation dependence,
    /// matching `getModFugacity(t, p, X, lnfug, sumne)`.
    ///
    /// Mutates `lnfug` in place and returns `sumne`.
    pub fn mod_fugacity(
        &self,
        temperature: f64,
        pressure: f64,
        x: &[f64],
        lnfug: &mut [f64],
    ) -> Result<f64> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        if x.len() != self.nc || lnfug.len() != self.nc {
            return Err(Error::InvalidInput(
                "composition and lnfug lengths must equal nc",
            ));
        }
        let ne = self.apparent_to_real_mole_numbers(x)?;
        let sumne: f64 = ne.iter().copied().sum();
        positive(sumne, "total real moles")?;
        for i in self.ncsym..self.nc {
            if x[i] <= 0.0 {
                return Err(Error::InvalidInput(
                    "positive apparent salt mole fraction required for mod_fugacity",
                ));
            }
            lnfug[i] += x[i].ln();
            for (k, &nek) in ne.iter().enumerate().skip(self.ncsym) {
                let v_ik = self.v_stoich[i * self.nce + k];
                if v_ik > 0.0 {
                    let xek = nek / sumne;
                    if xek <= 0.0 {
                        return Err(Error::InvalidInput(
                            "positive real ion mole fraction required for mod_fugacity",
                        ));
                    }
                    lnfug[i] -= v_ik * xek.ln();
                    lnfug[i] += v_ik * v_ik.ln();
                }
            }
        }
        Ok(sumne)
    }
}
