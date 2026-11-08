//! The anisotropic Kuwahara filter of Kyprianidis, Kang and Döllner, *Image and Video Abstraction
//! by Anisotropic Kuwahara Filtering* (Pacific Graphics 2009), with the polynomial weighting
//! functions of Kyprianidis, Semmo, Kang and Döllner, *Anisotropic Kuwahara Filtering with
//! Polynomial Weighting Functions* (EG UK TPCG 2010) as an alternative.
//!
//! Around each pixel the filter lays an ellipse whose major axis follows the local
//! orientation ([`orientation`]): a circle of radius `r` where the image is isotropic,
//! stretched up to `2r` by `r/2` along edges. The ellipse is cut into `N` overlapping sectors by
//! smooth weighting functions, defined on the unit disc and pulled back to the ellipse. Each
//! sector yields a weighted mean `mᵢ` and standard deviation `sᵢ`; the output is the mean of the
//! `mᵢ` weighted by `αᵢ = 1 / (1 + ‖sᵢ‖^q)`, so uniform sectors (those on one side of an edge)
//! dominate.

use std::f32::consts::{PI, TAU};

use crate::Image;
use crate::tensor::orientation;

/// How the ellipse is divided into sectors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Weighting {
    /// The 2009 paper (and Papari et al.): the characteristic function of each sector convolved
    /// with a Gaussian of deviation `σ_r / 3`, times a Gaussian of deviation `σ_r = 0.4` (on the
    /// unit disc). Sampled into a lookup table once per filter call, as the paper's texture map.
    Convolved,
    /// The 2010 paper: `[(x + ζ) − η y²]²₊` per sector, normalised to sum to the same Gaussian.
    /// No table and no convolution. `ζ` sets how much the sectors overlap at the centre (the
    /// paper's `2 / r` when `None`), `η` how much they overlap at their sides (when `None`, from the
    /// paper's choice that a sector's parabola meets the unit circle at `±3π / 2N`).
    Polynomial {
        /// Overlap at the centre.
        zeta: Option<f32>,
        /// Overlap at the sides.
        eta: Option<f32>,
    },
}

/// Parameters of [`anisotropic()`]. The defaults are the papers'.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Params {
    /// Radius `r` of the filter in isotropic regions, in pixels. The ellipse's axes range from
    /// `r/2` to `2r` (for `alpha = 1`). Default 6.
    pub radius: f32,
    /// Sharpness `q` of the choice between sectors: 0 averages them, large values pick the most
    /// uniform. Default 8.
    pub q: f32,
    /// Eccentricity control `α > 0`: the ellipse's axes are `r(α + A)/α` and `rα/(α + A)`.
    /// Infinity keeps a circle everywhere: the generalized Kuwahara filter. Default 1.
    pub alpha: f32,
    /// Standard deviation of the Gaussian that smooths the structure tensor. Default 2.
    pub sigma: f32,
    /// Number of sectors `N`. Default 8.
    pub sectors: usize,
    /// The papers define `αᵢ` for channels in `0..=255`; the deviations are multiplied by this
    /// before the power. 255 for images in `0..=1`. Default 255.
    pub scale: f32,
    /// The sector weighting functions. Default [`Weighting::Polynomial`] with the paper's values.
    pub weighting: Weighting,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            radius: 6.0,
            q: 8.0,
            alpha: 1.0,
            sigma: 2.0,
            sectors: 8,
            scale: 255.0,
            weighting: Weighting::Polynomial { zeta: None, eta: None },
        }
    }
}

/// Standard deviation of the radial Gaussian on the unit disc (both papers).
const SIGMA_R: f32 = 0.4;

/// The anisotropic Kuwahara filter. Rows are split across the available threads.
pub fn anisotropic(image: &Image, params: &Params) -> Image {
    let field = orientation(image, params.sigma);
    let sectors = Sectors::new(params);
    let (w, h) = (image.width, image.height);
    let mut out = Image::filled(w, h, [0.0; 3]);
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(h.max(1));
    let rows_per = h.div_ceil(threads.max(1)).max(1);
    std::thread::scope(|scope| {
        for (chunk, rows) in out.pixels.chunks_mut(rows_per * w).enumerate() {
            let (field, sectors) = (&field, &sectors);
            scope.spawn(move || {
                let mut scratch = Scratch::new(params.sectors);
                for (i, pixel) in rows.iter_mut().enumerate() {
                    let c = chunk * rows_per * w + i;
                    let o = field[c];
                    *pixel = filter_pixel(
                        image,
                        (c % w) as isize,
                        (c / w) as isize,
                        o.angle,
                        o.anisotropy,
                        params,
                        sectors,
                        &mut scratch,
                    );
                }
            });
        }
    });
    out
}

/// Per-thread accumulators: weighted sums of colour, squared colour and weight per sector.
struct Scratch {
    m: Vec<[f32; 3]>,
    s: Vec<[f32; 3]>,
    k: Vec<f32>,
    weights: Vec<f32>,
}

impl Scratch {
    fn new(n: usize) -> Self {
        Self { m: vec![[0.0; 3]; n], s: vec![[0.0; 3]; n], k: vec![0.0; n], weights: vec![0.0; n] }
    }
}

#[allow(clippy::too_many_arguments)]
fn filter_pixel(
    image: &Image,
    x0: isize,
    y0: isize,
    angle: f32,
    anisotropy: f32,
    params: &Params,
    sectors: &Sectors,
    scratch: &mut Scratch,
) -> [f32; 3] {
    let r = params.radius;
    // Semi-axes: a along the orientation, b across it (2010, §2).
    let (a, b) = if params.alpha.is_finite() {
        let e = (params.alpha + anisotropy) / params.alpha;
        (r * e, r / e)
    } else {
        (r, r)
    };
    let (sin, cos) = angle.sin_cos();
    // Half-extents of the ellipse's bounding box.
    let ex = (a * a * cos * cos + b * b * sin * sin).sqrt().ceil() as isize;
    let ey = (a * a * sin * sin + b * b * cos * cos).sqrt().ceil() as isize;

    let Scratch { m, s, k, weights } = scratch;
    m.fill([0.0; 3]);
    s.fill([0.0; 3]);
    k.fill(0.0);
    for dy in -ey..=ey {
        for dx in -ex..=ex {
            // v = S R₋φ d: the ellipse mapped onto the unit disc.
            let (fx, fy) = (dx as f32, dy as f32);
            let v = [(cos * fx + sin * fy) / a, (-sin * fx + cos * fy) / b];
            if v[0] * v[0] + v[1] * v[1] > 1.0 {
                continue;
            }
            sectors.weights(v, weights);
            let c = image.at(x0 + dx, y0 + dy);
            for i in 0..weights.len() {
                let wi = weights[i];
                if wi == 0.0 {
                    continue;
                }
                for ch in 0..3 {
                    m[i][ch] += wi * c[ch];
                    s[i][ch] += wi * c[ch] * c[ch];
                }
                k[i] += wi;
            }
        }
    }

    let (mut out, mut total) = ([0.0f32; 3], 0.0f32);
    for i in 0..m.len() {
        if k[i] <= 0.0 {
            continue;
        }
        let mean = m[i].map(|v| v / k[i]);
        let variance: f32 = (0..3).map(|ch| (s[i][ch] / k[i] - mean[ch] * mean[ch]).abs()).sum();
        let deviation = variance.sqrt() * params.scale;
        let alpha = 1.0 / (1.0 + deviation.powf(params.q));
        for ch in 0..3 {
            out[ch] += alpha * mean[ch];
        }
        total += alpha;
    }
    if total > 0.0 { out.map(|v| v / total) } else { image.at(x0, y0) }
}

/// The weighting functions `Kᵢ` on the unit disc.
enum Sectors {
    /// `K₀` sampled on a grid over `[-1, 1]²`; `Kᵢ = K₀ ∘ R₋₂πi/N`.
    Table { side: usize, k0: Vec<f32>, rotations: Vec<(f32, f32)> },
    /// The polynomial weights; per sector the unit vector it points along.
    Polynomial { zeta: f32, eta: f32, axes: Vec<(f32, f32)> },
}

impl Sectors {
    fn new(params: &Params) -> Self {
        let n = params.sectors.max(1);
        match params.weighting {
            Weighting::Convolved => {
                let side = 129;
                let rotations = (0..n).map(|i| (-TAU * i as f32 / n as f32).sin_cos()).collect();
                Sectors::Table { side, k0: convolved_k0(side, n), rotations }
            }
            Weighting::Polynomial { zeta, eta } => {
                let zeta = zeta.unwrap_or(2.0 / params.radius);
                let gamma = 3.0 * PI / (2.0 * n as f32);
                let eta = eta.unwrap_or((zeta + gamma.cos()) / (gamma.sin() * gamma.sin()));
                let axes = (0..n)
                    .map(|i| {
                        let (s, c) = (TAU * i as f32 / n as f32).sin_cos();
                        (c, s)
                    })
                    .collect();
                Sectors::Polynomial { zeta, eta, axes }
            }
        }
    }

    /// The weights of every sector at `v` in the unit disc.
    #[inline]
    fn weights(&self, v: [f32; 2], out: &mut [f32]) {
        match self {
            Sectors::Table { side, k0, rotations } => {
                for (o, &(s, c)) in out.iter_mut().zip(rotations) {
                    let u = [c * v[0] - s * v[1], s * v[0] + c * v[1]];
                    *o = bilinear(k0, *side, u);
                }
            }
            Sectors::Polynomial { zeta, eta, axes } => {
                let mut sum = 0.0;
                for (o, &(c, s)) in out.iter_mut().zip(axes) {
                    // Coordinates along the sector's axis and across it.
                    let along = c * v[0] + s * v[1];
                    let across = -s * v[0] + c * v[1];
                    let z = (along + zeta - eta * across * across).max(0.0);
                    *o = z * z;
                    sum += *o;
                }
                let g = (-(v[0] * v[0] + v[1] * v[1]) / (2.0 * SIGMA_R * SIGMA_R)).exp();
                let g = if sum > 0.0 { g / sum } else { 0.0 };
                out.iter_mut().for_each(|o| *o *= g);
            }
        }
    }
}

/// `K₀ = (χ₀ ⋆ G_σs) · G_σr` sampled on `side × side` points over `[-1, 1]²`, where `χ₀` is the
/// sector of arguments in `(−π/N, π/N]` and `σ_s = σ_r / 3`.
fn convolved_k0(side: usize, n: usize) -> Vec<f32> {
    let cell = 2.0 / (side - 1) as f32;
    let sigma_s = SIGMA_R / 3.0;
    let margin = (3.0 * sigma_s / cell).ceil() as usize;
    // χ₀ on a grid with a margin wide enough for the convolution.
    let big = side + 2 * margin;
    let half = PI / n as f32;
    let at = |i: usize| (i as f32 - margin as f32) * cell - 1.0;
    let mut chi: Vec<[f32; 1]> = Vec::with_capacity(big * big);
    for j in 0..big {
        for i in 0..big {
            let (x, y) = (at(i), at(j));
            // The fraction of the cell inside the sector, from 4 × 4 samples, so that the
            // sectors still add up to 1 at the centre and along their borders.
            let mut inside = 0;
            for sj in 0..4 {
                for si in 0..4 {
                    let (sx, sy) = (x + (si as f32 - 1.5) * cell / 4.0, y + (sj as f32 - 1.5) * cell / 4.0);
                    let arg = sy.atan2(sx);
                    inside += (arg > -half && arg <= half) as u32;
                }
            }
            chi.push([inside as f32 / 16.0]);
        }
    }
    let smooth = crate::tensor::gaussian(&chi, big, big, sigma_s / cell);
    let mut k0 = Vec::with_capacity(side * side);
    for j in 0..side {
        for i in 0..side {
            let (x, y) = (at(i + margin), at(j + margin));
            let g = (-(x * x + y * y) / (2.0 * SIGMA_R * SIGMA_R)).exp();
            k0.push(smooth[(j + margin) * big + i + margin][0] * g);
        }
    }
    k0
}

/// Bilinear lookup in a `side × side` table over `[-1, 1]²`.
#[inline]
fn bilinear(table: &[f32], side: usize, u: [f32; 2]) -> f32 {
    let scale = (side - 1) as f32 / 2.0;
    let (x, y) =
        (((u[0] + 1.0) * scale).clamp(0.0, (side - 1) as f32), ((u[1] + 1.0) * scale).clamp(0.0, (side - 1) as f32));
    let (i, j) = ((x as usize).min(side - 2), (y as usize).min(side - 2));
    let (fx, fy) = (x - i as f32, y - j as f32);
    let t = |i: usize, j: usize| table[j * side + i];
    let top = t(i, j) * (1.0 - fx) + t(i + 1, j) * fx;
    let bottom = t(i, j + 1) * (1.0 - fx) + t(i + 1, j + 1) * fx;
    top * (1.0 - fy) + bottom * fy
}

#[cfg(test)]
mod tests {
    use super::*;

    fn half_and_half() -> Image {
        let mut image = Image::filled(24, 24, [0.1, 0.2, 0.3]);
        for y in 0..24 {
            for x in 12..24 {
                image.pixels[y * 24 + x] = [0.9, 0.8, 0.7];
            }
        }
        image
    }

    #[test]
    fn keeps_a_flat_image() {
        for weighting in [Weighting::Convolved, Weighting::Polynomial { zeta: None, eta: None }] {
            let image = Image::filled(10, 10, [0.25, 0.5, 0.75]);
            let out = anisotropic(&image, &Params { weighting, ..Params::default() });
            for p in &out.pixels {
                for (v, expected) in p.iter().zip(image.pixels[0]) {
                    assert!((v - expected).abs() < 1e-5);
                }
            }
        }
    }

    #[test]
    fn keeps_an_edge_sharp() {
        for weighting in [Weighting::Convolved, Weighting::Polynomial { zeta: None, eta: None }] {
            let image = half_and_half();
            let out = anisotropic(&image, &Params { weighting, ..Params::default() });
            // The pixels either side of the edge stay close to their own side's colour.
            let (left, right) = (out.pixels[12 * 24 + 11], out.pixels[12 * 24 + 12]);
            assert!((left[0] - 0.1).abs() < 0.05, "{weighting:?} {left:?}");
            assert!((right[0] - 0.9).abs() < 0.05, "{weighting:?} {right:?}");
        }
    }

    #[test]
    fn sector_weights_partition_the_gaussian() {
        for weighting in [Weighting::Convolved, Weighting::Polynomial { zeta: None, eta: None }] {
            let sectors = Sectors::new(&Params { weighting, ..Params::default() });
            let mut w = vec![0.0; 8];
            for v in [[0.0, 0.0], [0.3, -0.2], [-0.5, 0.6], [0.1, 0.9]] {
                sectors.weights(v, &mut w);
                let g = (-(v[0] * v[0] + v[1] * v[1]) / (2.0 * SIGMA_R * SIGMA_R)).exp();
                let sum: f32 = w.iter().sum();
                assert!((sum - g).abs() < 0.02 * g.max(0.05), "{weighting:?} at {v:?}: {sum} vs {g}");
            }
        }
    }

    #[test]
    fn polynomial_eta_matches_the_paper() {
        let Sectors::Polynomial { zeta, eta, .. } = Sectors::new(&Params::default()) else { unreachable!() };
        assert!((zeta - 1.0 / 3.0).abs() < 1e-6);
        assert!((eta - 3.77).abs() < 0.01, "{eta}");
    }
}
