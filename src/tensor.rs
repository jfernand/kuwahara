//! Local orientation and anisotropy from the smoothed structure tensor (Kyprianidis, Kang and
//! Döllner 2009, §3; after Brox et al. and Weickert).
//!
//! The tensor sums the outer products of every channel's gradient, then a Gaussian averages it
//! over a neighbourhood. Its eigenvector of least change follows edges and stripes; the gap
//! between its eigenvalues says how strongly oriented the neighbourhood is.

use crate::Image;

/// The orientation of the neighbourhood of a pixel.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Orientation {
    /// Angle of the direction of least change (along edges), in radians, `x` right and `y` down.
    pub angle: f32,
    /// `(λ1 − λ2) / (λ1 + λ2)`: 0 where the neighbourhood is flat or isotropic, approaching 1
    /// along a clean edge.
    pub anisotropy: f32,
}

/// The orientation field of `image`: Sobel gradients, structure tensor, a Gaussian of standard
/// deviation `sigma` over the tensor, then its eigen-analysis.
pub fn orientation(image: &Image, sigma: f32) -> Vec<Orientation> {
    let (w, h) = (image.width, image.height);
    // The tensor (E, F, G) = Σ (fx·fx, fx·fy, fy·fy) over channels.
    let mut tensor = vec![[0.0f32; 3]; w * h];
    for y in 0..h as isize {
        for x in 0..w as isize {
            let p = |dx: isize, dy: isize| image.at(x + dx, y + dy);
            let (mut e, mut f, mut g) = (0.0, 0.0, 0.0);
            for k in 0..3 {
                let fx = (p(1, -1)[k] + 2.0 * p(1, 0)[k] + p(1, 1)[k] - p(-1, -1)[k] - 2.0 * p(-1, 0)[k] - p(-1, 1)[k])
                    / 4.0;
                let fy = (p(-1, 1)[k] + 2.0 * p(0, 1)[k] + p(1, 1)[k] - p(-1, -1)[k] - 2.0 * p(0, -1)[k] - p(1, -1)[k])
                    / 4.0;
                e += fx * fx;
                f += fx * fy;
                g += fy * fy;
            }
            tensor[y as usize * w + x as usize] = [e, f, g];
        }
    }
    let tensor = gaussian(&tensor, w, h, sigma);
    tensor
        .iter()
        .map(|&[e, f, g]| {
            let root = ((e - g) * (e - g) + 4.0 * f * f).sqrt();
            let (l1, l2) = ((e + g + root) / 2.0, (e + g - root) / 2.0);
            // Eigenvector of λ1 is the gradient direction; its normal (λ1 − E, −F) runs along the edge.
            let (tx, ty) = (l1 - e, -f);
            let angle = if tx == 0.0 && ty == 0.0 { 0.0 } else { ty.atan2(tx) };
            let anisotropy = if l1 + l2 > 0.0 { (l1 - l2) / (l1 + l2) } else { 0.0 };
            Orientation { angle, anisotropy }
        })
        .collect()
}

/// A separable Gaussian blur of standard deviation `sigma`, edges clamped.
pub(crate) fn gaussian<const N: usize>(data: &[[f32; N]], w: usize, h: usize, sigma: f32) -> Vec<[f32; N]> {
    if sigma <= 0.0 {
        return data.to_vec();
    }
    let radius = (3.0 * sigma).ceil() as isize;
    let kernel: Vec<f32> = (-radius..=radius).map(|i| (-(i * i) as f32 / (2.0 * sigma * sigma)).exp()).collect();
    let total: f32 = kernel.iter().sum();
    let kernel: Vec<f32> = kernel.iter().map(|k| k / total).collect();
    let pass = |src: &[[f32; N]], horizontal: bool| -> Vec<[f32; N]> {
        let mut out = vec![[0.0; N]; w * h];
        for y in 0..h as isize {
            for x in 0..w as isize {
                let mut acc = [0.0; N];
                for (k, &weight) in kernel.iter().enumerate() {
                    let d = k as isize - radius;
                    let (sx, sy) = if horizontal { (x + d, y) } else { (x, y + d) };
                    let (sx, sy) = (sx.clamp(0, w as isize - 1) as usize, sy.clamp(0, h as isize - 1) as usize);
                    let v = src[sy * w + sx];
                    for c in 0..N {
                        acc[c] += weight * v[c];
                    }
                }
                out[y as usize * w + x as usize] = acc;
            }
        }
        out
    };
    pass(&pass(data, true), false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::FRAC_PI_2;

    fn stripes(vertical: bool) -> Image {
        let mut image = Image::filled(32, 32, [0.0; 3]);
        for y in 0..32 {
            for x in 0..32 {
                let t = if vertical { x } else { y };
                image.pixels[y * 32 + x] = [((t as f32) * 0.7).sin() * 0.5 + 0.5; 3];
            }
        }
        image
    }

    #[test]
    fn follows_vertical_stripes() {
        let o = orientation(&stripes(true), 2.0)[16 * 32 + 16];
        assert!((o.angle.abs() - FRAC_PI_2).abs() < 0.05, "{o:?}");
        assert!(o.anisotropy > 0.9, "{o:?}");
    }

    #[test]
    fn follows_horizontal_stripes() {
        let o = orientation(&stripes(false), 2.0)[16 * 32 + 16];
        assert!(o.angle.sin().abs() < 0.05, "{o:?}");
        assert!(o.anisotropy > 0.9, "{o:?}");
    }

    #[test]
    fn flat_is_isotropic() {
        let o = orientation(&Image::filled(8, 8, [0.3; 3]), 2.0)[0];
        assert_eq!(o.anisotropy, 0.0);
    }
}
