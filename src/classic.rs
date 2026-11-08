//! The original filter of Kuwahara, Hachimura, Eiho and Kinoshita, *Processing of
//! RI-angiocardiographic images* (1976).
//!
//! Around each pixel lie four square regions of side `radius + 1`, one per quadrant, each
//! containing the pixel itself. The output is the mean colour of the region with the smallest
//! variance. Edges survive because a region straddling an edge has high variance and loses to one
//! lying on a single side of it; flat areas are smoothed. The squares leave blocky artefacts and
//! the hard choice between regions is unstable in noise, which the generalized and anisotropic
//! filters address.

use crate::Image;

/// The classic Kuwahara filter with square regions of side `radius + 1`.
///
/// A region's variance is the sum of its three channels' variances.
pub fn classic(image: &Image, radius: usize) -> Image {
    let r = radius as isize;
    let mut out = Image::filled(image.width, image.height, [0.0; 3]);
    for y in 0..image.height as isize {
        for x in 0..image.width as isize {
            let mut best = (f32::INFINITY, [0.0; 3]);
            for (sx, sy) in [(-1, -1), (1, -1), (-1, 1), (1, 1)] {
                let (mut m, mut s) = ([0.0f32; 3], [0.0f32; 3]);
                for j in 0..=r {
                    for i in 0..=r {
                        let c = image.at(x + sx * i, y + sy * j);
                        for k in 0..3 {
                            m[k] += c[k];
                            s[k] += c[k] * c[k];
                        }
                    }
                }
                let n = ((r + 1) * (r + 1)) as f32;
                let mean = m.map(|v| v / n);
                let variance: f32 = (0..3).map(|k| s[k] / n - mean[k] * mean[k]).sum();
                if variance < best.0 {
                    best = (variance, mean);
                }
            }
            out.pixels[y as usize * image.width + x as usize] = best.1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_a_sharp_edge() {
        // Left half black, right half white: every pixel has a region on its own side.
        let mut image = Image::filled(16, 8, [0.0; 3]);
        for y in 0..8 {
            for x in 8..16 {
                image.pixels[y * 16 + x] = [1.0; 3];
            }
        }
        assert_eq!(classic(&image, 3), image);
    }

    #[test]
    fn flattens_noise() {
        let mut image = Image::filled(9, 9, [0.5; 3]);
        image.pixels[4 * 9 + 4] = [1.0; 3];
        let out = classic(&image, 2);
        assert!(out.pixels[4 * 9 + 4][0] < 0.6);
    }
}
