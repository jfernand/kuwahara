//! Filters an image file and writes the result.
//!
//!     cargo run --release --example filter -- INPUT OUTPUT [--classic] [--radius R] [--q Q]
//!         [--alpha A] [--sigma S] [--isotropic] [--convolved] [--sectors N] [--passes N]
//!         [--orientation]
//!
//! By default it applies the anisotropic filter with the paper's parameters. `--classic` uses
//! Kuwahara's square regions instead, `--isotropic` the generalized (circular) filter, `--convolved` the 2009
//! paper's sector weights instead of the 2010 polynomial ones, and
//! `--orientation` writes the orientation field instead of filtering (hue = angle, value =
//! anisotropy). PNG and JPEG are read and written.

use std::path::PathBuf;
use std::time::Instant;

use kuwahara::{Image, Params, Weighting, anisotropic, classic, orientation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let mut paths = Vec::new();
    let mut params = Params::default();
    let (mut use_classic, mut show_orientation, mut passes) = (false, false, 1);
    while let Some(arg) = args.next() {
        let mut value = || -> Result<f32, Box<dyn std::error::Error>> {
            Ok(args.next().ok_or(format!("{arg} needs a value"))?.parse()?)
        };
        match arg.as_str() {
            "--classic" => use_classic = true,
            "--isotropic" => params.alpha = f32::INFINITY,
            "--convolved" => params.weighting = Weighting::Convolved,
            "--sectors" => params.sectors = value()? as usize,
            "--orientation" => show_orientation = true,
            "--radius" => params.radius = value()?,
            "--q" => params.q = value()?,
            "--alpha" => params.alpha = value()?,
            "--sigma" => params.sigma = value()?,
            "--passes" => passes = value()? as usize,
            _ if arg.starts_with("--") => return Err(format!("unknown option {arg}").into()),
            _ => paths.push(PathBuf::from(arg)),
        }
    }
    let [input, output] = <[PathBuf; 2]>::try_from(paths).map_err(|_| "usage: filter INPUT OUTPUT [options]")?;

    let rgb = ::image::open(&input)?.to_rgb8();
    let (w, h) = (rgb.width() as usize, rgb.height() as usize);
    let mut picture = Image::from_rgb8(w, h, rgb.as_raw());
    let start = Instant::now();
    if show_orientation {
        picture = orientation_image(&picture, params.sigma);
    } else {
        for _ in 0..passes {
            picture =
                if use_classic { classic(&picture, params.radius as usize) } else { anisotropic(&picture, &params) };
        }
    }
    let elapsed = start.elapsed();
    ::image::RgbImage::from_raw(w as u32, h as u32, picture.to_rgb8()).ok_or("bad buffer")?.save(&output)?;
    println!("{} -> {} ({w}x{h}, {:.2?})", input.display(), output.display(), elapsed);
    Ok(())
}

/// The orientation field as colour: hue is the angle (mod π), value the anisotropy.
fn orientation_image(picture: &Image, sigma: f32) -> Image {
    let field = orientation(picture, sigma);
    let pixels = field
        .iter()
        .map(|o| {
            let hue = (o.angle.rem_euclid(std::f32::consts::PI) / std::f32::consts::PI) * 6.0;
            let v = o.anisotropy;
            let x = 1.0 - (hue % 2.0 - 1.0).abs();
            let [r, g, b] = match hue as u32 {
                0 => [1.0, x, 0.0],
                1 => [x, 1.0, 0.0],
                2 => [0.0, 1.0, x],
                3 => [0.0, x, 1.0],
                4 => [x, 0.0, 1.0],
                _ => [1.0, 0.0, x],
            };
            [r * v, g * v, b * v]
        })
        .collect();
    Image { width: picture.width, height: picture.height, pixels }
}
