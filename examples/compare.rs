//! Writes `images/compare.png` for the README: a detail of a Kodak test photograph, then the
//! classic, generalized and anisotropic filters, and the anisotropic filter applied three times.
//!
//!     cargo run --release --example compare

use kuwahara::{Image, Params, anisotropic, classic};

/// The detail shown: x, y, width, height in `kodim23.png`.
const CROP: (u32, u32, u32, u32) = (380, 60, 300, 220);

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rgb = ::image::open("images/input/kodim23.png")?.to_rgb8();
    let (w, h) = (rgb.width() as usize, rgb.height() as usize);
    let input = Image::from_rgb8(w, h, rgb.as_raw());
    let isotropic = Params { alpha: f32::INFINITY, ..Params::default() };
    let once = anisotropic(&input, &Params::default());
    let thrice = (0..2).fold(once.clone(), |p, _| anisotropic(&p, &Params::default()));
    let panels = [input.clone(), classic(&input, 6), anisotropic(&input, &isotropic), once, thrice];

    let (cx, cy, cw, ch) = CROP;
    let gap = 4;
    let sheet_w = panels.len() as u32 * (cw + gap) - gap;
    let mut sheet = ::image::RgbImage::from_pixel(sheet_w, ch, ::image::Rgb([255, 255, 255]));
    for (k, panel) in panels.iter().enumerate() {
        let full = ::image::RgbImage::from_raw(w as u32, h as u32, panel.to_rgb8()).ok_or("bad buffer")?;
        let detail = ::image::imageops::crop_imm(&full, cx, cy, cw, ch).to_image();
        ::image::imageops::replace(&mut sheet, &detail, (k as u32 * (cw + gap)) as i64, 0);
    }
    sheet.save("images/compare.png")?;
    println!("wrote images/compare.png");
    Ok(())
}
