//! A plain RGB image of `f32`, the only type the filters read and write.

/// An RGB image stored row by row, top to bottom. Channels are usually in `0..=1`, but the filters
/// do not assume a range or a colour space: they average whatever they are given.
#[derive(Clone, Debug, PartialEq)]
pub struct Image {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// `width * height` pixels, row by row.
    pub pixels: Vec<[f32; 3]>,
}

impl Image {
    /// An image of `width * height` pixels of `colour`.
    pub fn filled(width: usize, height: usize, colour: [f32; 3]) -> Self {
        Self { width, height, pixels: vec![colour; width * height] }
    }

    /// An image from 8-bit RGB bytes (3 per pixel), mapped to `0..=1`.
    ///
    /// # Panics
    /// If `bytes` does not hold `3 * width * height` values.
    pub fn from_rgb8(width: usize, height: usize, bytes: &[u8]) -> Self {
        assert_eq!(bytes.len(), 3 * width * height, "expected 3 bytes per pixel");
        let pixels = bytes.as_chunks::<3>().0.iter().map(|p| p.map(|v| v as f32 / 255.0)).collect();
        Self { width, height, pixels }
    }

    /// The image as 8-bit RGB bytes, clamped to `0..=1` and rounded.
    pub fn to_rgb8(&self) -> Vec<u8> {
        self.pixels.iter().flat_map(|p| p.map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8)).collect()
    }

    /// The pixel at `(x, y)`, with coordinates clamped to the image (edge pixels repeat outward).
    #[inline]
    pub fn at(&self, x: isize, y: isize) -> [f32; 3] {
        let x = x.clamp(0, self.width as isize - 1) as usize;
        let y = y.clamp(0, self.height as isize - 1) as usize;
        self.pixels[y * self.width + x]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rgb8_round_trips() {
        let bytes: Vec<u8> = (0..2 * 3 * 3).map(|v| (v * 13) as u8).collect();
        assert_eq!(Image::from_rgb8(2, 3, &bytes).to_rgb8(), bytes);
    }

    #[test]
    fn at_clamps_to_the_edge() {
        let mut image = Image::filled(2, 2, [0.0; 3]);
        image.pixels[3] = [1.0; 3];
        assert_eq!(image.at(5, 7), [1.0; 3]);
        assert_eq!(image.at(-3, -1), [0.0; 3]);
    }
}
