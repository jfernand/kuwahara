//! Kuwahara filters: edge-preserving smoothing that turns photographs into flat, painterly patches.
//!
//! - [`classic()`]: Kuwahara et al. (1976). Four square regions per pixel; the least varied wins.
//! - [`anisotropic()`]: Kyprianidis, Kang and Döllner, *Image and Video Abstraction by Anisotropic
//!   Kuwahara Filtering* (Pacific Graphics 2009). Eight smoothly weighted sectors of an ellipse
//!   that follows the local orientation of the image ([`orientation()`]), blended by how uniform
//!   each sector is. With the anisotropy switched off it is the generalized Kuwahara filter of
//!   Papari, Petkov and Campisi (2007).
//!
//! The filters work on an [`Image`] of `f32` RGB and have no dependencies.

#![warn(missing_docs)]

pub mod anisotropic;
pub mod classic;
pub mod image;
pub mod tensor;

pub use anisotropic::{Params, Weighting, anisotropic};
pub use classic::classic;
pub use image::Image;
pub use tensor::{Orientation, orientation};

#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
