//! The Ziggurat exponential sampler of the unified engine.
//!
//! [`exponential_from_words`] draws a unit exponential from the caller's
//! random words (Marsaglia & Tsang's Ziggurat, adapted from Julia Random; see
//! THIRD_PARTY_NOTICES.md). Static Direct, next reaction, the mtDNA path and
//! generated loops use it for `ziggurat` waits. Release 0.7.0 removed the
//! network engine's alternative exact SSA methods (`crn-craft.ssa-run`,
//! `crn-craft.ssa-result`) that this module held; their selectors live on in
//! [`crate::rules`] as the Static Direct selectors.
mod ziggurat;
pub use ziggurat::exponential_from_words;
mod ziggurat_tables;
