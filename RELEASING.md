# Releasing

This M0 delivery is `0.1.0-dev.1`, with `publish = false` for both crates. It is a reviewable source project, not a simulator release. Do not publish a crate or Python wheel that implies implemented numerical features.

For a future release: complete the declared milestone, regenerate/version schemas, run the documented checks on supported toolchains/platforms, review dependency licenses and reference-data provenance, update changelog/capabilities, and run installed-package acceptance. Preserve lockfiles and evidence. Package availability must be checked at release time; earlier planning searches do not reserve names.

The owner decides when to create a repository, push, tag and publish. No release credentials or workflow that automatically publishes packages are included. MIT project code does not relicense dependencies or reference datasets.
