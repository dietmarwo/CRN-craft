# Releasing

This M2 delivery is `0.1.0-dev.4`, with `publish = false` for both crates. Direct SSA is implemented and locally validated; this remains a reviewable development source delivery. Public package release still requires owner approval and installed-package/platform acceptance.

For a future release: complete the declared milestone, regenerate/version schemas, run the documented checks on supported toolchains/platforms, review dependency licenses and reference-data provenance, update changelog/capabilities, and run installed-package acceptance. Preserve lockfiles and evidence. Package availability must be checked at release time; earlier planning searches do not reserve names.

The owner decides when to push, tag and publish. No release credentials or workflow that automatically publishes packages are included. MIT project code does not relicense dependencies or reference datasets.
