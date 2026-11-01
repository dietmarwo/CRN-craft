# Releasing

This M1b delivery is `0.1.0-dev.3`, with `publish = false` for both crates. It is a reviewable source project, not a simulator release. Publish a simulator release only after its named numerical acceptance gate passes.

For a future release: complete the declared milestone, regenerate/version schemas, run the documented checks on supported toolchains/platforms, review dependency licenses and reference-data provenance, update changelog/capabilities, and run installed-package acceptance. Preserve lockfiles and evidence. Package availability must be checked at release time; earlier planning searches do not reserve names.

The owner decides when to push, tag and publish. No release credentials or workflow that automatically publishes packages are included. MIT project code does not relicense dependencies or reference datasets.
