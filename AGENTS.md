# Working on CRN-craft

- Read README.md and docs/ai-context.md, then docs/architecture.md and docs/contracts.md. CRN-craft has one engine, the module `crn_craft::rules` of the crate `crn-craft`; release 0.7.0 merged the former network and rules engines into it (docs/campaign/merge-0.7/README.md). For model and method work also read docs/rules-language.md, docs/direct-ssa.md and docs/mtdna-rules.md. Optional workspaces (integrations, the parameter worker) have their own lockfiles and MSRVs; keep them out of the default dependency gate.
- The numerical engine must remain independent of job/campaign infrastructure.
- Preserve strict input handling and stable diagnostic codes. Unknown fields and duplicate JSON keys must never be ignored.
- Generate schemas from Rust types. A delivered contract shape change requires a new version and migration notes.
- Stored `rules-run` requests replay byte for byte (`crates/crn-craft/tests/stored/`), and generated-code fixtures (`crates/crn-craft/tests/generated/`) change only with the generator. Never edit either by hand. New behavior goes behind `crn-craft.model/1.0` or `crn-craft.run/1.0`.
- Match the named milestone. Do not claim a placeholder or planned method can execute.
- For numerical work, declare targets before running experiments. Retain failed rows. Do not loosen a statistical rule after seeing failures. A new sampling policy gets its own sampler identity and an exact-law check (`validation/merge/check_laws.py`).
- Use published algorithms. Retain attribution for any permissively licensed code adaptation; do not copy GPL simulator code into this MIT project.
- Run fmt, clippy with warnings denied, release tests with --locked, rustdoc with warnings denied, ruff, the validation scripts that CI names normally and under Python -O, and the dependency gate for relevant changes.
- Keep fixed-work numerical reproducibility separate from timing, hardware, paths and interrupted runs.
- Include a concise implementation log and validation limitations with each iteration.
- Do not tag, publish packages or push to a remote unless the user requests it.
- Profile models with the generic Static Direct path and its selectors, then the generated event loop where it applies, before adding a specialized kernel. Follow [the rules kernel review checklist](docs/developer/specialized-kernels.md).
