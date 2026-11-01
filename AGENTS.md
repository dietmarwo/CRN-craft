# Working on CRN-craft

- Read README.md and docs/ai-context.md. The current implementation is M2: Direct SSA, ordered ensembles, budgets/cancellation and verifiable archives, alongside the reviewed M1 analysis/exporters. Read docs/direct-ssa.md and docs/archive-contract.md.
- The numerical engine must remain independent of job/campaign infrastructure.
- Preserve strict input handling and stable diagnostic codes. Unknown fields and duplicate JSON keys must never be ignored.
- Generate schemas from Rust types. A delivered contract shape change requires a new version and migration notes.
- Match the named milestone. Do not claim a placeholder or planned method can execute.
- For numerical work, declare targets before running experiments. Retain failed rows. Do not loosen a statistical rule after seeing failures.
- Use published algorithms. Retain attribution for any permissively licensed code adaptation; do not copy GPL simulator code into this MIT project.
- Run fmt, clippy with warnings denied, release tests with --locked, the schema/CLI checker normally and under Python -O, and the dependency gate for relevant changes.
- Keep fixed-work numerical reproducibility separate from timing, hardware, paths and interrupted runs.
- Include a concise implementation log and validation limitations with each iteration.
- Do not tag, publish packages or push to a remote unless the user requests it.
