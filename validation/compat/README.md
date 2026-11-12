# Compatibility copy of the 0.6.0 manual

`manual-0.6.0.md` is `docs/mtdna-user-manual.md` exactly as it was at the
release tag `v0.6.0` (`git show v0.6.0:docs/mtdna-user-manual.md`). Never
edit it. CI runs every block with `python validation/check_manual.py --manual
validation/compat/manual-0.6.0.md`, normally and under `-O`, and compares the
printed output: code written for 0.6.0 against `crn_craft.rules` must keep
running with the same results.
