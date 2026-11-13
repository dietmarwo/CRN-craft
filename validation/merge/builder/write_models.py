"""Write the builder-form models of the E1 law check (validation/merge/RESULTS.md).

The models are `three_class()` of guide chapter 12, as
`validation/merge/manual-ratios/builder_form.py` defines it, at a capacity
small enough for a finite state projection. Run with any `crn_craft` wheel:

    python validation/merge/builder/write_models.py
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "manual-ratios"))

from builder_form import three_class  # noqa: E402

MODELS = {
    "mtdna-builder": {"takeover": None},
    "mtdna-builder-stop": {"takeover": 0.5},
}


def main():
    for name, options in MODELS.items():
        model = three_class(K=20, w0=20, mu=0.05, half_life_days=1, **options)
        (HERE / f"{name}.model.json").write_text(model.to_json() + "\n")
        print(name)


if __name__ == "__main__":
    main()
