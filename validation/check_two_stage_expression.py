"""Check the native FSP and the independent generating function."""

import sys
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

sys.path.insert(
    0, str(Path(__file__).resolve().parents[1] / "examples" / "two_stage_expression")
)
from model import joint_generator, protein_pmf  # noqa: E402
from study import native_check  # noqa: E402


def main():
    check = native_check()
    if max(check.values()) > 1e-9:
        raise SystemExit("Independent and native FSP disagree")
    law = joint_generator(0.7, 1.2, 2.0, 12, 30)
    init = np.zeros(law.shape[0])
    init[0] = 1.0
    p = expm_multiply(law, init)
    approximate = protein_pmf(0.7, 1.2, 2.0, 1.0, 512)
    gap = np.max(np.abs(p[:-1].reshape(13, 31).sum(axis=0) - approximate[:31]))
    if gap > 1e-9:
        raise SystemExit("Protein PMF disagrees with joint count law")
    print("native FSP, matrix exponential and generating function verified")


if __name__ == "__main__":
    main()
