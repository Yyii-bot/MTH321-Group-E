"""
run_all.py - one-command reproduction of all figures and console output
(Brief Section 7 acceptance standard).

On a machine with numpy/scipy/matplotlib installed, run
    python run_all.py
to regenerate all seven figures figures/fig01..fig07 from a clean state and to
print every validation experiment's complete console output (every number in
the report comes from this output).

Author: Rong Yin (Role 5: Testing & Validation)
"""

import sys
from pathlib import Path

# make the sibling modules importable regardless of the launch directory
sys.path.insert(0, str(Path(__file__).resolve().parent))

import reference
import validation_tests
import edge_cases
import convergence_study
import final_check


def main():
    steps = [
        ("Week 1  reference-solution three-level chain",
         reference.reference_report),
        ("Week 2  unit-test suite (PASS/FAIL)", validation_tests.main),
        ("Week 2-3  edge-case hunting experiments", edge_cases.main),
        ("Week 3  convergence-order verification",
         convergence_study.main),
        ("Week 4  pre-submission final check", final_check.main),
    ]
    n_ok = 0
    for title, fn in steps:
        print()
        print("#" * 72)
        print(f"# {title}")
        print("#" * 72)
        try:
            fn()
            n_ok += 1
        except Exception as exc:
            print(f"[run_all] this step raised an exception: {exc!r}")
    print()
    print("#" * 72)
    print(f"# run_all finished: {n_ok}/{len(steps)} steps succeeded")
    print("#" * 72)
    return n_ok == len(steps)


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
