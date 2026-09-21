"""Portable repository entry point: python code/run_all.py --help."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from robertson.experiments import main

if __name__ == "__main__":
    raise SystemExit(main())
