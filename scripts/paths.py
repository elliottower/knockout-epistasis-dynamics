"""The only module that computes a location from __file__. Everything else imports a name."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT_ROOT / "results"
AUDIT = RESULTS / "audit"
REFERENCE = PROJECT_ROOT / "reference"
PAPER = PROJECT_ROOT / "paper"
