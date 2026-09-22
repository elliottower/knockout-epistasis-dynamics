"""A stand-in that shadows the project's data_utils: it defines the same names, loaded from the
real file under $SHADOWED_ROOT, but is itself imported from here. Placed earlier on the import
path, it lets every import succeed, and the provenance check must still refuse to run."""

import importlib.util
import os
from pathlib import Path

_spec = importlib.util.spec_from_file_location("_shadowed_data_utils", Path(os.environ["SHADOWED_ROOT"]) / "data_utils.py")
_real = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_real)
globals().update({k: v for k, v in vars(_real).items() if not k.startswith("__")})
