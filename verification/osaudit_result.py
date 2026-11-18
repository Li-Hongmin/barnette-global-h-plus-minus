"""Run a verification module and print its machine-readable result as JSON.

Usage: python verification/osaudit_result.py MODULE
MODULE is a verification/check_*.py file. Its assertions must already have
executed on import (or in `result()`) and it must expose a result dictionary
or a zero-argument result function as `result`.
The module must use no external inputs beyond its source.
"""

from __future__ import annotations

import contextlib
import importlib.util
import json
import sys
from pathlib import Path


def main() -> None:
    path = Path(sys.argv[1])
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = module
    with contextlib.redirect_stdout(sys.stderr):
        spec.loader.exec_module(module)
        value = module.result() if callable(module.result) else module.result
    print(json.dumps(value))


if __name__ == "__main__":
    main()
