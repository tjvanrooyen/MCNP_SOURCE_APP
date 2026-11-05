"""Run the MCNP_SOURCE_APP regression suite.

From the project folder::

    python _Archive_Important/tests/run_tests.py

The suite finds the project folder by itself -- the nearest folder above this one
that contains ``msa_modules`` -- so it runs wherever the ``tests`` folder is
kept (archived as above, or copied back to the project root, where
``python -m tests.run_tests`` also works).

Exit status is 0 when every test passes and 1 otherwise, so the suite can gate
a release of the code or of a refreshed ICRP-107 data set.
"""

from __future__ import annotations

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def find_project_root(start: str) -> str:
    """The nearest folder at or above ``start`` that contains ``msa_modules``."""
    folder = start
    while not os.path.isdir(os.path.join(folder, "msa_modules")):
        parent = os.path.dirname(folder)
        if parent == folder:
            raise SystemExit(f"Cannot find the MCNP_SOURCE_APP project folder above {start}")
        folder = parent
    return folder


ROOT = find_project_root(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def main() -> int:
    # The folder above ``tests`` is the top level, so the package imports as
    # ``tests.<module>`` wherever it lives.
    suite = unittest.defaultTestLoader.discover(
        start_dir=HERE, top_level_dir=os.path.dirname(HERE)
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
