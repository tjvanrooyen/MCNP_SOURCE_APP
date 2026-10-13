"""Run the MCNP_SOURCE_APP regression suite.

    python -m tests.run_tests            (from the project root)

Exit status is 0 when every test passes and 1 otherwise, so the suite can gate
a release of the code or of a refreshed ICRP-107 data set.
"""

from __future__ import annotations

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def main() -> int:
    suite = unittest.defaultTestLoader.discover(
        start_dir=os.path.join(ROOT, "tests"), top_level_dir=ROOT
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
