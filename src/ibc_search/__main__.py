"""Allow ``python -m ibc_search`` as an alternative to the console script.

Useful when the package is run straight from a source checkout via
``PYTHONPATH=src``, without an install.
"""

from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
