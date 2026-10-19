"""Logging setup.

Logs go to stderr so that stdout stays usable for machine-readable output.
Warnings raised by the pipeline stages are *collected into the manifest* as
well as logged, so a run's caveats survive in the output rather than only in a
terminal scrollback.
"""

from __future__ import annotations

import logging
import sys

LOGGER_NAME = "ibc_search"


def setup_logging(level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger(LOGGER_NAME)
    logger.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s")
    )
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.propagate = False
    return logger


def get_logger(suffix: str = "") -> logging.Logger:
    return logging.getLogger(f"{LOGGER_NAME}.{suffix}" if suffix else LOGGER_NAME)
