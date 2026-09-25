"""Shared filesystem helpers: input-file validation and work-directory."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)


def require_file(path: Path, label: str = "File") -> Path:
    """Resolve path and raise if it is missing or not a regular file.

    ``label`` names the file's role in the error message (e.g. ``"GROMACS
    coordinate file"``) so failures point at the offending input.
    """
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if not path.is_file():
        raise ValueError(f"{label} path is not a file: {path}")
    return path


def resolve_work_dir(explicit: Path | None, *, prefix: str) -> Path:
    """Return ``explicit`` or a temporary directory.

    Auto-created temporary directories are logged so callers who omitted
    ``work_dir`` can still locate their outputs; they are not cleaned up
    automatically.
    """
    if explicit is not None:
        explicit.mkdir(parents=True, exist_ok=True)
        return explicit
    work_dir = Path(tempfile.mkdtemp(prefix=prefix))
    logger.info("No work_dir given; writing outputs to temporary directory %s", work_dir)
    return work_dir
