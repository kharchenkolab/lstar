"""``python -m lstar`` → the conversion CLI (same entry point as the ``lstar`` console script)."""
import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
