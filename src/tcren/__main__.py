"""Enable ``python -m tcren`` as an alias for the ``tcren`` console script."""

from .cli import app

if __name__ == "__main__":  # pragma: no cover
    app()
