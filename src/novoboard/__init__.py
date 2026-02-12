"""NovoBoard - Framework for evaluating de novo peptide sequencing."""

import logging

__version__ = "0.1.0"

# Set up library-level logging with null handler (best practice for libraries)
logging.getLogger(__name__).addHandler(logging.NullHandler())


def setup_logging(
    level: int = logging.INFO,
    format_str: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
) -> None:
    """Configure logging for NovoBoard.

    Call this from CLI or scripts to enable console logging.

    Args:
        level: Logging level (default: INFO)
        format_str: Log message format string
    """
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(format_str))

    logger = logging.getLogger(__name__)
    logger.setLevel(level)
    logger.addHandler(handler)
