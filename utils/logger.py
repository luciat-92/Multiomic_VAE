import logging
from pathlib import Path

def setup_logging(log_file="app.log"):
    """
    Set up logging configuration.
    """
    log_file_path = Path("logs") / log_file
    log_file_path.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=log_file_path,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

def get_logger(name="my_logger"):
    """
    Get a logger instance with the specified name.
    """
    return logging.getLogger(name)
