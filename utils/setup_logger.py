import logging

def setup_logger(name: str, log_file: str = None, level: int = logging.INFO) -> logging.Logger:
    """
    Sets up and returns a logger with the specified name, log file, and level.
    Logs to console by default, and optionally to a file.

    Args:
        name (str): Name of the logger.
        log_file (str, optional): Path to the log file. If None, file logging is disabled.
        level (int, optional): Logging level. Default is logging.INFO.

    Returns:
        logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler (optional)
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    # Avoid duplicate logs
    logger.propagate = False

    return logger