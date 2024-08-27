import os
import json
import logging

def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Set up logging configuration.

    Args:
        log_level: The logging level (default: logging.INFO).
        log_file: The log file path. If None, logs to stdout.
    """
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        filename=log_file,
        filemode='a',  # Append mode for logging to a file
    )
    
    if log_file is None:  # If no log file is specified, also log to stdout
        logging.getLogger().addHandler(logging.StreamHandler())

def create_output_directories(working_dir, output):
    output_dir = os.path.join(working_dir, output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Changed to makedirs to create intermediate directories if they don't exist
    temp_dir = os.path.join(output_dir, "temp")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)  # Changed to makedirs to ensure full directory creation
    return output_dir, temp_dir

def search(regex: str, df, case=False):
    textlikes = df.select_dtypes(include=[object, "object"])
    return df[
        textlikes.apply(
            lambda column: column.str.contains(regex, regex=True, case=case, na=False)
        ).any(axis=1)
    ]

def load_config(config_path=None):
    """
    Load configuration from a JSON file.
    
    Args:
        config_path: Path to the configuration JSON file.
        
    Returns:
        A dictionary containing configuration parameters, or None if no config path is provided.
    """
    if config_path is None:
        logging.info("No config file path provided, proceeding without configuration.")
        return None
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    try:
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)
        return config
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from the config file: {e}")
        raise
