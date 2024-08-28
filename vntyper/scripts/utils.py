import os
import json
import logging

def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Set up the logging configuration for the application.

    This function configures the root logger with the specified log level 
    and optional log file. If no log file is provided, logs are sent to stdout.

    Args:
        log_level (int): The logging level (default: logging.INFO).
        log_file (str): The log file path. If None, logs are sent to stdout.

    Returns:
        None
    """
    logger = logging.getLogger()  # Get the root logger
    logger.setLevel(log_level)  # Set the overall logging level
    
    # Clear existing handlers to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # Create formatter for log messages
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    if log_file:
        # Create file handler if log_file is specified
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    # Always add a single stream handler (console output)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

def create_output_directories(base_output_dir):
    """
    Create the necessary output directories for different analyses within the pipeline.
    
    This function ensures that directories are created for Kestrel, adVNTR, FASTQ/BAM processing, and alignment processing.
    
    Args:
        base_output_dir (str): The base output directory.
    
    Returns:
        dict: A dictionary containing paths to the specific subdirectories for each analysis.
    """
    # Define subdirectory paths
    dirs = {
        "base": base_output_dir,
        "kestrel": os.path.join(base_output_dir, "kestrel"),
        "advntr": os.path.join(base_output_dir, "advntr"),
        "fastq_bam_processing": os.path.join(base_output_dir, "fastq_bam_processing"),
        "alignment_processing": os.path.join(base_output_dir, "alignment_processing"),
    }
    
    # Create each directory if it doesn't exist
    for dir_path in dirs.values():
        try:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
                logging.info(f"Created directory: {dir_path}")
            else:
                logging.info(f"Directory already exists: {dir_path}")
        except Exception as e:
            logging.error(f"Failed to create directory {dir_path}: {e}")
            raise

    return dirs

def search(regex: str, df, case=False):
    """
    Search for a regex pattern across all text-like columns in a DataFrame.

    This function filters rows in the DataFrame that match the provided regular expression 
    in any text-like column.

    Args:
        regex (str): The regular expression pattern to search for.
        df (pd.DataFrame): The DataFrame to search within.
        case (bool): Whether the search should be case-sensitive (default: False).

    Returns:
        pd.DataFrame: A filtered DataFrame containing rows that match the regex.
    """
    logging.debug("Starting regex search in DataFrame.")
    try:
        # Select only text-like columns (object dtype) for searching
        textlikes = df.select_dtypes(include=[object, "object"])
        
        # Apply the regex search across the DataFrame and filter rows that match in any column
        result_df = df[
            textlikes.apply(
                lambda column: column.str.contains(regex, regex=True, case=case, na=False)
            ).any(axis=1)
        ]
        logging.debug("Regex search completed.")
        return result_df
    except Exception as e:
        logging.error(f"Error during regex search: {e}")
        raise

def load_config(config_path=None):
    """
    Load a configuration file in JSON format.

    This function reads a JSON configuration file and returns its contents 
    as a dictionary. If no config path is provided, it returns None.

    Args:
        config_path (str): Path to the JSON configuration file.

    Returns:
        dict: The configuration parameters loaded from the file.
        None: If no config file path is provided.

    Raises:
        FileNotFoundError: If the configuration file is not found.
        json.JSONDecodeError: If there is an error decoding the JSON file.
    """
    if config_path is None:
        logging.info("No config file path provided, proceeding without configuration.")
        return None
    
    if not os.path.exists(config_path):
        logging.error(f"Config file not found: {config_path}")
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    try:
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)
            logging.info(f"Configuration loaded from {config_path}")
            return config
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from the config file: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error loading config file {config_path}: {e}")
        raise
