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
    if not logger.handlers:  # Ensure no duplicate handlers are added
        # Configure the logger
        logging.basicConfig(
            level=log_level,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            filename=log_file,
            filemode='a',  # Append mode for logging to a file
        )
        
        if log_file is None:  # If no log file is specified, add a StreamHandler for stdout
            logger.addHandler(logging.StreamHandler())

def create_output_directories(working_dir, output):
    """
    Create the necessary output and temporary directories for the pipeline.

    This function ensures that both the main output directory and the 
    corresponding 'temp' subdirectory are created.

    Args:
        working_dir (str): The base working directory.
        output (str): The name of the output directory to create.

    Returns:
        tuple: A tuple containing the paths to the output directory and temp directory.
    """
    # Define the output directory path
    output_dir = os.path.join(working_dir, output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create the output directory if it doesn't exist
    
    # Define the temp subdirectory path within the output directory
    temp_dir = os.path.join(output_dir, "temp")
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)  # Create the temp directory if it doesn't exist
    
    return output_dir, temp_dir

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
    # Select only text-like columns (object dtype) for searching
    textlikes = df.select_dtypes(include=[object, "object"])
    
    # Apply the regex search across the DataFrame and filter rows that match in any column
    return df[
        textlikes.apply(
            lambda column: column.str.contains(regex, regex=True, case=case, na=False)
        ).any(axis=1)
    ]

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
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    try:
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)
        return config
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from the config file: {e}")
        raise
