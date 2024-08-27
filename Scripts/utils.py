import os
import logging

def setup_logging(log_file):
    logging.basicConfig(filename=log_file, level=logging.DEBUG, encoding='utf-8', format="%(asctime)s %(message)s")
    
def create_output_directories(working_dir, output):
    output_dir = os.path.join(working_dir, output)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    temp_dir = os.path.join(output_dir, "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    return output_dir, temp_dir

def search(regex: str, df, case=False):
    textlikes = df.select_dtypes(include=[object, "object"])
    return df[
        textlikes.apply(
            lambda column: column.str.contains(regex, regex=True, case=case, na=False)
        ).any(axis=1)
    ]
