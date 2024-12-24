# Installing adVNTR for VNtyper

## Overview

This document provides detailed instructions for installing adVNTR, which is used by VNtyper for genotyping Variable Number Tandem Repeats (VNTRs), specifically the MUC1 VNTR. You can choose between an automated installation using the provided helper script or a manual installation. Additionally, you have the option to specify a conda environment during installation.

---

## Installation Methods

- [Automated Installation with Helper Script](#automated-installation-with-helper-script)
- [Manual Installation](#manual-installation)
- [References](#references)

---

## Automated Installation with Helper Script

We have provided a helper script `install_advntr.sh` to automate the installation process. The script allows you to optionally specify a conda environment to activate before installation.

### Prerequisites:

- **Conda Environment**: Ensure that the conda environment exists if you plan to specify one. The script does **not** create conda environments.
- **Python Version**: adVNTR requires Python 2.7 or Python 3.7.

### Steps:

1. **Create the Conda Environment (Optional)**

   If you haven't already created a conda environment for adVNTR, you can do so using your existing `environment_envadvntr.yml` file:

   ```bash
   conda env create -f environment_envadvntr.yml
   ```
   Activate the environment:
   ```bash
   conda activate envadvntr
   ```

2. **Navigate to the Dependencies Directory**

   ```bash
   cd vntyper/dependencies/advntr
   ```

3. **Run the Helper Script**

   ```bash
   bash install_advntr.sh
   ```

   **Options:**

   - `-e`, `--env`: Name of the conda environment to activate (optional).
   - `-d`, `--install-dir`: Directory where adVNTR will be installed (default: `$PWD/adVNTR`).
   - `-o`, `--overwrite`: Overwrite the installation directory if it exists.
   - `-h`, `--help`: Display help message.

   **Examples:**

   - **Install adVNTR and activate a conda environment:**
     ```bash
     bash install_advntr.sh -e envadvntr
     ```

   - **Install adVNTR to a custom directory and activate a conda environment:**
     ```bash
     bash install_advntr.sh -e envadvntr -d /path/to/install/adVNTR
     ```

   - **Overwrite existing installation:**
     ```bash
     bash install_advntr.sh -o
     ```

4. **Verify Installation**

   After the script completes, adVNTR should be installed in the specified directory. You can verify the installation:
   ```bash
   advntr --help
   ```
   If `advntr` is not found, ensure that your conda environment is activated and that the `bin` directory of the environment is in your `PATH`.

---

## Manual Installation

If you prefer to install adVNTR manually, follow these steps.

### Prerequisites:

- **Conda Environment**: Ensure that your conda environment is created and activated.
  - Use your existing `environment_envadvntr.yml`:
    ```bash
    conda env create -f environment_envadvntr.yml
    conda activate envadvntr
    ```

### Steps:

1. **Clone the adVNTR Repository**

   ```bash
   git clone https://github.com/mehrdadbakhtiari/adVNTR.git --branch enhanced_hmm
   cd adVNTR
   ```

2. **Install adVNTR**

   ```bash
   python setup.py install
   ```

3. **Verify Installation**

   ```bash
   advntr --help
   ```

---

## References

We have **removed** reference installation from the `install_advntr.sh` script to avoid duplicating or mismatching reference data. 

To install **references** for VNtyper (including those used by adVNTR), use the dedicated command:
```bash
vntyper install-references -d /path/to/install/references --config-path /path/to/config.json
```
This leverages the `install_references.py` script and its associated `install_references_config.json`.

For more details on reference installation:
- See the script at **`vntyper/scripts/install_references.py`**  
- Or check **`install_references_config.json`** for configuration and URLs.

---

**Note:**  
For integration with VNtyper, ensure that the `advntr` command is accessible in your environment. You may need to update VNtyper’s configuration files to point to the correct adVNTR installation path if necessary.
