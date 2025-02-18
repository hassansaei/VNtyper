# vntyper/modules/__init__.py
# This file makes 'modules' a subpackage of 'vntyper'

import importlib

# Dictionary mapping module names to their package paths
AVAILABLE_MODULES = {
    "advntr": "vntyper.modules.advntr",
    # Add future modules here
    # 'moduleX': 'vntyper.modules.moduleX',
}


def is_module_available(module_name):
    """
    Check if a module is available in the AVAILABLE_MODULES.

    Args:
        module_name (str): Name of the module to check.

    Returns:
        bool: True if the module is available, False otherwise.
    """
    return module_name in AVAILABLE_MODULES


def load_module(module_name):
    """
    Dynamically load a module by its name.

    Args:
        module_name (str): Name of the module to load.

    Returns:
        module: The loaded module object.

    Raises:
        ImportError: If the module cannot be imported.
    """
    if not is_module_available(module_name):
        raise ImportError(f"Module '{module_name}' is not available.")

    module_path = AVAILABLE_MODULES[module_name]
    return importlib.import_module(module_path)
