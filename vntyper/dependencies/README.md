# Kestrel Dependency

## Overview

This directory contains the Kestrel JAR file required for [specific functionalities] within the `vntyper` environment. By including the JAR file directly in the repository, users can utilize Kestrel without needing separate installations.

## Included Files

- `kestrel.jar`: The executable JAR file for Kestrel.

# Comands to download and extract Kestrel
```bash
cd vntyper/dependencies
wget https://github.com/paudano/kestrel/releases/download/1.0.1/kestrel-1.0.1-linux.tar.gz
tar -xzf kestrel-1.0.1-linux.tar.gz
rm kestrel-1.0.1-linux.tar.gz
mv kestrel-1.0.1 kestrel
``` 