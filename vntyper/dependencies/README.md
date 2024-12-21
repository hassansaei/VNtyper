# Kestrel Dependency

## Overview

This directory contains the Kestrel JAR file required for genotyping within the `vntyper` environment. By including the JAR file directly in the repository, users can utilize Kestrel without needing separate installations.

## Note
We use version 1.0.1 as vntyper has been calibrated with this version and there are differences in the output of version 1.0.1 and 1.0.2.

## Included Files

- `kestrel.jar`: The executable JAR file for Kestrel v1.0.1.

## Commands to Download and Extract Kestrel

```bash
cd vntyper/dependencies
wget https://github.com/paudano/kestrel/releases/download/1.0.1/kestrel-1.0.1-linux.tar.gz
tar -xzf kestrel-1.0.1-linux.tar.gz
rm kestrel-1.0.1-linux.tar.gz
mv kestrel-1.0.1 kestrel
```
