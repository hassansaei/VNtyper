#!/bin/bash
# Update and install dependecies
echo "Update and install dependecies"
sudo apt-get update
sudo apt-get install make
sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev \
    squashfs-tools libseccomp-dev wget pkg-config git cryptsetup debootstrap
# install the Go language
echo "installing the Go language"
wget https://dl.google.com/go/go1.13.linux-amd64.tar.gz
sudo tar --directory=/usr/local -xzvf go1.13.linux-amd64.tar.gz
export PATH=/usr/local/go/bin:$PATH
# Singularity source code 
echo "Downloading Singularity source code"
wget https://github.com/singularityware/singularity/releases/download/v3.5.3/singularity-3.5.3.tar.gz
tar -xzvf singularity-3.5.3.tar.gz
# build and install!
echo "build and install"
cd singularity
./mconfig
cd builddir
make
sudo make install
# Tab completion
. etc/bash_completion.d/singularity
sudo cp etc/bash_completion.d/singularity /etc/bash_completion.d/
# Running as test
echo "Running singularity..."
singularity run library://godlovedc/funny/lolcow
