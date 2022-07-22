#!/bin/bash

set -ex
set -o pipefail

# Update and install dependecies
echo "Update and install dependecies"
sudo apt-get install make
sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev \
    squashfs-tools libseccomp-dev wget pkg-config git cryptsetup debootstrap
sudo apt-get install bwa
sudo apt-get install 

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

WORKDIR=${PWD}
WORKFLOW_INPUT_DIR=${WORKDIR}/Files

if [ -d ${WORKFLOW_INPUT_DIR} ]; then
    echo "Directory Files exists..."
else
    mkdir -p "${WORKDIR}" "${WORKFLOW_INPUT_DIR}"
    cd ${WORKFLOW_INPUT_DIR}
    wget_download https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip vntr_data_genic_loci.zip
    unzip vntr_data_genic_loci.zip
    wget_download https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz chr1.fa.gz 
    gunzip chr1.fa.gz
    #wget_download
    # Building ref index with BWA index
    echo "Building BWA indexs..."
    bwa index -a bwtsw chr1.fa
fi

cd ${WORKFLOW_INPUT_DIR}

if [ ! -f ${WORKFLOW_INPUT_DIR}/chr1.fa ]; then
   cd ${WORKFLOW_INPUT_DIR}
   wget_download https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip vntr_data_genic_loci.zip
   unzip vntr_data_genic_loci.zip
   wget_download https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz chr1.fa.gz
   gunzip chr1.fa.gz
    #wget_download
    # Building ref index with BWA index
    echo "Building BWA indexs..."
    bwa index -a bwtsw chr1.fa
else
   echo "databases already exist..."

fi

# Installing singularity 

if ! [ -x "$(command -v singularity)" ]; then
  echo 'Singularity is not installed.'
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
else
    echo "Singularity is already instaled..."
fi

SCRIPT=${WORKDIR}/Scripts

# Installing Kestrel 
if [ ! -d "${SCRIPT}"/kestrel-1.0.1 ]; then
    cd ${SCRIPT}
    echo "Downloading Kestrel v1.0.1 ..."
    wget_download https://github.com/paudano/kestrel/releases/download/1.0.1/kestrel-1.0.1-linux.tar.gz kestrel-1.0.1-linux.tar.gz
    tar -xvzf kestrel-1.0.1-linux.tar.gz
    mv kestrel-1.0.1-linux/kestrel-1.0.1 . && rm -r kestrel-1.0.1-linux/ && rm kestrel-1.0.1-linux.tar.gz
    wget http://opengene.org/fastp/fastp
    chmod u+x ./fastp
else
    echo "Kestrel already exists..."

fi

# Building singularity image for code-adVNTR
if [ ! -f "${SCRIPT}"/code-adVNTR.sif ]; then
    cd ${SCRIPT}
    echo "Building singularity image file for code-adVNTR..."
    sudo singularity build code-adVNTR.sif adVNTR.def
else
    echo "adVNTR signularity image already exists..."

fi

# Picard tools
if [ ! -f "${SCRIPT}"/picard.jar ]; then
    cd ${SCRIPT}
    echo "downloading picard tools..."
    wget_download https://github.com/broadinstitute/picard/releases/download/2.27.3/picard.jar picard.jar 
else
    echo "Picard tools already exists..."

fi
