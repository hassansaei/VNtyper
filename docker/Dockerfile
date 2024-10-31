
# Use a specific Miniforge3 version for reproducibility
FROM condaforge/miniforge3:4.12.0-0

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/vntyper/input \
    DEFAULT_OUTPUT_DIR=/opt/vntyper/output \
    REFERENCE_DIR=/opt/vntyper/references

# Update and install necessary system packages in a single RUN command to reduce layers
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zip \
        curl \
        && rm -rf /var/lib/apt/lists/*

# Install mamba for faster conda package installation
RUN conda install -y mamba -n base -c conda-forge

# Define build arguments with default values
ARG REPO_URL=https://github.com/hassansaei/VNtyper.git
ARG REPO_DIR=/opt/vntyper

# Set default repository URL and directory
ENV REPO_URL=${REPO_URL}
ENV REPO_DIR=${REPO_DIR}

# Clone the VNtyper repository using environment variables
RUN git clone ${REPO_URL} ${REPO_DIR}

# Set the working directory
WORKDIR ${REPO_DIR}

# Copy VNtyper environment files using RUN cp
RUN cp conda/environment_vntyper.yml /tmp/environment_vntyper.yml && \
    cp conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml

# Create and install conda environments
RUN mamba env create -f /tmp/environment_vntyper.yml && \
    mamba env create -f /tmp/environment_envadvntr.yml && \
    conda clean -afy && \
    rm /tmp/environment_vntyper.yml /tmp/environment_envadvntr.yml

# Install VNtyper using setup.py within the vntyper environment
RUN conda run -n vntyper pip install .

# Make adVNTR installation script executable
RUN chmod +x vntyper/dependencies/advntr/install_advntr.sh

# Install adVNTR environment
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate envadvntr && \
    bash vntyper/dependencies/advntr/install_advntr.sh"

# Install VNtyper references during build
RUN conda run -n vntyper vntyper install-references \
    --output-dir $REFERENCE_DIR \
    --config-path /opt/vntyper/config/config.json \
    --skip-indexing

# Install FastAPI, Uvicorn, Celery, Redis, and other dependencies
RUN conda run -n vntyper pip install fastapi==0.115.3 uvicorn==0.32.0 redis==5.2.0 celery==5.4.0 python-multipart==0.0.12 fastapi-limiter==0.1.6

# Set up default input, output, and reference directories
RUN mkdir -p $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR

# Create a non-root user for running the application
RUN useradd -m appuser

# Set up entrypoint using RUN cp commands instead of COPY
RUN cp docker/entrypoint.sh /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

# Copy FastAPI app code to the container using RUN cp
RUN cp -r docker/app /opt/vntyper/app

# Set ownership of directories and app to appuser
RUN chown -R appuser:appuser $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR /opt/vntyper/app

# Expose FastAPI port
EXPOSE 8000

# Switch to a non-root user
USER appuser

# Define the entry point to either run VNtyper, Celery worker, Celery Beat, or start the FastAPI server
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

# Add HEALTHCHECK to ensure the container is running correctly
# HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
#   CMD curl -f http://localhost:8000/health/ || exit 1
=======
from python:3.9-bullseye
LABEL org.opencontainers.image.authors="hassan.saeiahan@gmail.com"

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install apt-utils
RUN apt-get -y install openssh-server
RUN apt-get -y --force-yes install libssl-dev
RUN apt-get -y install libz-dev samtools muscle make gcc g++ wget git unzip python-dev libbz2-dev liblzma-dev libz-dev bwa nano
RUN apt-get -y install htop less
RUN apt-get -y install sqlite3 libsqlite3-dev
RUN apt-get -y install screen readline-common
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install python-tk

RUN apt-get install python2.7
RUN cd /tmp; wget https://bootstrap.pypa.io/pip/2.7/get-pip.py ; python2 get-pip.py
RUN pip install -U setuptools
RUN python2 -m pip install biopython==1.76 pysam>=0.9.1.4 networkx==1.11 scipy joblib scikit-learn==0.20.4 cython==0.29.15 pysqlite

RUN apt-get install -y parallel 
RUN apt -y update 
RUN apt install -y openjdk-11-jdk-headless
RUN cd /tmp; wget https://github.com/paudano/kestrel/releases/download/1.0.1/kestrel-1.0.1-linux.tar.gz ; cd /usr/local/lib ;  tar -xvzf /tmp/kestrel-1.0.1-linux.tar.gz ; ln -s /usr/local/lib/kestrel-1.0.1/kestrel /usr/local/bin/kestrel ; rm /tmp/kestrel-1.0.1-linux.tar.gz
RUN cd /tmp ; wget http://opengene.org/fastp/fastp ; chmod ugo+x fastp; mv fastp /usr/local/bin
RUN mkdir /usr/local/lib/picard_2.27.4 ; cd /usr/local/lib/picard_2.27.4 ; wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar

RUN pip3 install --upgrade pip
RUN pip3 install grpcio
RUN pip3 install "regex<2022.1.18"
RUN pip3 install pandas
RUN pip3 install biopython
RUN pip3 install pysam

RUN apt-get install python3

RUN cd /tmp; git clone https://github.com/mehrdadbakhtiari/adVNTR.git --branch enhanced_hmm ; cd adVNTR ;  make; make install ; python2 setup.py install

COPY VNtyper_1.3.0.tar /tmp
RUN mkdir /SOFT; cd /SOFT; tar -xvf /tmp/VNtyper_1.3.0.tar; mv VNtyper_1.3.0 VNtyper
RUN cd /SOFT; wget https://github.com/biod/sambamba/releases/download/v0.6.8/sambamba-0.6.8-linux-static.gz ; gunzip sambamba-0.6.8-linux-static.gz ; mv sambamba-0.6.8-linux-static sambamba-0.6.8 ; chmod +x sambamba-0.6.8 ; ln -s /opt/sambamba-0.6.8 /usr/local/bin
RUN cd /SOFT/VNtyper/Files ; wget https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip ; unzip vntr_data_genic_loci.zip ; cd vntr_data/ ; mv hg19_genic_VNTRs.db ../.


RUN echo 'tzdata tzdata/Areas select Europe' | debconf-set-selections
RUN echo 'tzdata tzdata/Zones/Europe select Paris' | debconf-set-selections

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile
#RUN groupadd -g <groupID numeric ex 3002> <groupName string ex gruser>
#RUN useradd -m -u <userID numeric ex 6003> -g <groupID numeric ex 3002> -s /bin/bash <userName string ex jdoe>
RUN mkdir /var/run/sshd
#RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

#RUN echo '<userName string ex jdoe>:<encrypted password>' | chpasswd -e
EXPOSE 22

ENTRYPOINT ["python3","/SOFT/VNtyper/VNtyper.py"]
#CMD ["--bam", "-p", "/SOFT/VNtyper/", "-ref", "/SOFT/VNtyper/Files/chr1.fa", "-ref_VNTR", "/SOFT/VNtyper/Files/MUC1-VNTR.fa", "-m", "/SOFT/VNtyper/Files/hg19_genic_VNTRs.db"]

### docker build -t vntyper:1.3.0 .
### docker run -d -p <local address IP>:<local PORT>:<internal port ex 22> -v <PATH to local data>:/WORKSPACE --rm --name VNtyper_ssh vntyper:1.3.0

