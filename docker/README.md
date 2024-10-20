# VNtyper Docker Container

A lightweight Docker container for **VNtyper**, enabling easy execution of the tool with customizable configurations and external data handling.

## **Building the Docker Image**

1. **Clone the VNtyper Repository:**

   ```bash
   git clone https://github.com/berntpopp/VNtyper.git
   cd VNtyper
   ```

2. **Create the entrypoint.sh Script:**

Save the entrypoint.sh content provided above into a file named entrypoint.sh in the VNtyper directory.

3. **(Optional) Create the container_config.json File:**

If you wish to use a configuration file, save the provided JSON content into container_config.json.

4. **Build the Docker Image:**

   ```bash
    docker build --build-arg REPO_URL=https://github.com/berntpopp/VNtyper.git \
                --build-arg REPO_DIR=/opt/vntyper \
                -t vntyper:2.0.0-alpha.5 .
   ```

## **Running the Docker Container**

### **Basic Help Command**
Display VNtyper help:

   ```bash
   docker run --rm vntyper:2.0.0-alpha.5 --help
   ```