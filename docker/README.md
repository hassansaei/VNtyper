# VNtyper Docker Container

A Docker container for **VNtyper**, enabling easy execution of the tool with customizable configurations and external data handling.

## **Building the Docker Image**

1. **Clone the VNtyper Repository:**

   ```bash
   git clone https://github.com/hassansaei/VNtyper.git
   cd VNtyper
   ```

2. **Create the `entrypoint.sh` Script:**

   Save the `entrypoint.sh` content provided above into a file named `entrypoint.sh` in the VNtyper directory.

3. **(Optional) Create the `config.json` File:**

   If you wish to use a configuration file, save the provided JSON content into `docker/config.json`.

4. **Build the Docker Image:**

   ```bash
   docker build --no-cache --build-arg REPO_URL=https://github.com/hassansaei/VNtyper.git \
               --build-arg REPO_DIR=/opt/vntyper \
               -t vntyper:latest .
   ```
5. **Pull the Docker Image from Docker Hub:**

    ```bash
    docker pull saei/vntyper:latest
    ```
6. **Generate apptainer Image from Docker Image:**

    ```bash
    apptainer pull docker://saei/vntyper:latest
    ```

## **Running the Docker Container**

### **CLI Usage**

Run docker interactively:

```bash
   docker run -w /opt/vntyper --rm \
    -v /local/input/folder/:/opt/vntyper/input \
    -v /local/output/folder/:/opt/vntyper/output \
    vntyper:latest \
    vntyper pipeline --bam /local/input/folder/filename.bam \
    -o /local/output/folder/filename/
```
Run apptainer interactively:

```bash
    apptainer run --pwd /opt/vntyper \
    -B /local/input/folder/:/opt/vntyper/input \
    -B /local/output/folder/:/opt/vntyper/output \
    vntyper_2.0.0.sif vntyper pipeline \
    --bam /opt/vntyper/input/filename.bam \
    -o /opt/vntyper/output/filename/ 
```

### **API Usage**

#### **1. Run the API Server**

Start the FastAPI server by running the container:

```bash
docker run -d -p 8000:8000 \
    -v /local/input/folder/:/opt/vntyper/input \
    -v /local/output/folder/:/opt/vntyper/output \
    vntyper:latest
```

#### **2. Submit a Job via API**

Use `curl` to upload the BAM file and pass the required parameters:

```bash
curl -X POST "http://localhost:8000/run-job/" \
    -F "file=@/local/input/folder/filename.bam" \
    -F "thread=4" \
    -F "reference_assembly=hg38" \
    -F "fast_mode=true" \
    -F "keep_intermediates=true" \
    -F "archive_results=true"
```

**Response:**

```json
{
  "message": "Job started",
  "output_dir": "/download/filename"
}
```

#### **3. Download Results**

After the job completes, download the zipped results:

```bash
curl -O "http://localhost:8000/download/filename.zip"
```

## **Notes**

- **Volume Mounts:**
  - Ensure that the local directories `/local/input/folder/` and `/local/output/folder//` exist and have appropriate read/write permissions.
  
- **API Parameters:**
  - `thread`: Number of threads to use (e.g., `4`).
  - `reference_assembly`: Reference genome assembly (e.g., `hg38`).
  - `fast_mode`: Enable fast mode (`true` or `false`).
  - `keep_intermediates`: Retain intermediate files (`true` or `false`).
  - `archive_results`: Archive results as a ZIP file (`true` or `false`).

- **Accessing the API:**
  - **Submit a Job:** Use the `/run-job/` endpoint to submit a BAM file for processing.
  - **Download Results:** Use the `/download/{output_dir}` endpoint to retrieve the processed results.

- **Health Check:**
  - The Docker container includes a health check to ensure it's running correctly. You can monitor the container's health status using:
    ```bash
    docker ps
    ```
  
- **Logs:**
  - To view the container logs for debugging:
    ```bash
    docker logs <container_id>
    ```
