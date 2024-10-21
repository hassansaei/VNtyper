#!/bin/bash
set -e

# If the first argument is 'vntyper', execute VNtyper locally
if [ "$1" = "vntyper" ]; then
  exec conda run -n vntyper "$@"
else
  # Otherwise, start the FastAPI web service
  exec uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
fi
