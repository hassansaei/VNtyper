#!/bin/bash
set -e

# Source conda to ensure the environment is activated
source /opt/conda/etc/profile.d/conda.sh

# If the first argument is 'vntyper', execute VNtyper locally
if [ "$1" = "vntyper" ]; then
  exec conda run -n vntyper "$@"
# If the first argument is 'celery', execute Celery worker
elif [ "$1" = "celery" ]; then
  exec conda run -n vntyper celery -A app.celery_app worker --loglevel=info
# If the first argument is 'beat', execute Celery Beat
elif [ "$1" = "beat" ]; then
  exec conda run -n vntyper celery -A app.celery_app beat --loglevel=info
else
  # Otherwise, start the FastAPI web service
  exec conda run -n vntyper uvicorn app.main:app --host 0.0.0.0 --port 8000
fi
