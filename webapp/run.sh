#!/usr/bin/env bash
# Launch the veggenes web app on the veggenes conda env.
set -e
PORT="${1:-8050}"
PY=/usersoftware/peerd/kuov/miniforge3/envs/veggenes/bin/python
cd "$(dirname "$0")"
echo "veggenes -> http://127.0.0.1:${PORT}"
exec "$PY" app.py --port "$PORT"
