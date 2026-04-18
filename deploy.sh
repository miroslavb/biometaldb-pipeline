#!/bin/bash
set -e

DEPLOY_DIR="/root/.hermes-agent2/biometaldb"
PYTHON="/root/.hermes-agent2/hermes-agent/venv/bin/python3"
PORT=8502

cd "$DEPLOY_DIR"

echo "=== Installing system dependencies ==="
apt-get install -y openbabel python3-openbabel avogadro 2>/dev/null || true

echo "=== Pulling latest from GitHub ==="
git pull origin master

echo "=== Stopping old server ==="
PID=$(ss -tlnp | grep ":$PORT" | grep -oP 'pid=\K[0-9]+' || true)
if [ -n "$PID" ]; then
    kill "$PID" 2>/dev/null || true
    sleep 1
    # Force kill if still running
    kill -9 "$PID" 2>/dev/null || true
    echo "Stopped PID $PID"
else
    echo "No running server found on port $PORT"
fi

echo "=== Starting server ==="
nohup "$PYTHON" mol_server.py > /tmp/mol_server.log 2>&1 &
sleep 2

# Verify
NEW_PID=$(ss -tlnp | grep ":$PORT" | grep -oP 'pid=\K[0-9]+' || true)
if [ -n "$NEW_PID" ]; then
    echo "✅ Server running on PID $NEW_PID (port $PORT)"
else
    echo "❌ Server failed to start — check /tmp/mol_server.log"
    exit 1
fi
