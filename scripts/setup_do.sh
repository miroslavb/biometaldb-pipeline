#!/bin/bash
# Bootstrap script for DigitalOcean droplet (Ubuntu 24.04)
# Installs Python, RDKit, clones repo, prepares for batch enrichment.
set -euo pipefail

echo "=== BiometalDB Batch Enrichment - DO Setup ==="

# System deps
apt-get update -qq
apt-get install -y -qq python3 python3-pip python3-venv git sqlite3

# Create venv
python3 -m venv /opt/biometaldb-env
source /opt/biometaldb-env/bin/activate

# Python deps
pip install --quiet rdkit-pypi

# Clone repo
cd /opt
git clone https://github.com/miroslavb/biometaldb-pipeline.git
cd biometaldb-pipeline

# Create data dir and expect DB to be uploaded there
mkdir -p data

echo "=== Setup complete ==="
echo "Upload DB: scp biometaldb.sqlite root@<IP>:/opt/biometaldb-pipeline/data/"
echo "Run: source /opt/biometaldb-env/bin/activate && cd /opt/biometaldb-pipeline && python scripts/batch_enrich_do.py --workers 4"
