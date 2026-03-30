#!/bin/bash
# Launch AF3 Pulldown Analysis App
cd "$(dirname "$0")"
source .venv/bin/activate
streamlit run app.py "$@"
