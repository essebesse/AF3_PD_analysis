#!/bin/bash
# Launch PulldownRanker
cd "$(dirname "$0")"
source .venv/bin/activate
streamlit run app.py "$@"
