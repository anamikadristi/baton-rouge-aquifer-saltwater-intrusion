#!/usr/bin/env bash
# Run the full pipeline from the repository root:  bash scripts/run_all.sh
# Optional: export BRAQ_RAW_DIR=/some/folder to keep raw downloads outside the repo.
set -euo pipefail
cd "$(dirname "$0")"
python 01_download_data.py
python 02_prepare_wells_levels.py
python 04_fault_and_withdrawals.py   # fault line is needed by step 3
python 03_prepare_chloride.py
python 05_water_level_trends.py
python 06_chloride_trends_maps.py
python 07_intrusion_illustration.py
python 08_overview_figures_gis.py
