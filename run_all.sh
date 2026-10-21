#!/usr/bin/env bash
# Run the whole technical-validation pipeline and leave the reports in outputs/.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p outputs
for script in 00_data_fingerprint.py 01_inventory_stats.py 02_qc_database.py 03_area_statistics.py \
              04_compare_jom.py 05_compare_tosevski.py \
              06_depth_rule_records.py 07_qc_topology.py \
              09_validation_tables.py 10_specification_findings.py \
              11_landslide_areas.py 12_overlap_map.py; do
    echo
    echo "### $script"
    python3 "$script"
done
echo
echo "Reports and tables written to $(pwd)/outputs"
