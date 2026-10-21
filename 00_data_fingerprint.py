#!/usr/bin/env python3
"""
00 - Fingerprint the input data.

The inventory is edited between analysis runs. Without a record of which
version of the shapefiles produced a given number, two reports can disagree
for no reason other than that the data changed in between: during the
correction round the source/deposit overlaps in DRV_GLIM_O_I were corrected
while the checks were being run, and the same check produced five findings and
then none.

This script writes the content hash, size and feature count of every input
file, plus a single digest for the whole state. Diff the CSV against the one in
the reference outputs to see exactly which layer differs.

Output: outputs/00_data_fingerprint.txt
        outputs/00_data_fingerprint.csv
"""

import os

import glim
from glim import Report


def main():
    rep = Report("00_data_fingerprint", "00 - Input data fingerprint")
    rows = []
    for label, path in glim._input_files():
        stat = os.stat(path)
        digest = glim.sha256(path)
        features = ""
        if label.endswith(".shp"):
            directory = os.path.dirname(path)
            try:
                features = len(glim.read_layer(directory, os.path.basename(path)))
            except Exception as error:  # a layer we cannot open is itself a finding
                features = f"ERROR: {error}"
        rows.append([label, stat.st_size, digest[:16], features])
    rep.h2("Input files")
    rep.table(["file", "bytes", "sha256 (first 16)", "features"], rows)
    rep.blank()
    rep.kv("files fingerprinted", len(rows))
    rep.kv("data state digest", glim.data_state()["digest"])
    rep.blank()
    rep.line("Comparing this table with the one in the reference outputs shows which")
    rep.line("layers, if any, differ from the published data.")
    glim.write_csv("00_data_fingerprint",
                   ["file", "bytes", "sha256_16", "features"], rows)
    rep.save()


if __name__ == "__main__":
    main()
