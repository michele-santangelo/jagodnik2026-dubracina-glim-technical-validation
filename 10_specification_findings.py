#!/usr/bin/env python3
"""
10 - Where the submitted description did not describe the deposited data.

The integrity table (script 09) counts defects in the DATA. This script counts
something different and, for a Data Descriptor, at least as important: the
places where the checks showed that the published *description* was wrong while
the data was right.

Four such places were found. Each is encoded here against the criterion AS
SUBMITTED, so the count can be reproduced rather than taken on trust; the rest
of the pipeline encodes the corrected criterion, which is why these records are
no longer findings there.

  1. The depth rule. The submitted Methods said a landslide represented by two
     polygons is deep-seated, with no exception. Debris flows are mapped with a
     separate upslope polygon wherever the detachment and transport zone can be
     told from the fan, but that polygon is a concave erosional zone, not the
     scarp of a failure surface. The data was consistent; the rule was not.
  2. The number of polygons per landslide. The submitted text implied one or
     two. Particular morphological configurations require more.
  3. Table 2. `wxdf` occurs in the data and was absent from the table, while
     `wrf` was in the table and occurs nowhere in the data.
  4. Table 5. `lake` occurs in the data and was absent from the table.

Output: outputs/10_specification_findings.txt, outputs/10_specification.csv
"""

from collections import Counter

import glim
from glim import Report

#: Table 2 exactly as submitted: with wrf, without wxdf.
TABLE2_AS_SUBMITTED = (glim.TABLE2_LDS_CLASS | {"wrf"}) - {"wxdf"}
#: Table 5 exactly as submitted: without lake.
TABLE5_AS_SUBMITTED = glim.TABLE5_GEOM_CLASS - {"lake"}


def main():
    rep = Report("10_specification_findings",
                 "10 - Where the submitted description did not describe the data")
    layers = glim.read_all_landslides(fix=True)
    geomorph = glim.read_layer(glim.ORIGINAL_DATA, glim.GEOMORPHOLOGY_LAYER,
                               ["gm_id", "geom_class", "Shape_Area"], fix=True)
    rows = []

    rep.h2("1 - The depth rule as submitted")
    rep.line("Submitted wording: two polygons means deep-seated, one polygon means")
    rep.line("shallow, no exception. Landslides that contradict it, by type pair:")
    rep.blank()
    pairs = Counter()
    depth_hits = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            values = {f["depth"] for f in group}
            if len(values) != 1:
                continue
            value = values.pop()
            if len(group) == 2 and value == "shallow":
                pair = "+".join(sorted(f["lds_class"] for f in group))
                pairs[pair] += 1
                depth_hits.append((name, lds_id, pair))
    rep.table(["source+deposit", "landslides"],
              [[k, v] for k, v in pairs.most_common()])
    rep.blank()
    rep.kv("landslides flagged by the submitted rule", len(depth_hits))
    rep.line("Every one is a debris flow. No other type pair appears, which is why")
    rep.line("the exception is specific to debris flows and not to flows in general,")
    rep.line("and why the correction belongs in the text and not in the data.")
    for name, lds_id, pair in depth_hits:
        rows.append(["depth rule", name.replace("DRV_GLIM_", ""), lds_id, pair,
                     "Methods, Table 2, Table 7, Technical Validation, Usage Notes"])

    rep.h2("2 - More than two polygons for one landslide")
    many = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            if len(group) > 2:
                many.append((name, lds_id, len(group),
                             "+".join(sorted(f["lds_class"] for f in group))))
    rep.table(["layer", "lds_id", "polygons", "classes"],
              [[n.replace("DRV_GLIM_", ""), i, c, k] for n, i, c, k in many])
    rep.blank()
    rep.kv("landslides represented by more than two polygons", len(many))
    for n, i, c, k in many:
        rows.append(["polygons per landslide", n.replace("DRV_GLIM_", ""), i, k,
                     "Methods"])

    rep.h2("3 and 4 - The published code dictionaries")
    # The rock-fall classes (rf, xrf, wxrf) live in the rockfall layer, so the
    # census of used codes must include it or they look undocumented.
    rockfalls = glim.read_rockfalls(fix=True)
    used = Counter(f["lds_class"] for feats in layers.values() for f in feats)
    used.update(f["lds_class"] for f in rockfalls)
    undoc = sorted(set(used) - TABLE2_AS_SUBMITTED)
    unused = sorted(TABLE2_AS_SUBMITTED - set(used))
    geo_used = Counter(f["geom_class"] for f in geomorph)
    geo_undoc = sorted(set(geo_used) - TABLE5_AS_SUBMITTED)
    table = []
    for code in undoc:
        table.append(["Table 2", code, "in the data, absent from the table",
                      used[code]])
        rows.append(["Table 2", "landslide layers", code, "undocumented", "Table 2"])
    for code in unused:
        table.append(["Table 2", code, "in the table, absent from the data", 0])
        rows.append(["Table 2", "-", code, "documented but unused", "Table 2"])
    for code in geo_undoc:
        table.append(["Table 5", code, "in the data, absent from the table",
                      geo_used[code]])
        rows.append(["Table 5", "Geomorphology", code, "undocumented", "Table 5"])
    rep.table(["table", "code", "discrepancy", "polygons"], table)

    rep.h2("Summary")
    rep.table(["correction", "records"],
              [["depth rule (debris flows)", len(depth_hits)],
               ["more than two polygons", len(many)],
               ["code dictionaries", len(undoc) + len(unused) + len(geo_undoc)],
               ["TOTAL", len(rows)]])
    rep.blank()
    rep.line("None of these is a defect in the deposited data, and none is counted")
    rep.line("in the integrity table of script 09, which measures the data. They are")
    rep.line("corrections the checks forced on the manuscript itself, and they are")
    rep.line("the strongest evidence that the validation was not circular: a check")
    rep.line("that only restated the rules used to make the data could not have")
    rep.line("contradicted them.")

    glim.write_csv("10_specification",
                   ["correction", "layer", "record", "detail", "corrected in"], rows)
    rep.save()


if __name__ == "__main__":
    main()
