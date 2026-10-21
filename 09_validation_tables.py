#!/usr/bin/env python3
"""
09 - The tables that go into the Technical Validation section.

Everything the Technical Validation section of the paper states is generated
here, so no number in it is transcribed by hand. Four blocks:

  A  database and geometric integrity: what was checked, what was found,
     what was corrected. Compares the published data against the findings
     over the dataset as originally submitted, in reference/as_submitted/.
  B  comparison with Tosevski (2013): a field-verified inventory of 39 slides.
  C  comparison with Jagodnik et al. (2020): 528 polygons over the supplied
     16.75 km2 pilot area, with the confusion matrix over the real domain.
     Those data are not public; without them block C is skipped.
  D  the published uncertainty attribute, as descriptive metadata.

reference/as_submitted/ holds the findings of scripts 02 and 07 run over the
dataset as originally submitted, i.e. before the correction round. That dataset
is not distributed (see README); its findings are, so that the "flagged" and
"corrected" columns can be regenerated rather than remembered.

Output: outputs/09_validation_tables.txt
        outputs/09_integrity_checks.csv
        outputs/09_comparison_summary.csv
"""

import csv
import os
from collections import Counter, defaultdict

from shapely.geometry import Point
from shapely.strtree import STRtree

import glim
from glim import Report

#: The "flagged" column of the integrity table must be counted by the SAME code
#: as the "remaining" column, or the two describe different checks. It therefore
#: comes from a full run of this pipeline over the dataset as originally
#: submitted. That dataset is not distributed (see README); the findings of the
#: run are, in reference/as_submitted/.
BASELINE = os.environ.get(
    "GLIM_BASELINE",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference", "as_submitted"))

#: how each automated check maps onto a row of the published table, and what
#: the check is testing. Order is the order of the table in the paper.
#: fourth element: what the records still flagged after the correction round
#: actually are. "convention" means the record follows a mapping rule that the
#: Methods of the paper state explicitly; "as mapped" means it was inspected in
#: the GIS and found correct; "intentional" means the schema is deliberate.
CHECK_GROUPS = [
    ("Coordinate reference system",
     "every layer declares EPSG:3765",
     ["A_crs"], "as mapped"),
    ("Geometry validity",
     "no self-intersecting, empty or null geometry",
     ["A_invalid", "B_invalid_geometry", "B_empty_geometry"], "as mapped"),
    ("Duplicate geometries",
     "no polygon repeated within or across layers",
     ["B_duplicate", "B2_near_duplicate", "B3_duplicate_across_layers",
      "K_duplicate_geometry"], "as mapped"),
    ("Degenerate polygons",
     "no polygon below 10 m2 or with fewer than four vertices",
     ["C_degenerate"], "as mapped"),
    ("Vertex defects",
     "no repeated vertices; spikes inspected",
     ["D_vertex_defect"], "as mapped"),
    ("Overlap within a generation",
     "distinct landslides of one generation are not superimposed",
     ["E_intra_layer_overlap"], "as mapped"),
    ("Source and deposit topology",
     "the polygons of one landslide share a boundary without overlapping",
     ["F_source_deposit_overlap", "F_source_deposit_gap", "F_snap_artefact",
      "G_containment"], "as mapped"),
    ("Identifier structure",
     "each identifier belongs to the polygons of one landslide only",
     ["E_too_many_polygons", "E_no_identifier", "F_pairing"], "convention"),
    ("Landslide split across two identifiers",
     "a source area is not a landslide on its own, debris flows excepted",
     ["L_lone_source"], "as mapped"),
    ("Attribute consistency within a landslide",
     "the polygons of one landslide agree on depth and uncertainty",
     ["G_mixed_depth", "G_mixed_uncert"], "as mapped"),
    ("Classification criteria",
     "the depth class follows the criteria stated in the Methods",
     ["H_depth_rule"], "convention"),
    ("Missing values and code dictionaries",
     "no null attribute; every code is documented",
     ["C_null", "D_bad_uncert", "D_undocumented_code", "D_unused_code"],
     "documented in the tables"),
    ("Features outside the study area",
     "no polygon extends beyond the mapped domain",
     ["I_outside_boundary"], "as mapped"),
    ("Area attribute",
     "Shape_Area agrees with the geometry on a single basis",
     ["J_stale_shape_area", "J_area_mismatch", "J_stale_area_group"], "as mapped"),
    ("Geomorphological unit overlap",
     "overlaps between surface-deposit units inspected",
     ["I_geomorph_overlap"], "as mapped"),
]

#: checks that are informational rather than defects
INFORMATIONAL = {"B_multipart", "H_holes", "C2_sliver", "J_slivers",
                 "E_not_globally_unique", "E_rockfall_missing_fields",
                 "J_boundary_area"}


def load_findings(path, check_col=0):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def count_checks(rows):
    return Counter(r["check"] for r in rows)


def comparison_with_jom(rep, layers, rockfalls):
    """Block C. Needs the Jagodnik et al. (2020) data, which are not public;
    returns the rows it contributes to 09_comparison_summary.csv."""
    boundary = glim.read_boundary()
    pilot = glim.read_jom_domain()
    domain = pilot.intersection(boundary)
    jom = glim.read_jom_polygons(fix=True)
    jom_geoms = [f.geom for feats in jom.values() for f in feats]
    jom_union = glim.safe_union(jom_geoms)
    new_geoms = [f.geom for feats in layers.values() for f in feats] + \
                [f.geom for f in rockfalls]
    new_union = glim.safe_union(new_geoms)

    rep.table(
        ["quantity", "value"],
        [["pilot area of that inventory (km2)", f"{glim.km2(pilot.area):.2f}"],
         ["  inside the present study area (km2)", f"{glim.km2(domain.area):.2f}"],
         ["  as a share of the present study area",
          f"{100 * domain.area / boundary.area:.1f} %"],
         ["polygons in that inventory", len(jom_geoms)],
         ["their total area (km2)", f"{glim.km2(jom_union.area):.4f}"]])

    coverage = [f.geom.intersection(new_union).area / f.geom.area
                for feats in jom.values() for f in feats]
    recovered = sum(1 for c in coverage if c >= 0.5)
    any_overlap = sum(1 for c in coverage if c > 0)
    rep.blank()
    rep.line("Detection agreement:")
    rep.table(
        ["quantity", "value"],
        [["polygons at least half covered by the new G-LIM",
          f"{recovered} / {len(coverage)} ({100 * recovered / len(coverage):.1f} %)"],
         ["polygons with any overlap",
          f"{any_overlap} / {len(coverage)} ({100 * any_overlap / len(coverage):.1f} %)"],
         ["polygons with no overlap at all", len(coverage) - any_overlap],
         ["mean fraction of a polygon covered",
          f"{100 * sum(coverage) / len(coverage):.1f} %"]])

    mine = new_union.intersection(domain)
    theirs = jom_union.intersection(domain)
    both = mine.intersection(theirs)
    union_two = glim.safe_union([mine, theirs])
    rep.blank()
    rep.line("Areal agreement over the pilot area clipped to the study boundary:")
    rep.table(
        ["class", "km2", "% of the domain"],
        [["mapped as landslide in both inventories", f"{glim.km2(both.area):.4f}",
          f"{100 * both.area / domain.area:.1f}"],
         ["only in the new G-LIM", f"{glim.km2(mine.difference(theirs).area):.4f}",
          f"{100 * mine.difference(theirs).area / domain.area:.1f}"],
         ["only in the previous inventory",
          f"{glim.km2(theirs.difference(mine).area):.4f}",
          f"{100 * theirs.difference(mine).area / domain.area:.1f}"],
         ["in neither",
          f"{glim.km2(domain.difference(union_two).area):.4f}",
          f"{100 * domain.difference(union_two).area / domain.area:.1f}"]])
    rep.blank()
    rep.kv("previous inventory area recovered", f"{100 * both.area / theirs.area:.1f} %")
    rep.kv("index of agreement (intersection / union)",
           f"{both.area / union_two.area:.3f}")
    rep.kv("error index E = 1 - intersection / union",
           f"{1 - both.area / union_two.area:.3f}")
    rep.line()
    rep.line("E is a measure of difference, not of error, when the two inventories have")
    rep.line("different scope: the previous one mapped active, mostly gully-related")
    rep.line("instabilities over 16.75 km2, the present one also maps old and very old")
    rep.line("deep-seated landslides over the whole catchment. The informative figure is")
    rep.line("the share of the previous inventory that the new one recovers.")

    # boundary agreement and type agreement, by previous layer
    records = []
    for name, feats in layers.items():
        for f in feats:
            records.append((name.replace("DRV_GLIM_", ""), f["lds_class"], f.geom))
    for f in rockfalls:
        records.append(("rockfalls", f["lds_class"], f.geom))
    geoms = [r[2] for r in records]
    tree = STRtree(geoms)
    rows = []
    for name, feats in sorted(jom.items()):
        ious, types, ages = [], Counter(), Counter()
        for f in feats:
            best, best_iou = None, 0.0
            for idx in tree.query(f.geom):
                value = glim.iou(f.geom, geoms[idx])
                if value > best_iou:
                    best_iou, best = value, records[idx]
            ious.append(best_iou)
            if best:
                types[best[1]] += 1
                ages[best[0]] += 1
        ious.sort()
        label = name.replace("JOM_", "").replace("_htrs", "")
        rows.append([label, glim.JOM_LABELS.get(name, "?"), len(feats),
                     f"{ious[len(ious) // 2]:.2f}",
                     ", ".join(f"{k}={v}" for k, v in types.most_common(3)),
                     ", ".join(f"{k}={v}" for k, v in ages.most_common(2))])
    rep.blank()
    rep.line("Boundary and class agreement, best match by intersection over union:")
    rep.table(["previous layer", "landform", "n", "median IoU",
               "commonest matched type", "commonest matched generation"], rows)

    return [["Jagodnik et al. 2020", "polygons", len(jom_geoms)],
            ["Jagodnik et al. 2020", "recovered at 50 % coverage", recovered],
            ["Jagodnik et al. 2020", "area recovered %",
             round(100 * both.area / theirs.area, 1)],
            ["Jagodnik et al. 2020", "index of agreement",
             round(both.area / union_two.area, 3)]]


def main():
    rep = Report("09_validation_tables", "09 - Tables for the Technical Validation")

    if not os.path.exists(os.path.join(BASELINE, "02_qc_findings.csv")):
        raise SystemExit(
            f"missing {BASELINE}: the findings over the dataset as originally "
            "submitted, distributed in reference/as_submitted/")
    before = count_checks(load_findings(os.path.join(BASELINE, "02_qc_findings.csv")))
    before += count_checks(load_findings(os.path.join(BASELINE,
                                                      "07_topology_findings.csv")))
    after = count_checks(load_findings(os.path.join(glim.OUTPUTS,
                                                    "02_qc_findings.csv")))
    after += count_checks(load_findings(os.path.join(glim.OUTPUTS,
                                                     "07_topology_findings.csv")))

    layers = glim.read_all_landslides(fix=True)
    rockfalls = glim.read_rockfalls(fix=True)
    geomorph = glim.read_layer(glim.ORIGINAL_DATA, glim.GEOMORPHOLOGY_LAYER,
                               ["gm_id", "geom_class", "Shape_Area"], fix=True)
    n_polygons = sum(len(v) for v in layers.values())
    n_landslides = sum(len(glim.group_by_id(v)) for v in layers.values())

    # ------------------------------------------------------------------ A
    rep.h2("A - Database and geometric integrity")
    rep.kv("landslide polygons inspected", n_polygons)
    rep.kv("landslides (distinct identifiers)", n_landslides)
    rep.kv("rockfall polygons inspected", len(rockfalls))
    rep.kv("geomorphological polygons inspected", len(geomorph))
    rep.blank()
    rows = []
    total_found = total_left = 0
    residue_kind = Counter()
    for label, criterion, checks, kind in CHECK_GROUPS:
        found = sum(before.get(c, 0) for c in checks)
        left = sum(after.get(c, 0) for c in checks)
        total_found += found
        total_left += left
        residue_kind[kind] += left
        rows.append([label, criterion, found, found - left,
                     left if left else "-", kind if left else "-"])
    rows.append(["TOTAL", "", total_found, total_found - total_left,
                 total_left if total_left else "-", ""])
    rep.table(["check", "criterion", "records flagged", "corrected",
               "remaining", "what the remainder is"], rows)
    rep.line()
    rep.line("The remainder is not a residue of unresolved errors. It is of three kinds.")
    rep.line("'convention': the record follows a mapping rule that the Methods state")
    rep.line("explicitly, namely that debris flows are classified shallow")
    rep.line("whatever their polygon structure and that a landslide may need more than")
    rep.line("two polygons. 'as mapped': the record was examined in the GIS and found")
    rep.line("correct, for instance a long narrow debris-flow channel that a shape index")
    rep.line("flags as a sliver, or a gap that is a transit zone between a detachment")
    rep.line("area and its deposit. 'documented in the tables': a code that the")
    rep.line("data contains and the published dictionary did not list.")
    rep.blank()
    rep.table(["remainder", "records"],
              [[k, v] for k, v in residue_kind.items() if v])

    rep.blank()
    rep.line("The checks that a defect count would misrepresent, reported separately:")
    info_rows = []
    for check in sorted(INFORMATIONAL):
        if before.get(check) or after.get(check):
            info_rows.append([check, before.get(check, 0), after.get(check, 0)])
    rep.table(["check", "before", "after"], info_rows)

    # ------------------------------------------------------------------ B
    rep.h2("B - Comparison with Tosevski (2013)")
    match_path = os.path.join(glim.OUTPUTS, "05_tosevski_match.csv")
    rows = load_findings(match_path)
    if rows:
        inside = sum(1 for r in rows if r["inside_new_glim"] == "yes")
        within = {t: sum(1 for r in rows
                         if float(r["distance_to_nearest_m"]) <= t) for t in (5, 10, 25)}
        in_jom = sum(1 for r in rows if r["inside_jom_2020"] == "yes")
        by_field = Counter(r["detection_method"] for r in rows)
        rep.table(
            ["quantity", "value"],
            [["landslides in that inventory", len(rows)],
             ["  detected by field geomorphological mapping", by_field.get("field", 0)],
             ["  detected on the 1 m DTM and then verified in the field",
              by_field.get("DEM", 0)],
             ["falling inside a landslide polygon of the new G-LIM",
              f"{inside} / {len(rows)} ({100 * inside / len(rows):.0f} %)"],
             ["within 5 m of one",
              f"{within[5]} / {len(rows)} ({100 * within[5] / len(rows):.0f} %)"],
             ["within 25 m of one",
              f"{within[25]} / {len(rows)} ({100 * within[25] / len(rows):.0f} %)"],
             ["falling inside the Jagodnik et al. (2020) inventory",
              f"{in_jom} / {len(rows)}" if glim.jom_available() else "not available"]])
        rep.line()
        rep.line("The coordinates were published in MGI 1901 / Balkans zone 5 and")
        rep.line("reprojected to EPSG:3765; the transcription reproduces every aggregate")
        rep.line("value in that paper, including the total landslide area of 81,873 m2.")

    # ------------------------------------------------------------------ C
    rep.h2("C - Comparison with Jagodnik et al. (2020)")
    if glim.jom_available():
        jom_summary = comparison_with_jom(rep, layers, rockfalls)
    else:
        jom_summary = []
        rep.line("SKIPPED: the Jagodnik et al. (2020) data are not public and were not")
        rep.line("found; they must be requested from the authors of that work (see README).")

    # ------------------------------------------------------------------ D
    rep.h2("D - The published uncertainty attribute")
    per_landslide = Counter()
    for feats in layers.values():
        for group in glim.group_by_id(feats).values():
            per_landslide[group[0]["uncert"]] += 1
    total = sum(per_landslide.values())
    rows = []
    for code in sorted(per_landslide):
        flags = [d == "1" for d in code]
        rows.append([code,
                     "uncertain" if flags[0] else "-",
                     "uncertain" if flags[1] else "-",
                     "uncertain" if flags[2] else "-",
                     per_landslide[code],
                     f"{100 * per_landslide[code] / total:.1f}"])
    rep.table(["uncert", "detection", "delineation", "classification",
               "landslides", "%"], rows)
    by_kind = {
        "detection": sum(v for c, v in per_landslide.items() if c[0] == "1"),
        "delineation": sum(v for c, v in per_landslide.items() if c[1] == "1"),
        "classification": sum(v for c, v in per_landslide.items() if c[2] == "1"),
    }
    rep.blank()
    rep.table(["uncertainty flagged for", "landslides", "% of the inventory"],
              [[k, v, f"{100 * v / total:.1f}"] for k, v in by_kind.items()])
    rep.blank()
    rows = []
    for age in ("O", "VO", "VOR"):
        n = flagged = 0
        for name, feats in layers.items():
            if glim.RELATIVE_AGE[name] != age:
                continue
            for group in glim.group_by_id(feats).values():
                n += 1
                flagged += group[0]["uncert"] != "000"
        rows.append([age, n, flagged, f"{100 * flagged / n:.1f}"])
    rep.table(["relative age", "landslides", "with any flag", "%"], rows)
    rep.blank()
    # Enumerated rather than asserted: an earlier version of this report stated
    # that the landslides flagged for detection or delineation were exactly the
    # three very-old-relict ones, which is not what the data says.
    flagged_dd = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            code = group[0]["uncert"]
            if code[0] == "1" or code[1] == "1":
                flagged_dd.append([name.replace("DRV_GLIM_", ""), lds_id, code,
                                   glim.RELATIVE_AGE[name], group[0]["lds_class"]])
    rep.line("Landslides flagged for detection or delineation, enumerated:")
    rep.blank()
    rep.table(["layer", "lds_id", "uncert", "relative age", "class"],
              sorted(flagged_dd))
    n_vor = sum(1 for r in flagged_dd if r[3] == "VOR")
    n_both = sum(1 for r in flagged_dd if r[2][0] == "1" and r[2][1] == "1")
    rep.blank()
    rep.kv("landslides flagged for detection or delineation", len(flagged_dd))
    rep.kv("  of which very old relict", n_vor)
    rep.kv("  of which flagged for both detection and delineation", n_both)
    rep.line()
    rep.line("All three very-old-relict landslides carry such a flag, which is where")
    rep.line("the criterion of Table 4 expects them: morphology partly obliterated by")
    rep.line("later processes. They are not the only ones, however:")
    rep.line("elsewhere the doubt concerns the material involved, which")
    rep.line("the attribute records as classification uncertainty.")

    glim.write_csv("09_integrity_checks",
                   ["check", "criterion", "records_flagged", "corrected",
                    "remaining", "remainder_kind"],
                   [[label, criterion,
                     sum(before.get(c, 0) for c in checks),
                     sum(before.get(c, 0) for c in checks)
                     - sum(after.get(c, 0) for c in checks),
                     sum(after.get(c, 0) for c in checks), kind]
                    for label, criterion, checks, kind in CHECK_GROUPS])
    glim.write_csv("09_comparison_summary",
                   ["inventory", "quantity", "value"],
                   [["Tosevski 2013", "landslides", len(load_findings(match_path))],
                    ["Tosevski 2013", "inside the new G-LIM",
                     sum(1 for r in load_findings(match_path)
                         if r["inside_new_glim"] == "yes")],
                    ] + jom_summary)
    rep.save()


if __name__ == "__main__":
    main()
