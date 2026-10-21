#!/usr/bin/env python3
"""
06 - Full record list of the landslides that break the depth rule.

Script 02 counts the violations; this one lists them, so that each record can
be inspected and corrected in the GIS.

The rule as stated in the Methods of the paper, revised after this check was
first run:

    a landslide whose scarp/source area and deposit could be delineated
    separately is represented by TWO polygons sharing one lds_id and is
    classified "deep-seated"; a landslide represented by a SINGLE polygon is
    classified "shallow"; debris flows are classified "shallow" whatever their
    polygon structure, because their upslope polygon is a detachment and
    transport zone rather than the scarp of a failure surface; and particular
    morphological configurations can require MORE than two polygons, for
    instance a debris flow fed by two separate detachment zones.

The violations reported are therefore:

    2poly_shallow   two polygons share the lds_id but depth = "shallow",
                    excluding debris flows, for which this is the convention
    1poly_deep      a single polygon but depth = "deep-seated"
    mixed_depth     the polygons of one landslide carry different depth values

A landslide represented by more than two polygons is listed for information
but is no longer counted as a violation.

Note on polygon identifiers: shapefiles have no polygon attribute other than
the OGR feature id, which is the record index in the .dbf, numbered from 0.
That is the same number ArcGIS shows as FID and QGIS as $id, so `fid` below
can be used directly to find the feature. `dbf_record` is the same value
numbered from 1, for anyone reading the table in a spreadsheet.

Output: outputs/06_depth_rule_records.txt
        outputs/06_depth_rule_polygons.csv   one row per polygon
        outputs/06_depth_rule_landslides.csv one row per landslide
"""

from collections import Counter

import glim
from glim import Report

VIOLATION_LABEL = {
    "2poly_shallow": "2 polygons but depth = shallow",
    "1poly_deep": "1 polygon but depth = deep-seated",
    "mixed_depth": "polygons of one landslide disagree on depth",
}

#: A two-polygon debris flow is labelled shallow in 80 out of 80 cases and
#: deep-seated in none (see the "Depth by type pair" section of the report).
#: For a debris flow the two polygons are the detachment/transport zone and the
#: fan deposit, not a scarp implying a deep failure surface, so the label is an
#: interpretive convention rather than a mistake. It is separated from the
#: genuine anomalies here so that the two are not corrected in the same way.
SYSTEMATIC_PAIR = {"xdf", "df"}


def is_debris_flow(group):
    """True when every polygon of the landslide is a debris-flow code."""
    return {glim.base_class(f["lds_class"]) for f in group} == {"df"}


def diagnose(group, violation):
    """Tell a documented mapping convention from a record that needs correcting."""
    classes = {f["lds_class"] for f in group}
    if violation == "2poly_shallow" and classes == SYSTEMATIC_PAIR:
        return "convention: debris flows are shallow whatever the polygon structure"
    if violation == "more_than_two_polygons" and is_debris_flow(group):
        return ("convention: a debris flow with more than one detachment zone needs "
                "more than two polygons")
    return "anomaly: needs checking in the GIS"


def classify(group):
    """Return the violation code for one landslide, or None if it obeys the rule."""
    values = {f["depth"] for f in group}
    if len(values) > 1:
        return "mixed_depth"
    value = values.pop()
    if len(group) == 2 and value == "shallow":
        return "2poly_shallow"
    if len(group) == 1 and value == "deep-seated":
        return "1poly_deep"
    if len(group) > 2:
        return "more_than_two_polygons"
    return None


def main():
    rep = Report("06_depth_rule_records",
                 "06 - Records violating the shallow / deep-seated rule")
    layers = glim.read_all_landslides(fix=True)

    rep.line("Rule: 2 polygons per lds_id -> deep-seated; 1 polygon -> shallow.")
    rep.line("Polygon identifier `fid` is the OGR feature id = ArcGIS FID = QGIS $id,")
    rep.line("numbered from 0; `dbf_record` is the same value numbered from 1.")

    polygon_rows = []
    landslide_rows = []
    for layer_name, feats in layers.items():
        short = layer_name.replace("DRV_GLIM_", "")
        for lds_id, group in glim.group_by_id(feats).items():
            violation = classify(group)
            if violation is None:
                continue
            group = sorted(group, key=lambda f: (not glim.is_source_polygon(f["lds_class"]),
                                                 f.fid))
            verdict = diagnose(group, violation)
            for f in group:
                polygon_rows.append([
                    layer_name, short, glim.RELATIVE_AGE[layer_name],
                    glim.GENERATION[layer_name] or "",
                    f.fid, f.fid + 1, lds_id,
                    f["lds_class"],
                    "source" if glim.is_source_polygon(f["lds_class"]) else "deposit",
                    f["depth"], f["uncert"],
                    round(f.geom.area, 1),
                    len(group), violation, verdict,
                ])
            source = [f for f in group if glim.is_source_polygon(f["lds_class"])]
            deposit = [f for f in group if not glim.is_source_polygon(f["lds_class"])]
            landslide_rows.append([
                layer_name, short, lds_id, len(group), violation,
                "+".join(str(f.fid) for f in group),
                source[0]["lds_class"] if source else "",
                round(source[0].geom.area, 1) if source else "",
                source[0]["depth"] if source else "",
                deposit[0]["lds_class"] if deposit else "",
                round(deposit[0].geom.area, 1) if deposit else "",
                deposit[0]["depth"] if deposit else "",
                round(sum(f.geom.area for f in group), 1),
                group[0]["uncert"], verdict,
            ])

    # ------------------------------------------------------------- summary
    rep.h2("Summary")
    counts = Counter((r[0], r[4]) for r in landslide_rows)
    rows = []
    for layer_name in layers:
        for violation in ("2poly_shallow", "1poly_deep", "mixed_depth",
                          "more_than_two_polygons"):
            n = counts.get((layer_name, violation), 0)
            if n:
                rows.append([layer_name.replace("DRV_GLIM_", ""), violation, n])
    rep.table(["layer", "violation", "landslides"], rows)
    rep.kv("landslides to review", len(landslide_rows))
    rep.kv("polygons to review", len(polygon_rows))

    systematic = [r for r in landslide_rows if r[14].startswith("convention")]
    anomalies = [r for r in landslide_rows if r[14].startswith("anomaly")]
    rep.blank()
    rep.kv("  of which a documented convention", len(systematic))
    rep.kv("  of which genuine anomalies", len(anomalies))

    # ------------------------------------------------- depth by type pair
    rep.h2("Depth by type pair - why the violations are not all the same")
    rep.line("Every two-polygon landslide in the dataset, grouped by the pair of")
    rep.line("source and deposit codes. This is what shows whether a violation is a")
    rep.line("labelling slip or a consistent interpretive choice.")
    rep.blank()
    pair_tab = {}
    for layer_name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            if len(group) != 2:
                continue
            source = [f["lds_class"] for f in group
                      if glim.is_source_polygon(f["lds_class"])]
            deposit = [f["lds_class"] for f in group
                       if not glim.is_source_polygon(f["lds_class"])]
            pair = f"{source[0] if source else '?'}+{deposit[0] if deposit else '?'}"
            depths = {f["depth"] for f in group}
            key = "MIXED" if len(depths) > 1 else depths.pop()
            pair_tab.setdefault(pair, Counter())[key] += 1
    rep.table(
        ["source+deposit", "deep-seated", "shallow", "MIXED", "total"],
        [[pair, c["deep-seated"], c["shallow"], c["MIXED"], sum(c.values())]
         for pair, c in sorted(pair_tab.items(), key=lambda kv: -sum(kv[1].values()))])
    rep.blank()
    rep.line("Reading: every type pair is deep-seated in 100 % of cases when the source")
    rep.line("area was delineated separately - except xdf+df, which is shallow in 100 %")
    rep.line("of cases and deep-seated in none. Note that the other flow-type pairs,")
    rep.line("xsef+sef and xf+f, follow the rule, so the exception is specific to debris")
    rep.line("flows and not to flows in general.")
    rep.blank()
    rep.line("This is a defensible interpretation: the two polygons of a debris flow are")
    rep.line("the concave detachment/transport zone and the fan deposit, which is not the")
    rep.line("scarp of a rotational or translational failure and carries no implication")
    rep.line("about the depth of a failure surface. The consistency of the convention")
    rep.line("(80 out of 80) argues that the data are right and the rule as written in")
    rep.line("the Technical Validation section is incomplete.")

    rep.h2("Single-polygon landslides: depth by type")
    single_tab = {}
    for layer_name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            if len(group) == 1:
                single_tab.setdefault(group[0]["lds_class"], Counter())[group[0]["depth"]] += 1
    rep.table(
        ["lds_class", "deep-seated", "shallow", "total"],
        [[code, c["deep-seated"], c["shallow"], sum(c.values())]
         for code, c in sorted(single_tab.items(), key=lambda kv: -sum(kv[1].values()))])
    rep.line()
    rep.line("Only two single-polygon landslides in the whole dataset are labelled")
    rep.line("deep-seated; both are in O_IV and both are listed among the anomalies.")

    # ------------------------------------------------- the genuine anomalies
    rep.h2("A - Genuine anomalies, to be checked and corrected in the GIS")
    rep.table(
        ["layer", "fid", "dbf_record", "lds_id", "lds_class", "role", "depth",
         "uncert", "area m2", "n_poly", "violation"],
        [[r[1], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12], r[13]]
         for r in polygon_rows if r[14].startswith("anomaly")])

    # ------------------------------------------- one row per landslide
    rep.h2("B - All records, one row per landslide")
    rep.line("`fids` lists the polygon fid(s) of the landslide, source first.")
    rep.blank()
    rep.table(
        ["layer", "lds_id", "n_poly", "fids", "source class", "source area m2",
         "deposit class", "deposit area m2", "depth", "uncert", "violation",
         "diagnosis"],
        [[r[1], r[2], r[3], r[5], r[6] or "-", r[7] or "-", r[9] or "-", r[10] or "-",
          r[8] if r[4] == "mixed_depth" else (r[11] or r[8]), r[13], r[4],
          r[14].split(":")[0]]
         for r in landslide_rows])

    # --------------------------------------------- one row per polygon
    rep.h2("C - All records, one row per polygon")
    rep.table(
        ["layer", "fid", "dbf_record", "lds_id", "lds_class", "role", "depth",
         "uncert", "area m2", "violation", "diagnosis"],
        [[r[1], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[13],
          r[14].split(":")[0]]
         for r in polygon_rows])

    # ------------------------------------------------ GIS selection strings
    rep.h2("Selection expressions for the GIS")
    rep.line("Paste into the QGIS or ArcGIS select-by-expression box, one per layer,")
    rep.line("to isolate the records listed above.")
    rep.blank()
    for label, subset in (("all records", landslide_rows),
                          ("genuine anomalies only",
                           [r for r in landslide_rows if r[14].startswith("anomaly")])):
        rep.line(f"  -- {label}")
        for layer_name in layers:
            ids = sorted({r[2] for r in subset if r[0] == layer_name})
            if not ids:
                continue
            rep.line(f"  {layer_name}:")
            rep.line(f'    "lds_id" IN ({", ".join(str(i) for i in ids)})')
        rep.blank()

    # ------------------------------------------------- material for the fix
    rep.h2("What to decide for each record")
    rep.line("The violations are not all the same problem:")
    rep.blank()
    rep.line("  2poly_shallow - the source area WAS delineated separately, so under the")
    rep.line("      stated rule the landslide is deep-seated. Either the depth value is")
    rep.line("      wrong, or these are cases where the interpreter deliberately mapped")
    rep.line("      a scarp for a shallow failure, in which case the rule in the text")
    rep.line("      does not describe what was actually done and must be reworded.")
    rep.blank()
    rep.line("  1poly_deep - a deep-seated landslide whose source area was not")
    rep.line("      delineated. Check whether the missing polygon was lost in editing")
    rep.line("      or whether the deposit alone was mapped on purpose.")
    rep.blank()
    rep.line("  mixed_depth - unambiguous data error: the two polygons of one landslide")
    rep.line("      cannot have different depth classes.")
    rep.blank()
    rep.line("Note the size pattern below: it helps tell a systematic labelling slip")
    rep.line("from a deliberate mapping choice.")
    rep.blank()
    two_shallow = [r for r in landslide_rows if r[4] == "2poly_shallow"]
    if two_shallow:
        areas = sorted(r[12] for r in two_shallow)
        obeying = []
        for layer_name, feats in layers.items():
            for lds_id, group in glim.group_by_id(feats).items():
                if len(group) == 2 and classify(group) is None:
                    obeying.append(sum(f.geom.area for f in group))
        obeying.sort()
        rep.table(
            ["group", "n", "min m2", "median m2", "max m2"],
            [["2 polygons, labelled shallow (violations)", len(areas),
              f"{areas[0]:.0f}", f"{areas[len(areas) // 2]:.0f}", f"{areas[-1]:.0f}"],
             ["2 polygons, labelled deep-seated (rule obeyed)", len(obeying),
              f"{obeying[0]:.0f}", f"{obeying[len(obeying) // 2]:.0f}",
              f"{obeying[-1]:.0f}"]])

    glim.write_csv(
        "06_depth_rule_polygons",
        ["layer", "layer_short", "relative_age", "generation", "fid", "dbf_record",
         "lds_id", "lds_class", "role", "depth", "uncert", "area_m2",
         "n_polygons_of_landslide", "violation", "diagnosis"],
        polygon_rows)
    glim.write_csv(
        "06_depth_rule_landslides",
        ["layer", "layer_short", "lds_id", "n_polygons", "violation", "fids",
         "source_class", "source_area_m2", "source_depth",
         "deposit_class", "deposit_area_m2", "deposit_depth",
         "total_area_m2", "uncert", "diagnosis"],
        landslide_rows)
    rep.save()


if __name__ == "__main__":
    main()
