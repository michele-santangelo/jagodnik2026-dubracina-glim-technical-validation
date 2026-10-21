#!/usr/bin/env python3
"""
07 - Geometric and topological quality control, at record level.

Script 02 checks the database side (attributes, dictionaries, identifiers) and
the coarse geometry flags. This one goes after the geometry itself, and reports
every offending record with its fid so that it can be found and fixed in the
GIS one by one.

Checks
  A  validity, with the coordinate of each defect
  B  exact and near-duplicate geometries, within and across layers
  C  degenerate and sliver polygons (tiny area, extreme shape index, few vertices)
  D  vertex-level defects: repeated vertices, sub-decimetre segments, spikes
  E  overlaps between different landslides of the SAME generation layer
  F  overlap between the source area and the deposit of the same landslide
  G  source area contained in its deposit, or the reverse
  H  interior rings (holes)
  I  overlaps inside the geomorphological layer
  J  the study-area boundary

Geometries are read unrepaired throughout.

Output: outputs/07_qc_topology.txt
        outputs/07_topology_findings.csv   one row per offending record
"""

import math
import os
from collections import Counter, defaultdict

from shapely.geometry import Polygon
from shapely.strtree import STRtree
from shapely.validation import explain_validity

import glim
from glim import Report

#: two polygons are treated as intentionally adjacent, not overlapping, when the
#: shared area is below this. Digitising at 1:10,000 from a 1 m DTM leaves
#: slivers of a fraction of a square metre along shared boundaries.
OVERLAP_TOL_M2 = 1.0
#: segments shorter than this are redundant nodes for a 1 m DTM at 1:10,000
SHORT_SEGMENT_M = 0.10
#: interior angle below which a vertex is reported as a spike
SPIKE_DEG = 3.0
#: a gap narrower than this between two polygons is a floating-point artefact of
#: snapping, not a mapped gap: the observed cases are of the order of 1e-10 m
SNAP_TOL_M = 0.01
#: polygons smaller than this are degenerate rather than small
DEGENERATE_M2 = 10.0
#: shape index above which a polygon is reported as a sliver.
#: index = perimeter / (2 * sqrt(pi * area)); a circle is 1.0, a 1:100 sliver ~5.7
SLIVER_INDEX = 6.0
#: near-duplicate threshold on intersection over union
NEAR_DUPLICATE_IOU = 0.98

findings = []


def record(check, layer, feature, description, severity="check"):
    """Log one finding, keyed by a stable attribute key plus a position.

    `feature` may be a glim.Feature or, where no feature applies, a string.
    The fid is recorded too but only as a hint: it is positional and changes
    as soon as a record is deleted.
    """
    if isinstance(feature, str):
        key, fid, cx, cy = feature, "-", "", ""
    else:
        key = glim.stable_key(feature)
        fid = feature.fid
        cx, cy = glim.centroid_xy(feature)
    findings.append([check, layer, key, severity, description, cx, cy, fid])


def shape_index(geom):
    if geom.area <= 0:
        return float("inf")
    return geom.length / (2.0 * math.sqrt(math.pi * geom.area))


def rings(geom):
    for part in glim.parts(geom):
        if isinstance(part, Polygon):
            yield part.exterior
            for hole in part.interiors:
                yield hole


def vertex_defects(geom):
    """Return (repeated vertices, short segments, spikes) for one geometry."""
    repeated = short = spikes = 0
    for ring in rings(geom):
        coords = list(ring.coords)[:-1]  # drop the closing duplicate
        n = len(coords)
        if n < 3:
            continue
        for i in range(n):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % n]
            length = math.hypot(bx - ax, by - ay)
            if length == 0.0:
                repeated += 1
                continue
            if length < SHORT_SEGMENT_M:
                short += 1
            cx, cy = coords[(i + 2) % n]
            length2 = math.hypot(cx - bx, cy - by)
            if length2 == 0.0:
                continue
            # angle at vertex b between segments b->a and b->c
            dot = ((ax - bx) * (cx - bx) + (ay - by) * (cy - by)) / (length * length2)
            angle = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
            if angle < SPIKE_DEG:
                spikes += 1
    return repeated, short, spikes


def n_vertices(geom):
    return sum(len(list(r.coords)) - 1 for r in rings(geom))


def main():
    rep = Report("07_qc_topology", "07 - Geometric and topological quality control")
    layers = glim.read_all_landslides(fix=False)
    rockfalls = glim.read_rockfalls(fix=False)
    geomorph = glim.read_layer(glim.ORIGINAL_DATA, glim.GEOMORPHOLOGY_LAYER,
                               ["gm_id", "geom_class", "Shape_Area"], fix=False)
    polygon_sets = list(layers.items()) + [(glim.ROCKFALL_LAYER, rockfalls),
                                           (glim.GEOMORPHOLOGY_LAYER, geomorph)]

    def ident_of(feature):
        return feature.attrs.get("lds_id", feature.attrs.get("gm_id"))

    def short(name):
        return name.replace("DRV_GLIM_", "")

    rep.line("Every offending record is keyed by a STABLE attribute identifier, because")
    rep.line("shapefile feature ids are positional: deleting one record renumbers every")
    rep.line("record after it. The keys are (lds_id/lds_class) for the landslide layers")
    rep.line("and gm_id for the geomorphological layer. The rockfall layer has no usable")
    rep.line("key, since Id is 0 for all 486 records, so it is located by centroid.")
    rep.line("Columns marked fid* are positional hints only, valid for the current files.")
    rep.line(f"Tolerances: overlap > {OVERLAP_TOL_M2} m2, segment < {SHORT_SEGMENT_M} m,")
    rep.line(f"spike angle < {SPIKE_DEG} deg, degenerate area < {DEGENERATE_M2} m2,")
    rep.line(f"sliver shape index > {SLIVER_INDEX}, near-duplicate IoU > {NEAR_DUPLICATE_IOU}.")

    # ------------------------------------------------------------ A validity
    rep.h2("A - Invalid geometries")
    rows = []
    for name, feats in polygon_sets:
        for f in feats:
            if not f.geom.is_valid:
                reason = explain_validity(f.geom)
                rows.append([short(name), glim.stable_key(f),
                             f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                             f"{f.geom.area:.1f}", f.fid, reason])
                record("A_invalid", name, f, reason, "fix")
    rep.table(["layer", "key", "class", "area m2", "fid*", "reason"],
              rows or [["-", "-", "-", "-", "-", "none"]])
    if rows:
        rep.line()
        rep.line("A ring self-intersection is a figure-of-eight in the outline. In QGIS,")
        rep.line("Vector > Geometry Tools > Check Validity locates the node; the fix is to")
        rep.line("move or delete the crossing vertex, not to run a global buffer(0), which")
        rep.line("would silently reshape the polygon.")

    # -------------------------------------------------------- B duplicates
    rep.h2("B - Exact duplicate geometries")
    rows = []
    for name, feats in polygon_sets:
        seen = defaultdict(list)
        for f in feats:
            key = (f.geom if f.geom.is_valid else f.geom.buffer(0)).wkb
            seen[key].append(f)
        for group in seen.values():
            if len(group) > 1:
                for f in group:
                    rows.append([short(name), glim.stable_key(f),
                                 f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                                 f"{f.geom.area:.1f}", f.fid,
                                 ", ".join(glim.stable_key(g) for g in group
                                           if g.fid != f.fid)])
                    record("B_duplicate", name, f, "identical geometry to " +
                           ", ".join(glim.stable_key(g) for g in group
                                     if g.fid != f.fid), "fix")
    rep.table(["layer", "key", "class", "area m2", "fid*", "identical to"],
              rows or [["-", "-", "-", "-", "-", "none"]])

    rep.h2("B2 - Near-duplicate geometries within a layer")
    rep.line(f"Pairs with IoU > {NEAR_DUPLICATE_IOU}: the same landform digitised twice.")
    rep.blank()
    rows = []
    for name, feats in polygon_sets:
        geoms = [f.geom if f.geom.is_valid else f.geom.buffer(0) for f in feats]
        tree = STRtree(geoms)
        done = set()
        for i, f in enumerate(feats):
            for j in tree.query(geoms[i]):
                if j <= i or (i, j) in done:
                    continue
                done.add((i, j))
                value = glim.iou(geoms[i], geoms[j])
                if value > NEAR_DUPLICATE_IOU:
                    rows.append([short(name), glim.stable_key(f),
                                 glim.stable_key(feats[j]),
                                 f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                                 feats[j].attrs.get("lds_class",
                                                    feats[j].attrs.get("geom_class", "")),
                                 f"{value:.3f}", f"{geoms[i].area:.1f}",
                                 f"{f.fid}/{feats[j].fid}"])
                    record("B2_near_duplicate", name, f,
                           f"IoU {value:.3f} with {glim.stable_key(feats[j])}", "check")
    rep.table(["layer", "key A", "key B", "class A", "class B", "IoU", "area m2",
               "fid* A/B"], rows or [["-"] * 8])

    rep.h2("B3 - Identical geometries across different layers")
    rep.line("A landslide copied from one generation layer into another.")
    rep.blank()
    across = defaultdict(list)
    for name, feats in polygon_sets:
        for f in feats:
            key = (f.geom if f.geom.is_valid else f.geom.buffer(0)).wkb
            across[key].append((name, f))
    rows = []
    for group in across.values():
        names = {n for n, _ in group}
        if len(group) > 1 and len(names) > 1:
            for name, f in group:
                rows.append([short(name), glim.stable_key(f),
                             f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                             f"{f.geom.area:.1f}", f.fid,
                             ", ".join(f"{short(n)} {glim.stable_key(g)}"
                                       for n, g in group
                                       if not (g.fid == f.fid and n == name))])
                record("B3_duplicate_across_layers", name, f,
                       "identical geometry present in " +
                       ", ".join(f"{short(n)} {glim.stable_key(g)}" for n, g in group
                                 if not (g.fid == f.fid and n == name)), "fix")
    rep.table(["layer", "key", "class", "area m2", "fid*", "also in"],
              rows or [["-", "-", "-", "-", "-", "none"]])

    # ------------------------------------------------- C degenerate / slivers
    rep.h2("C - Degenerate polygons")
    rep.line(f"Area below {DEGENERATE_M2} m2, or fewer than four distinct vertices.")
    rep.blank()
    rows = []
    for name, feats in polygon_sets:
        for f in feats:
            nv = n_vertices(f.geom)
            if f.geom.area < DEGENERATE_M2 or nv < 4:
                cx, cy = glim.centroid_xy(f)
                rows.append([short(name), glim.stable_key(f),
                             f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                             f"{f.geom.area:.3f}", nv,
                             f"{f.attrs.get('Shape_Area', 0) or 0:.1f}",
                             f"{cx} {cy}", f.fid])
                record("C_degenerate", name, f,
                       f"area {f.geom.area:.3f} m2, {nv} vertices, "
                       f"Shape_Area attribute {f.attrs.get('Shape_Area', 0) or 0:.1f}",
                       "fix")
    rep.table(["layer", "key", "class", "area m2", "vertices",
               "Shape_Area attr", "centroid E N", "fid*"], rows or [["-"] * 8])
    if rows:
        rep.line()
        rep.line("Every one of these is in the geomorphological layer, all but one is talus,")
        rep.line("and each carries a Shape_Area of thousands of square metres against a")
        rep.line("geometry of a fraction of one. That value is shared with a group of larger")
        rep.line("polygons: one original unit was subdivided and the field was never")
        rep.line("recomputed, so no polygon in the dataset actually has the area the field")
        rep.line("reports. Section J of report 02 lists the groups. Do not look for a")
        rep.line("polygon of 178,967 m2 - none exists; that number is a leftover carried by")
        rep.line("twelve talus polygons whose geometries sum to 52,844 m2.")

    rep.h2("C2 - Sliver polygons")
    rep.line(f"Shape index above {SLIVER_INDEX} (a circle is 1.0). Long, narrow polygons")
    rep.line("are expected for debris-flow channels and talus strips, so this list needs")
    rep.line("visual triage rather than blanket correction.")
    rep.blank()
    rows = []
    for name, feats in polygon_sets:
        for f in feats:
            geom = f.geom if f.geom.is_valid else f.geom.buffer(0)
            index = shape_index(geom)
            if index > SLIVER_INDEX:
                rows.append([short(name), glim.stable_key(f),
                             f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                             f"{geom.area:.1f}", f"{geom.length:.1f}", f"{index:.1f}",
                             f.fid])
                record("C2_sliver", name, f,
                       f"shape index {index:.1f}, area {geom.area:.1f} m2", "check")
    rows.sort(key=lambda r: -float(r[5]))
    rep.table(["layer", "key", "class", "area m2", "perimeter m", "shape index", "fid*"],
              rows[:40] or [["-"] * 7])
    if len(rows) > 40:
        rep.line(f"    ... {len(rows) - 40} more in the CSV, sorted by shape index.")

    # ------------------------------------------------------ D vertex defects
    rep.h2("D - Vertex-level defects")
    rows = []
    totals = Counter()
    for name, feats in polygon_sets:
        for f in feats:
            repeated, short_seg, spikes = vertex_defects(f.geom)
            if repeated or spikes:
                rows.append([short(name), glim.stable_key(f),
                             f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                             repeated, short_seg, spikes, n_vertices(f.geom), f.fid])
                record("D_vertex_defect", name, f,
                       f"{repeated} repeated vertices, {short_seg} short segments, "
                       f"{spikes} spikes", "fix" if repeated else "check")
            totals["repeated"] += repeated
            totals["short"] += short_seg
            totals["spikes"] += spikes
            totals["vertices"] += n_vertices(f.geom)
    rep.kv("total vertices in the polygon layers", totals["vertices"])
    rep.kv("repeated (zero-length) vertices", totals["repeated"])
    rep.kv(f"segments shorter than {SHORT_SEGMENT_M} m", totals["short"])
    rep.kv(f"spikes sharper than {SPIKE_DEG} deg", totals["spikes"])
    rep.blank()
    rows.sort(key=lambda r: (-r[3], -r[5]))
    rep.table(["layer", "key", "class", "repeated", "short seg", "spikes",
               "vertices", "fid*"], rows[:40] or [["-"] * 8])
    if len(rows) > 40:
        rep.line(f"    ... {len(rows) - 40} more in the CSV.")

    # -------------------------------------- E overlaps within one generation
    rep.h2("E - Overlaps between different landslides of the same generation")
    rep.line("Two landslides in one generation layer should be adjacent at most, never")
    rep.line("superimposed: a superposition is what defines a different generation.")
    rep.line("Overlaps here are therefore either digitising errors or landslides that")
    rep.line("belong in different generation layers.")
    rep.blank()
    rows = []
    for name, feats in layers.items():
        geoms = [f.geom if f.geom.is_valid else f.geom.buffer(0) for f in feats]
        tree = STRtree(geoms)
        seen = set()
        for i, f in enumerate(feats):
            for j in tree.query(geoms[i]):
                if j <= i:
                    continue
                if feats[j]["lds_id"] == f["lds_id"]:
                    continue  # source/deposit of one landslide: check F
                key = (i, j)
                if key in seen:
                    continue
                seen.add(key)
                area = geoms[i].intersection(geoms[j]).area
                if area > OVERLAP_TOL_M2:
                    smaller = min(geoms[i].area, geoms[j].area)
                    rows.append([short(name), glim.stable_key(f),
                                 glim.stable_key(feats[j]),
                                 f"{area:.1f}", f"{100 * area / smaller:.1f}",
                                 f"{f.fid}/{feats[j].fid}"])
                    record("E_intra_layer_overlap", name, f,
                           f"overlaps {glim.stable_key(feats[j])} by {area:.1f} m2",
                           "check")
    rows.sort(key=lambda r: -float(r[3]))
    rep.table(["layer", "key A", "key B", "overlap m2", "% of smaller", "fid* A/B"],
              rows[:40] or [["-"] * 6])
    if len(rows) > 40:
        rep.line(f"    ... {len(rows) - 40} more in the CSV, sorted by overlap area.")
    rep.blank()
    rep.kv("overlapping pairs within one generation layer", len(rows))

    # -------------------------------------- F source / deposit relationship
    rep.h2("F - Source area and deposit of the same landslide")
    rep.line("The two polygons should share a boundary: no overlap, no gap. An overlap")
    rep.line("double counts the area; a gap leaves a strip of unmapped ground between")
    rep.line("the scarp and the deposit of one landslide.")
    rep.blank()
    overlap_rows = []
    gap_rows = []
    snap_rows = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            if len(group) != 2:
                continue
            a, b = group
            ga = a.geom if a.geom.is_valid else a.geom.buffer(0)
            gb = b.geom if b.geom.is_valid else b.geom.buffer(0)
            inter = ga.intersection(gb).area
            if inter > OVERLAP_TOL_M2:
                overlap_rows.append([short(name), lds_id,
                                     a["lds_class"], b["lds_class"], f"{inter:.1f}",
                                     f"{100 * inter / min(ga.area, gb.area):.1f}",
                                     f"{a.fid}/{b.fid}"])
                record("F_source_deposit_overlap", name, a,
                       f"overlaps {glim.stable_key(b)}, the other polygon of the same "
                       f"landslide, by {inter:.1f} m2", "fix")
            elif not ga.intersects(gb):
                distance = ga.distance(gb)
                if distance < SNAP_TOL_M:
                    snap_rows.append([short(name), lds_id,
                                      a["lds_class"], b["lds_class"],
                                      f"{distance:.2e}", f"{a.fid}/{b.fid}"])
                    record("F_snap_artefact", name, a,
                           f"separated from {glim.stable_key(b)} by {distance:.2e} m, "
                           "a snapping artefact rather than a mapped gap", "cosmetic")
                else:
                    gap_rows.append([short(name), lds_id,
                                     a["lds_class"], b["lds_class"], f"{distance:.2f}",
                                     f"{a.fid}/{b.fid}"])
                    record("F_source_deposit_gap", name, a,
                           f"does not touch {glim.stable_key(b)}, the other polygon of "
                           f"the same landslide: gap {distance:.2f} m", "check")
    rep.line("F1 - overlapping source and deposit")
    overlap_rows.sort(key=lambda r: -float(r[4]))
    rep.table(["layer", "lds_id", "class A", "class B", "overlap m2", "% of smaller",
               "fid* A/B"], overlap_rows[:30] or [["-"] * 7])
    if len(overlap_rows) > 30:
        rep.line(f"    ... {len(overlap_rows) - 30} more in the CSV.")
    rep.blank()
    rep.line("F2 - a real gap between source and deposit")
    gap_rows.sort(key=lambda r: -float(r[4]))
    rep.table(["layer", "lds_id", "class A", "class B", "gap m", "fid* A/B"],
              gap_rows[:30] or [["-"] * 6])
    if len(gap_rows) > 30:
        rep.line(f"    ... {len(gap_rows) - 30} more in the CSV, sorted by gap width.")
    rep.blank()
    rep.line(f"F3 - separated by less than {SNAP_TOL_M} m: a snapping artefact, not a gap")
    rep.table(["layer", "lds_id", "class A", "class B", "gap m", "fid* A/B"],
              snap_rows or [["-"] * 6])
    rep.line()
    rep.line("These four are of the order of 1e-10 m: the shared boundary was digitised")
    rep.line("twice and the two copies differ in the last bits of the coordinate. Nothing")
    rep.line("to correct by hand; a snap-to-vertex pass at 0.001 m would remove them.")
    rep.blank()
    rep.kv("landslides with overlapping source and deposit", len(overlap_rows))
    rep.kv("landslides with a real gap between the two polygons", len(gap_rows))
    rep.kv("landslides separated only by a snapping artefact", len(snap_rows))

    # ----------------------------------------------------- G containment
    rep.h2("G - One polygon of a landslide contained in the other")
    rows = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            if len(group) != 2:
                continue
            a, b = group
            ga = a.geom if a.geom.is_valid else a.geom.buffer(0)
            gb = b.geom if b.geom.is_valid else b.geom.buffer(0)
            for inner, outer, fi, fo in ((ga, gb, a, b), (gb, ga, b, a)):
                if outer.buffer(0.01).contains(inner):
                    rows.append([short(name), lds_id,
                                 fi["lds_class"], fo["lds_class"],
                                 f"{inner.area:.1f}", f"{outer.area:.1f}",
                                 f"{fi.fid}/{fo.fid}"])
                    record("G_containment", name, fi,
                           f"contained in {glim.stable_key(fo)} of the same landslide",
                           "fix")
    rep.table(["layer", "lds_id", "inner class", "outer class", "inner area m2",
               "outer area m2", "fid* inner/outer"], rows or [["-"] * 7])

    # ---------------------------------------------------------- H holes
    rep.h2("H - Interior rings (holes)")
    rep.line("A hole in a landslide polygon means an area inside the landslide that was")
    rep.line("deliberately excluded. Confirm each one is intentional and not a")
    rep.line("digitising artefact.")
    rep.blank()
    rows = []
    for name, feats in polygon_sets:
        for f in feats:
            holes = [h for part in glim.parts(f.geom) if isinstance(part, Polygon)
                     for h in part.interiors]
            if holes:
                areas = sorted(Polygon(h).area for h in holes)
                cx, cy = glim.centroid_xy(f)
                rows.append([short(name), glim.stable_key(f),
                             f.attrs.get("lds_class", f.attrs.get("geom_class", "")),
                             len(holes), f"{areas[0]:.1f}", f"{areas[-1]:.1f}",
                             f"{sum(areas):.1f}", f"{cx} {cy}"])
                record("H_holes", name, f,
                       f"{len(holes)} interior ring(s), total {sum(areas):.1f} m2",
                       "check")
    rows.sort(key=lambda r: -r[3])
    rep.table(["layer", "key", "class", "holes", "smallest m2", "largest m2",
               "total m2", "centroid E N"], rows or [["-"] * 8])

    # ------------------------------------------- I geomorphology overlaps
    rep.h2("I - Overlaps inside the geomorphological layer")
    rep.line("Surface-deposit units of different classes describe different materials at")
    rep.line("the same place, so overlaps need a decision: either the units are nested by")
    rep.line("design, and this must be documented, or the boundaries are wrong.")
    rep.blank()
    geoms = [f.geom if f.geom.is_valid else f.geom.buffer(0) for f in geomorph]
    tree = STRtree(geoms)
    rows = []
    seen = set()
    for i, f in enumerate(geomorph):
        for j in tree.query(geoms[i]):
            if j <= i or (i, j) in seen:
                continue
            seen.add((i, j))
            area = geoms[i].intersection(geoms[j]).area
            if area > OVERLAP_TOL_M2:
                rows.append([f["gm_id"], geomorph[j]["gm_id"],
                             f["geom_class"], geomorph[j]["geom_class"], f"{area:.1f}",
                             f"{100 * area / min(geoms[i].area, geoms[j].area):.1f}",
                             f"{f.fid}/{geomorph[j].fid}"])
                record("I_geomorph_overlap", glim.GEOMORPHOLOGY_LAYER, f,
                       f"overlaps gm_id {geomorph[j]['gm_id']} "
                       f"({geomorph[j]['geom_class']}) by {area:.1f} m2", "check")
    rows.sort(key=lambda r: -float(r[4]))
    rep.table(["gm_id A", "gm_id B", "class A", "class B", "overlap m2",
               "% of smaller", "fid* A/B"], rows[:30] or [["-"] * 7])
    if len(rows) > 30:
        rep.line(f"    ... {len(rows) - 30} more in the CSV.")
    rep.blank()
    rep.kv("overlapping pairs of geomorphological units", len(rows))
    pair_counts = Counter()
    for r in rows:
        pair_counts[tuple(sorted((r[2], r[3])))] += 1
    if pair_counts:
        rep.line("    by class pair: " +
                 ", ".join(f"{a}/{b}={n}" for (a, b), n in pair_counts.most_common()))

    # ------------------------------------------------------- J boundary
    rep.h2("J - Study-area boundary")
    boundary_feats = glim.read_layer(glim.ORIGINAL_DATA, glim.BOUNDARY_LAYER,
                                     ["Id", "area_km2"], fix=False)
    b = boundary_feats[0]
    repeated, short_seg, spikes = vertex_defects(b.geom)
    rep.table(["property", "value"],
              [["features", len(boundary_feats)],
               ["valid", "yes" if b.geom.is_valid else explain_validity(b.geom)],
               ["geometry type", b.geom.geom_type],
               ["vertices", n_vertices(b.geom)],
               ["interior rings", len(b.geom.interiors)
                if isinstance(b.geom, Polygon) else "n/a"],
               ["grid area km2", f"{glim.km2(b.geom.area):.4f}"],
               ["ellipsoidal area km2",
                f"{glim.km2(glim.geodesic_area(b.geom)):.4f}"
                if glim.HAVE_PYPROJ else "pyproj not installed"],
               ["area_km2 attribute", f"{b['area_km2']:.6f}"],
               ["attribute basis",
                "ellipsoidal" if glim.HAVE_PYPROJ
                and abs(b["area_km2"] - glim.km2(glim.geodesic_area(b.geom))) < 1e-3
                else "grid"],
               ["repeated vertices", repeated],
               [f"segments < {SHORT_SEGMENT_M} m", short_seg],
               [f"spikes < {SPIKE_DEG} deg", spikes]])
    geodesic = glim.geodesic_area(b.geom) if glim.HAVE_PYPROJ else None
    matches = (abs(b["area_km2"] - glim.km2(b.geom.area)) < 1e-3
               or (geodesic is not None
                   and abs(b["area_km2"] - glim.km2(geodesic)) < 1e-3))
    if not matches:
        record("J_boundary_area", glim.BOUNDARY_LAYER, "boundary",
               "area_km2 agrees with neither the grid nor the ellipsoidal area",
               "fix")

    # ---------------------------------------------------------- summary
    rep.h2("Summary by check")
    by_check = Counter(f[0] for f in findings)
    severity = {f[0]: f[3] for f in findings}
    rep.table(["check", "records", "severity"],
              [[k, n, severity.get(k, "")] for k, n in sorted(by_check.items())])
    rep.kv("total records flagged", len(findings))
    rep.blank()
    rep.line("Severity 'fix' marks a defect that is wrong in any reading: an invalid")
    rep.line("geometry, a duplicate, a fragment, an overlap between the two polygons of")
    rep.line("one landslide. Severity 'check' marks something that needs a mapping")
    rep.line("decision and may well be correct as it stands.")

    rep.h2("Worklist, ordered")
    rep.line("Grouped by what has to happen, not by check, so it can be worked through")
    rep.line("from top to bottom. Counts are records, not pairs.")
    rep.blank()
    order = ["fix", "check", "cosmetic"]
    for level in order:
        subset = [f for f in findings if f[3] == level]
        if not subset:
            continue
        label = {"fix": "1 - Wrong in any reading: correct before resubmitting",
                 "check": "2 - Needs a mapping decision: inspect, then correct or document",
                 "cosmetic": "3 - Cosmetic: harmless at 1:10,000, fix only if convenient"}[level]
        rep.line(label)
        by_layer = defaultdict(list)
        for f in subset:
            by_layer[f[1]].append(f)
        for layer_name, items in sorted(by_layer.items()):
            keys = sorted({i[2] for i in items if i[2] not in ("-", "layer-wide")})
            checks = ", ".join(sorted({i[0] for i in items}))
            rep.line(f"    {short(layer_name):14s} {len(items):3d} records   {checks}")
            if keys and len(keys) <= 40:
                label = "gm_id" if layer_name == glim.GEOMORPHOLOGY_LAYER else "key"
                rep.line(f"                   {label}: {', '.join(keys)}")
        rep.blank()

    rep.h2("Selection expressions")
    rep.line("Built on the attribute keys, so they keep working after records are deleted.")
    rep.blank()
    by_layer = defaultdict(set)
    for f in findings:
        if f[2] not in ("-", "layer-wide", "boundary") and f[3] != "cosmetic":
            by_layer[f[1]].add(f[2])
    for layer_name, keys in sorted(by_layer.items()):
        rep.line(f"  {layer_name}")
        if layer_name == glim.GEOMORPHOLOGY_LAYER:
            ids = sorted(int(k) for k in keys)
            rep.line(f'    "gm_id" IN ({", ".join(str(i) for i in ids)})')
        elif layer_name in glim.LANDSLIDE_LAYERS:
            pairs = sorted((int(k.split("/")[0]), k.split("/")[1]) for k in keys)
            ids = sorted({p[0] for p in pairs})
            rep.line(f'    "lds_id" IN ({", ".join(str(i) for i in ids)})')
            rep.line("    per polygon: " +
                     " OR ".join(f'("lds_id"={i} AND "lds_class"=\'{c}\')'
                                 for i, c in pairs[:6]) +
                     (" OR ..." if len(pairs) > 6 else ""))
        else:
            rep.line("    no stable key: located by centroid, see the CSV")
        rep.blank()

    glim.write_csv("07_topology_findings",
                   ["check", "layer", "key", "severity", "description",
                    "centroid_e", "centroid_n", "fid_volatile"], findings)
    rep.save()


if __name__ == "__main__":
    main()
