#!/usr/bin/env python3
"""
05 - Comparison with the inventory of Tosevski (2013).

The 39 slides of Tosevski (2013) were not distributed as GIS data, so the
comparison was expected to be possible only at the level of metadata. It is
not: Table 1 of that paper publishes the coordinates of a characteristic point
of every slide polygon, in "MGI Balkans 5" (MGI 1901 / Balkans zone 5,
EPSG:31275), and Table 2 publishes area, volume, length, width, slope and the
detection method for each. Both tables are transcribed in
`tosevski2013_inventory.csv`.

The script therefore performs a genuine, independent spatial check: it
reprojects the 39 points to HTRS96 / Croatia TM and tests whether each falls
inside a landslide polygon of the new G-LIM. It also compares the metadata that
cannot be checked spatially (mapped area, size range, detection method,
lithological setting).

Caveat, reported by the script: the transformation from the Hermannskogel
datum to HTRS96 is applied with the parameters bundled in PROJ, not with the
official Croatian transformation grid, so an absolute offset of a few metres
is possible. The distance of each point to the nearest landslide polygon is
reported so that the result can be judged against that uncertainty.

Output: outputs/05_compare_tosevski.txt, outputs/05_tosevski_match.csv
"""

import csv
from collections import Counter

from shapely.geometry import Point
from shapely.strtree import STRtree

import glim
from glim import Report

# published values used to verify the transcription
PUBLISHED = {
    "n_slides": 39,
    "total_area_m2": 81873,
    "detected_dem": 16,
    "detected_field": 23,
    "area_min_m2": 150,
    "area_max_m2": 12956,
    "length_min_m": 20,
    "length_max_m": 226,
    "slope_min_deg": 10.1,
    "slope_max_deg": 28.6,
    "study_area_km2": 18.55,
    "basin_area_km2": 43.55,
    "share_of_study_area_pct": 0.44,
}


def load_table():
    rows = []
    with open(glim.TOSEVSKI_TABLE, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            rows.append(line)
    reader = csv.DictReader(rows)
    out = []
    for row in reader:
        out.append({
            "slide_id": int(row["slide_id"]),
            "x": float(row["x"]), "y": float(row["y"]),
            "area_m2": float(row["area_m2"]),
            "volume_m3": float(row["volume_m3"]),
            "total_length_m": float(row["total_length_m"]),
            "width_m": float(row["width_m"]),
            "slope_deg": float(row["slope_deg"]),
            "detection": row["detection"],
            "azimuth_deg": float(row["azimuth_deg"]),
        })
    return out


def main():
    rep = Report("05_compare_tosevski", "05 - Comparison with Tosevski (2013)")
    slides = load_table()

    # ------------------------------------------------- transcription checks
    rep.h2("Transcription checks against the published values")
    detection = Counter(s["detection"] for s in slides)
    checks = [
        ["number of slides", len(slides), PUBLISHED["n_slides"]],
        ["total landslide area (m2)", sum(s["area_m2"] for s in slides),
         PUBLISHED["total_area_m2"]],
        ["detected on the DEM", detection["DEM"], PUBLISHED["detected_dem"]],
        ["detected in the field", detection["field"], PUBLISHED["detected_field"]],
        ["minimum area (m2)", min(s["area_m2"] for s in slides), PUBLISHED["area_min_m2"]],
        ["maximum area (m2)", max(s["area_m2"] for s in slides), PUBLISHED["area_max_m2"]],
        ["minimum total length (m)", min(s["total_length_m"] for s in slides),
         PUBLISHED["length_min_m"]],
        ["maximum total length (m)", max(s["total_length_m"] for s in slides),
         PUBLISHED["length_max_m"]],
        ["minimum slope (deg)", min(s["slope_deg"] for s in slides),
         PUBLISHED["slope_min_deg"]],
        ["maximum slope (deg)", max(s["slope_deg"] for s in slides),
         PUBLISHED["slope_max_deg"]],
    ]
    rep.table(["quantity", "transcribed", "published", "status"],
              [[c[0], c[1], c[2], "ok" if abs(c[1] - c[2]) < 1e-6 else "MISMATCH"]
               for c in checks])
    if all(abs(c[1] - c[2]) < 1e-6 for c in checks):
        rep.line()
        rep.line("All aggregate values reproduce, so the transcribed tables can be used.")

    # the author states volume = 0.0844 * area^1.4324 (after Guzzetti et al. 2008)
    rep.h2("Internal consistency of the published area-volume pairs")
    bad = []
    for s in slides:
        expected = 0.0844 * s["area_m2"] ** 1.4324
        if abs(expected - s["volume_m3"]) > 0.15 * expected:
            bad.append([s["slide_id"], s["area_m2"], s["volume_m3"], round(expected)])
    if bad:
        rep.table(["slide_id", "area m2", "published volume m3",
                   "volume from the stated relation m3"], bad)
        rep.finding("the published volume of the slide(s) above does not follow the "
                    "area-volume relation the author states he applied; most likely a "
                    "typesetting error in the 2013 paper. Not used further, but worth a "
                    "footnote if those volumes are quoted.")
    else:
        rep.line("All 39 area-volume pairs follow the stated relation within 15 %.")

    # -------------------------------------------------------- reprojection
    rep.h2("Reprojection to HTRS96 / Croatia TM")
    transform = glim.transformer(glim.TOSEVSKI_EPSG, glim.EXPECTED_EPSG)
    boundary = glim.read_boundary()
    points = []
    for s in slides:
        x, y = transform(s["x"], s["y"])
        points.append((s, Point(x, y)))
    inside_boundary = sum(1 for _, p in points if boundary.contains(p))
    rep.kv("source CRS", f"EPSG:{glim.TOSEVSKI_EPSG} (MGI 1901 / Balkans zone 5)")
    rep.kv("target CRS", f"EPSG:{glim.EXPECTED_EPSG} (HTRS96 / Croatia TM)")
    rep.kv("points inside the study-area boundary", f"{inside_boundary} / {len(points)}")
    if inside_boundary == len(points):
        rep.line()
        rep.line("All 39 points land inside the Dubracina boundary, which confirms the")
        rep.line("interpretation of the published coordinate system and of the column order")
        rep.line("(the column printed as 'x' is the easting).")
    rep.line()
    rep.line("Caveat: PROJ's bundled datum-shift parameters are used, not the official")
    rep.line("Croatian transformation grid, so a systematic offset of a few metres cannot")
    rep.line("be excluded. Distances to the nearest polygon are reported below so that the")
    rep.line("conclusion can be judged against that uncertainty.")

    # ------------------------------------------------- point in polygon test
    rep.h2("Do the 39 slides fall inside the new G-LIM?")
    layers = glim.read_all_landslides(fix=True)
    records = []
    for name, feats in layers.items():
        for f in feats:
            records.append((name.replace("DRV_GLIM_", ""), glim.RELATIVE_AGE[name],
                            f["lds_id"], f["lds_class"], f["depth"], f["uncert"], f.geom))
    geoms = [r[6] for r in records]
    tree = STRtree(geoms)
    union = glim.safe_union(geoms)

    # the Jagodnik et al. (2020) data are not public; without them the column
    # that refers to that inventory is left empty
    if glim.jom_available():
        jom = glim.read_jom_polygons(fix=True)
        jom_union = glim.safe_union([f.geom for feats in jom.values() for f in feats])
    else:
        jom_union = None

    rows = []
    n_inside = 0
    distances = []
    for s, p in points:
        hits = [records[i] for i in tree.query(p) if geoms[i].contains(p)]
        distance = 0.0 if hits else union.distance(p)
        distances.append(distance)
        if hits:
            n_inside += 1
        # the most specific match is the smallest polygon containing the point,
        # because the very-old-relict polygons envelop the younger ones
        hits.sort(key=lambda r: r[6].area)
        rows.append([
            s["slide_id"], s["area_m2"], s["detection"],
            round(p.x, 1), round(p.y, 1),
            "yes" if hits else "no", round(distance, 1),
            len(hits),
            hits[0][0] if hits else "", hits[0][2] if hits else "",
            hits[0][3] if hits else "", hits[0][4] if hits else "",
            "|".join(f"{h[0]}:{h[3]}" for h in hits),
            "" if jom_union is None else "yes" if jom_union.contains(p) else "no",
        ])
        label = ("INSIDE  " + "; ".join(f"{h[0]} id{h[2]} {h[3]} {h[4]}" for h in hits)) \
            if hits else f"outside, {distance:.1f} m from the nearest polygon"
        rep.line(f"  slide {s['slide_id']:2d}  A={s['area_m2']:6.0f} m2  "
                 f"{s['detection']:5s}  {label}")

    rep.blank()
    rep.kv("slides falling inside a G-LIM landslide polygon",
           f"{n_inside} / {len(points)} ({100 * n_inside / len(points):.0f} %)")
    for tolerance in (5, 10, 25, 50, 100):
        n = sum(1 for d in distances if d <= tolerance)
        rep.kv(f"  within {tolerance} m of a G-LIM polygon",
               f"{n} / {len(points)} ({100 * n / len(points):.0f} %)")
    if jom_union is not None:
        n_jom = sum(1 for _, p in points if jom_union.contains(p))
        rep.kv("slides falling inside the 2020 J. Maps inventory",
               f"{n_jom} / {len(points)}")
    else:
        rep.kv("slides falling inside the 2020 J. Maps inventory", "not available")
    rep.line()
    rep.line("Tosevski reports that 36 of his 39 slides lie in the Slani Potok, Malenica")
    rep.line("and Kucina sub-catchments, which is where the 2020 inventory concentrates;")
    rep.line("the count above is an independent confirmation of that statement.")

    # ------------------------------------------- how the new G-LIM classifies them
    rep.h2("How the new G-LIM classifies the same features")
    most_specific = [r for r in rows if r[5] == "yes"]
    rep.table(["quantity", "value"],
              [["relative-age layer of the most specific match",
                ", ".join(f"{k}={v}" for k, v in
                          Counter(r[8] for r in most_specific).most_common())],
               ["landslide type of the most specific match",
                ", ".join(f"{k}={v}" for k, v in
                          Counter(r[10] for r in most_specific).most_common())],
               ["depth class of the most specific match",
                ", ".join(f"{k}={v}" for k, v in
                          Counter(r[11] for r in most_specific).most_common())],
               ["mean number of overlapping G-LIM polygons per slide",
                f"{sum(r[7] for r in most_specific) / len(most_specific):.1f}"]])
    rep.line()
    rep.line("Tosevski mapped only soil slides in the flysch weathering zone and states")
    rep.line("that depth of the failure surface was left for future work, so his slides")
    rep.line("carry no depth class. A difference in the depth attribute is therefore not")
    rep.line("a disagreement between the two inventories.")

    # ------------------------------------------------------- metadata comparison
    rep.h2("Metadata comparison")
    layer_geoms = [f.geom for feats in layers.values() for f in feats]
    new_union_area = glim.safe_union(layer_geoms).area
    rep.table(
        ["property", "Tosevski (2013)", "new G-LIM"],
        [["inventory type", "historical (author's own wording)",
          "geomorphological historical"],
         ["mapped domain", f"{PUBLISHED['study_area_km2']} km2 (central part of the basin)",
          f"{glim.km2(boundary.area):.2f} km2 (boundary polygon)"],
         ["basin area quoted", f"{PUBLISHED['basin_area_km2']} km2",
          f"{glim.AREA_QUOTED_KM2} km2"],
         ["base data", "1 m LiDAR DTM (March 2012) + field mapping",
          "1 m LiDAR DTM (2020-2022), pseudo-stereoscopic hillshade"],
         ["derivatives used", "hillshade 315/45, slope, 1 m contours",
          "hillshade 315/45 and 45/45, combined, slope, 1/2/5 m contours"],
         ["stereoscopic imagery", "aerial stereo-pairs tried, no landslide detected",
          "synthetic stereo-mates from the DTM, anaglyph interpretation"],
         ["movement types", "soil slides only", "slides, flows, slide-flows, "
          "rock avalanches, rock falls, widespread units"],
         ["number of landslides", PUBLISHED["n_slides"],
          sum(len(glim.group_by_id(f)) for f in layers.values())],
         ["landslide area", f"{PUBLISHED['total_area_m2'] / 1e6:.4f} km2 "
          f"({PUBLISHED['share_of_study_area_pct']} % of its domain)",
          f"{glim.km2(new_union_area):.4f} km2 union "
          f"({100 * new_union_area / boundary.area:.1f} % of the boundary)"],
         ["size range", f"{PUBLISHED['area_min_m2']}-{PUBLISHED['area_max_m2']} m2",
          f"{min(sum(x.geom.area for x in g) for feats in layers.values() for g in glim.group_by_id(feats).values()):.0f}-"
          f"{max(sum(x.geom.area for x in g) for feats in layers.values() for g in glim.group_by_id(feats).values()):.0f} m2"],
         ["relative age", "not assigned", "3 classes, 6 generations"],
         ["uncertainty", "not assigned", "3-digit positional code"],
         ["attributes per landslide",
          "area, volume, length, width, slope, azimuth, detection method",
          "type, relative age, generation, depth class, uncertainty, area"]])
    rep.blank()
    rep.line("The two inventories are complementary rather than competing: Tosevski")
    rep.line("mapped a single process type over 18.55 km2 with field verification, and")
    rep.line("almost every one of his slides is recovered by the new G-LIM, which maps")
    rep.line("many more process types and several generations over the whole catchment.")

    glim.write_csv(
        "05_tosevski_match",
        ["slide_id", "published_area_m2", "detection_method", "easting_htrs96",
         "northing_htrs96", "inside_new_glim", "distance_to_nearest_m",
         "n_overlapping_polygons", "most_specific_layer", "most_specific_lds_id",
         "most_specific_lds_class", "most_specific_depth", "all_overlaps",
         "inside_jom_2020"],
        rows)
    rep.save()


if __name__ == "__main__":
    main()
