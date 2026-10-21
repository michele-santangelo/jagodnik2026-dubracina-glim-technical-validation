#!/usr/bin/env python3
"""
03 - Area statistics: cumulative area, union footprint, overlap and density.

The sum of the polygon areas is not the land area affected by landsliding, because
polygons of different generations overlap and because each deep-seated
landslide is represented by two polygons. The script reports, separately:

  * the cumulative polygon area,
  * the area per landslide (polygons aggregated by lds_id),
  * the union / dissolved footprint, i.e. the land area actually affected,
  * the overlap between relative-age classes and between generations,
  * landslide density, computed against both the measured boundary polygon and
    the 43.6 km2 hydrological basin area,
  * the size distribution, with the count of records below the 400 m2
    completeness threshold.

Output: outputs/03_area_statistics.txt, outputs/03_area_by_class.csv,
        outputs/03_overlap_matrix.csv
"""

from collections import defaultdict

import glim
from glim import Report


def main():
    rep = Report("03_area_statistics", "03 - Area statistics, overlap and density")
    layers = glim.read_all_landslides(fix=True)
    rockfalls = glim.read_rockfalls(fix=True)
    boundary = glim.read_boundary()

    # ---------------------------------------------------------- area basis
    rep.h2("How areas are measured")
    rep.table(["property", "value"],
              [[k.replace("_", " "), v] for k, v in glim.AREA_BASIS.items()])
    rep.line()
    rep.line("This paragraph belongs in the Data Records: without it a reader who")
    rep.line("recomputes the areas from the projected coordinates gets numbers about")
    rep.line("300 ppm larger and has no way to tell whether the published field is wrong.")
    rep.blank()
    rep.line("Every figure below is given on both bases. Note that all ratios - density,")
    rep.line("percentage of the study area - are identical on the two bases to the")
    rep.line("precision reported, because the scale factor cancels in a ratio of two")
    rep.line("areas measured the same way. Only absolute areas differ, in the fourth")
    rep.line("significant figure.")

    # ------------------------------------------------------------- study area
    rep.h2("Study area")
    rep.kv("DRV_Boundary measured area (km2)", f"{glim.km2(boundary.area):.4f}")
    rep.kv("DRV_Boundary interior rings (holes)", len(boundary.interiors))
    rep.kv("study area, ellipsoidal (km2)", f"{glim.AREA_STUDY_KM2:.4f}")
    rep.kv("study area, planar grid (km2)", f"{glim.AREA_STUDY_GRID_KM2:.4f}")
    boundary_attr = glim.read_layer(glim.ORIGINAL_DATA, glim.BOUNDARY_LAYER,
                                    ["Id", "area_km2"])[0]["area_km2"]
    rep.kv("area_km2 attribute in DRV_Boundary", f"{boundary_attr:.6f}")
    if abs(boundary_attr - glim.AREA_STUDY_KM2) > 1e-3:
        rep.finding("the area_km2 attribute of DRV_Boundary still holds the planar grid "
                    f"area ({boundary_attr:.4f} km2). Every other layer now stores "
                    "ellipsoidal areas, so recompute it with the same setting "
                    f"({glim.AREA_STUDY_KM2:.4f} km2) for consistency.")
    rep.kv("area quoted in the submitted manuscript (km2)", glim.AREA_QUOTED_KM2)
    rep.kv("difference (km2)", f"{glim.AREA_QUOTED_KM2 - glim.AREA_STUDY_KM2:+.4f}")
    rep.line()
    rep.line("The study area is the measured area of DRV_Boundary, 38.4 km2, which is the")
    rep.line("domain that was actually mapped. The 43.6 km2 of the submitted manuscript is")
    rep.line("the hydrological basin area (43.55 km2 in Tosevski 2013); the difference is")
    rep.line("the sub-horizontal karst surfaces, deliberately excluded. Every density and")
    rep.line("coverage figure below refers to the mapped domain.")

    # ------------------------------------------------- cumulative vs union
    rep.h2("Cumulative polygon area versus union footprint")
    all_geoms = [f.geom for feats in layers.values() for f in feats]
    cumulative = sum(g.area for g in all_geoms)
    union = glim.safe_union(all_geoms)
    rf_geoms = [f.geom for f in rockfalls]
    rf_union = glim.safe_union(rf_geoms)
    union_all = glim.safe_union([union, rf_union])

    def pct(area):
        return (f"{100 * area / boundary.area:.1f} % of the study area / "
                f"{100 * area / (glim.AREA_QUOTED_KM2 * 1e6):.1f} % of 43.6 km2")

    def geo(geom_or_list):
        """Ellipsoidal area in km2 of a geometry or of a list of geometries."""
        if not glim.HAVE_PYPROJ:
            return None
        if isinstance(geom_or_list, list):
            return sum(glim.geodesic_area(g) for g in geom_or_list) / 1e6
        return glim.geodesic_area(geom_or_list) / 1e6

    rep.table(
        ["quantity", "km2", "share"],
        [["cumulative polygon area, landslide layers", f"{glim.km2(cumulative):.4f}",
          pct(cumulative)],
         ["union footprint, landslide layers", f"{glim.km2(union.area):.4f}",
          pct(union.area)],
         ["cumulative polygon area, rockfall layer", f"{glim.km2(sum(g.area for g in rf_geoms)):.4f}",
          pct(sum(g.area for g in rf_geoms))],
         ["union footprint, rockfall layer", f"{glim.km2(rf_union.area):.4f}",
          pct(rf_union.area)],
         ["union footprint, all layers together", f"{glim.km2(union_all.area):.4f}",
          pct(union_all.area)]])
    overlap = cumulative - union.area
    rep.blank()
    rep.kv("area double counted by the cumulative sum (km2)", f"{glim.km2(overlap):.4f}")
    rep.kv("  as a share of the cumulative area", f"{100 * overlap / cumulative:.1f} %")
    rep.blank()
    rep.finding(
        f"the manuscript reports 20.6 km2 = 47 % of the study area. That is the cumulative "
        f"sum divided by 43.6 km2 ({100 * cumulative / (glim.AREA_QUOTED_KM2 * 1e6):.1f} %). "
        f"The land area actually affected by landsliding is the union, "
        f"{glim.km2(union.area):.2f} km2, i.e. {100 * union.area / (glim.AREA_QUOTED_KM2 * 1e6):.1f} % "
        f"of 43.6 km2 or {100 * union.area / boundary.area:.1f} % of the mapped boundary.")
    rep.blank()
    rep.kv("union entirely inside the boundary",
           "yes" if union.difference(boundary).area < 1.0 else "NO")

    if glim.HAVE_PYPROJ:
        rep.blank()
        rep.line("The same quantities on the ellipsoidal basis, which is the basis the")
        rep.line("dataset's Shape_Area field uses and the one the paper quotes:")
        rep.table(
            ["quantity", "grid km2", "ellipsoidal km2", "% of study area"],
            [["study area", f"{glim.km2(boundary.area):.4f}", f"{geo(boundary):.4f}",
              "100.0"],
             ["cumulative polygon area", f"{glim.km2(cumulative):.4f}",
              f"{geo(all_geoms):.4f}", f"{100 * geo(all_geoms) / geo(boundary):.1f}"],
             ["union footprint, landslides", f"{glim.km2(union.area):.4f}",
              f"{geo(union):.4f}", f"{100 * geo(union) / geo(boundary):.1f}"],
             ["union footprint, rockfalls", f"{glim.km2(rf_union.area):.4f}",
              f"{geo(rf_union):.4f}", f"{100 * geo(rf_union) / geo(boundary):.1f}"],
             ["union footprint, all layers", f"{glim.km2(union_all.area):.4f}",
              f"{geo(union_all):.4f}",
              f"{100 * geo(union_all) / geo(boundary):.1f}"]])

    # --------------------------------------------------- per relative age
    rep.h2("Area per relative-age class and generation")
    rows = []
    per_layer_union = {}
    for name, feats in layers.items():
        geoms = [f.geom for f in feats]
        u = glim.safe_union(geoms)
        per_layer_union[name] = u
        rows.append([name.replace("DRV_GLIM_", ""), glim.RELATIVE_AGE[name],
                     glim.GENERATION[name] or "-", len(feats),
                     f"{glim.km2(sum(g.area for g in geoms)):.4f}",
                     f"{glim.km2(u.area):.4f}",
                     f"{100 * u.area / boundary.area:.2f}"])
    rep.table(["layer", "age", "gen", "polygons", "cumulative km2", "union km2",
               "% of boundary"], rows)

    by_age = defaultdict(list)
    for name, feats in layers.items():
        by_age[glim.RELATIVE_AGE[name]].extend(f.geom for f in feats)
    rep.blank()
    rows = []
    for age in ("O", "VO", "VOR"):
        u = glim.safe_union(by_age[age])
        rows.append([age, len(by_age[age]),
                     f"{glim.km2(sum(g.area for g in by_age[age])):.4f}",
                     f"{glim.km2(u.area):.4f}",
                     f"{100 * u.area / boundary.area:.2f}"])
    rep.table(["relative age", "polygons", "cumulative km2", "union km2",
               "% of boundary"], rows)

    # ------------------------------------------------------- overlap matrix
    rep.h2("Overlap between generation layers (km2)")
    names = list(layers)
    short = [n.replace("DRV_GLIM_", "") for n in names]
    matrix = []
    for i, a in enumerate(names):
        row = []
        for j, b in enumerate(names):
            if j < i:
                row.append("")
            else:
                inter = per_layer_union[a].intersection(per_layer_union[b]).area
                row.append(f"{glim.km2(inter):.4f}")
        matrix.append([short[i]] + row)
    rep.table(["layer"] + short, matrix)
    rep.line()
    rep.line("The diagonal is the union area of the layer itself; off-diagonal cells are")
    rep.line("the areas where two generations overlap. This is the quantitative basis for")
    rep.line("stating that generation polygons are stacked rather than adjacent.")
    glim.write_csv("03_overlap_matrix", ["layer"] + short, matrix)

    # ------------------------------------------------- recurrence / stacking
    rep.h2("Stacking depth (how many generation layers cover the same ground)")
    rep.line("Computed by successive intersection of the per-layer unions: the area")
    rep.line("covered by at least n layers.")
    stack_rows = []
    remaining = per_layer_union.copy()
    covered = []
    for n in range(1, len(names) + 1):
        # area covered by at least n of the layer unions, via pairwise intersections
        from itertools import combinations
        at_least = None
        for combo in combinations(names, n):
            inter = per_layer_union[combo[0]]
            for other in combo[1:]:
                inter = inter.intersection(per_layer_union[other])
                if inter.is_empty:
                    break
            if inter.is_empty:
                continue
            at_least = inter if at_least is None else glim.safe_union([at_least, inter])
        if at_least is None or at_least.is_empty:
            break
        stack_rows.append([n, f"{glim.km2(at_least.area):.4f}",
                           f"{100 * at_least.area / boundary.area:.2f}"])
        covered.append(at_least)
    rep.table(["covered by >= n layers", "km2", "% of boundary"], stack_rows)

    # ------------------------------------------------------------- per landslide
    rep.h2("Size distribution")
    per_landslide = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            per_landslide.append(sum(f.geom.area for f in group))
    per_landslide.sort()
    polygons = sorted(g.area for g in all_geoms)

    def quantile(values, q):
        return values[min(len(values) - 1, int(q * len(values)))]

    rep.table(
        ["statistic", "per landslide (aggregated by lds_id)", "per polygon"],
        [["n", len(per_landslide), len(polygons)],
         ["min m2", f"{per_landslide[0]:.1f}", f"{polygons[0]:.1f}"],
         ["5th percentile m2", f"{quantile(per_landslide, 0.05):.1f}",
          f"{quantile(polygons, 0.05):.1f}"],
         ["median m2", f"{quantile(per_landslide, 0.5):.1f}",
          f"{quantile(polygons, 0.5):.1f}"],
         ["mean m2", f"{sum(per_landslide) / len(per_landslide):.1f}",
          f"{sum(polygons) / len(polygons):.1f}"],
         ["95th percentile m2", f"{quantile(per_landslide, 0.95):.1f}",
          f"{quantile(polygons, 0.95):.1f}"],
         ["max m2", f"{per_landslide[-1]:.1f}", f"{polygons[-1]:.1f}"]])

    thr = glim.COMPLETENESS_THRESHOLD_M2
    below_l = sum(1 for a in per_landslide if a < thr)
    below_p = sum(1 for a in polygons if a < thr)
    rep.blank()
    rep.kv(f"landslides below the {thr:.0f} m2 threshold",
           f"{below_l} ({100 * below_l / len(per_landslide):.1f} %)")
    rep.kv(f"polygons below the {thr:.0f} m2 threshold",
           f"{below_p} ({100 * below_p / len(polygons):.1f} %)")
    rep.finding(
        f"the smallest landslide in the database measures {per_landslide[0]:.1f} m2 and "
        f"{below_l} landslides fall below the {thr:.0f} m2 completeness threshold. The threshold is "
        "the smallest *systematically* mapped unit, not the smallest record.")

    # ----------------------------------------------------------------- density
    rep.h2("Landslide density")
    rep.table(
        ["denominator", "km2", "landslides per km2"],
        [["measured boundary polygon", f"{glim.km2(boundary.area):.4f}",
          f"{len(per_landslide) / glim.km2(boundary.area):.1f}"],
         ["43.6 km2 quoted in the text", f"{glim.AREA_QUOTED_KM2:.4f}",
          f"{len(per_landslide) / glim.AREA_QUOTED_KM2:.1f}"]])
    rep.line()
    rep.line("The manuscript reports 25 landslides per km2, which corresponds to the")
    rep.line("43.6 km2 denominator and to the earlier count of 1,094 landslides.")

    # ------------------------------------------------------------------ export
    rows = []
    for name, feats in layers.items():
        geoms = [f.geom for f in feats]
        u = per_layer_union[name]
        rows.append([name, glim.RELATIVE_AGE[name], glim.GENERATION[name] or "",
                     len(feats), len(glim.group_by_id(feats)),
                     round(sum(g.area for g in geoms), 1), round(u.area, 1)])
    rows.append([glim.ROCKFALL_LAYER, "", "", len(rockfalls), "",
                 round(sum(g.area for g in rf_geoms), 1), round(rf_union.area, 1)])
    rows.append(["ALL LANDSLIDE LAYERS", "", "", len(all_geoms), len(per_landslide),
                 round(cumulative, 1), round(union.area, 1)])
    rows.append(["ALL LAYERS INCL. ROCKFALLS", "", "",
                 len(all_geoms) + len(rf_geoms), "",
                 round(cumulative + sum(g.area for g in rf_geoms), 1),
                 round(union_all.area, 1)])
    glim.write_csv("03_area_by_class",
                   ["layer", "relative_age", "generation", "polygons", "landslides",
                    "cumulative_area_m2", "union_area_m2"], rows)
    rep.save()


if __name__ == "__main__":
    main()
