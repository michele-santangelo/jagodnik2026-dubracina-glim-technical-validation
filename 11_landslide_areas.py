#!/usr/bin/env python3
"""
11 - One row per landslide, with its identifier and its total area.

A table in which every landslide appears once with its total area, because a landslide can be represented by more than one polygon
and summing Shape_Area over the polygon table therefore does not give the area
of the landslides: it double-counts nothing, but it answers a different
question, and a user who groups by lds_id has to know to do so.

The dataset now carries `lds_area`, the sum of Shape_Area over the polygons of
one landslide, so the answer is in the data. This script writes it out as one
CSV per landslide layer, and checks the stored value against the sum recomputed
here and against the geometry, so the published field is verified rather than
trusted.

Rockfalls, the geomorphological layer and the boundary are excluded: the first
has one polygon per record, and the other two are not landslides.

Output: outputs/11_landslide_areas_<layer>.csv, one per layer,
        outputs/11_landslide_areas.txt
"""

import glim
from glim import Report

TOL_M2 = 0.05


def main():
    rep = Report("11_landslide_areas",
                 "11 - One row per landslide, with identifier and total area")
    layers = glim.read_all_landslides(fix=True)

    rep.line("One CSV per landslide layer. Columns:")
    rep.blank()
    rep.line("  lds_id          identifier, unique within the layer")
    rep.line("  lds_class       type; the source-area polygon carries the x prefix")
    rep.line("  n_polygons      how many polygons represent this landslide")
    rep.line("  source_area_m2  area of the scarp/source polygons, blank where none")
    rep.line("  deposit_area_m2 area of the deposit polygons")
    rep.line("  total_area_m2   total area, the published lds_area field")
    rep.line("  sum_shape_area  the same total recomputed from Shape_Area here")
    rep.line("  geometry_m2     planar area of the polygons, for reference")
    rep.line("  depth, uncert   thematic attributes, constant within a landslide")
    rep.line("  cluster_id      overlap cluster, NA where not applicable")
    rep.blank()
    rep.line("lds_area is ellipsoidal, as Shape_Area is; geometry_m2 is planar on")
    rep.line("the layer CRS and runs about 300 ppm larger. See script 03.")
    rep.blank()

    summary = []
    for name, feats in layers.items():
        groups = glim.group_by_id(feats)
        rows = []
        mismatch = 0
        for lds_id, group in sorted(groups.items()):
            deposit = [f for f in group
                       if not glim.is_source_polygon(f["lds_class"])]
            shown = (deposit or group)[0]
            recomputed = sum(f["Shape_Area"] for f in group)
            stored = group[0]["lds_area"]
            if stored is None:
                stored = recomputed
            elif abs(stored - recomputed) > TOL_M2:
                mismatch += 1
            # the source and deposit areas separately as well as the total,
            # so that a user need not perform a dissolve to get them.
            src = sum(f["Shape_Area"] for f in group
                      if glim.is_source_polygon(f["lds_class"]))
            dep = sum(f["Shape_Area"] for f in group
                      if not glim.is_source_polygon(f["lds_class"]))
            rows.append([lds_id, shown["lds_class"], len(group),
                         round(src, 2) if src else "",
                         round(dep, 2) if dep else "",
                         round(stored, 2), round(recomputed, 2),
                         round(sum(f.geom.area for f in group), 2),
                         shown["depth"], shown["uncert"],
                         group[0]["cluster_id"] or "NA"])
        short = name.replace("DRV_GLIM_", "")
        glim.write_csv(f"11_landslide_areas_{short}",
                       ["lds_id", "lds_class", "n_polygons",
                        "source_area_m2", "deposit_area_m2", "total_area_m2",
                        "sum_shape_area", "geometry_m2", "depth", "uncert",
                        "cluster_id"], rows)
        total = sum(r[5] for r in rows)
        summary.append([short, len(feats), len(rows), f"{glim.km2(total):.4f}",
                        mismatch])

    rep.table(["layer", "polygons", "landslides", "total area km2",
               "lds_area disagreeing with the sum"], summary)
    rep.blank()
    rep.kv("landslides, all layers", sum(r[2] for r in summary))
    rep.kv("polygons, all layers", sum(r[1] for r in summary))
    bad = sum(r[4] for r in summary)
    rep.kv("records where the published lds_area does not match the sum", bad)
    if bad:
        rep.finding(f"{bad} landslides carry an lds_area that is not the sum of "
                    "the Shape_Area of their polygons.")
    rep.blank()
    rep.line("The GIS instruction for the Usage Notes, in QGIS expression form,")
    rep.line("for users who prefer to recompute it rather than read the field:")
    rep.blank()
    rep.line('  sum("Shape_Area", group_by:="lds_id")')
    rep.blank()
    rep.line("and, to obtain one record per landslide, dissolve on lds_id or use")
    rep.line("Processing > Vector geometry > Dissolve with lds_id as the field.")
    rep.save()


if __name__ == "__main__":
    main()
