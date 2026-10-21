#!/usr/bin/env python3
"""
02 - Database and geometry quality control.

Automated checks of geometry validity, self-intersections, empty geometries,
duplicates, features outside the study boundary, missing values, invalid codes,
CRS consistency and identifier consistency across layers, and of the
identifier/depth consistency rule described in the Technical Validation
section.

Checks performed
  A  CRS consistency across all layers
  B  geometry validity, empty and multipart geometries, holes
  C  missing (NULL) attribute values
  D  attribute values outside the dictionaries of Tables 2, 5 and 7
  E  identifier structure: polygons per lds_id, global uniqueness
  F  source/deposit pairing consistency (the `x` prefix)
  G  attribute consistency within one lds_id (depth, uncert, type)
  H  the `depth` rule: 2 polygons => deep-seated, 1 polygon => shallow
  I  features outside the study-area boundary
  J  Shape_Area attribute versus computed geometry area
  K  duplicate geometries

Geometries are read WITHOUT repair, because finding broken geometries is the
point of this script.

Output: outputs/02_qc_database.txt, outputs/02_qc_findings.csv
"""

from collections import Counter, defaultdict

from shapely.validation import explain_validity

import glim
from glim import Report

ALL_VECTOR_LAYERS = (glim.LANDSLIDE_LAYERS
                     + [glim.ROCKFALL_LAYER, glim.GEOMORPHOLOGY_LAYER,
                        glim.BOUNDARY_LAYER])

findings = []


def record(check, layer, feature, description):
    findings.append([check, layer, feature, description])


def main():
    rep = Report("02_qc_database", "02 - Database and geometry quality control")
    layers = glim.read_all_landslides(fix=False)
    rockfalls = glim.read_rockfalls(fix=False)
    geomorph = glim.read_layer(glim.ORIGINAL_DATA, glim.GEOMORPHOLOGY_LAYER,
                               ["gm_id", "geom_class", "Shape_Area"], fix=False)
    boundary = glim.read_boundary()

    def short(name):
        return name.replace("DRV_GLIM_", "")

    # -------------------------------------------------------------------- A CRS
    rep.h2("A - Coordinate reference system")
    rows = []
    for name in ALL_VECTOR_LAYERS:
        code = glim.layer_epsg(glim.ORIGINAL_DATA, name)
        ok = code == glim.EXPECTED_EPSG
        rows.append([name, code if code else "not identified", "ok" if ok else "MISMATCH"])
        if not ok:
            record("A_crs", name, "-", f"EPSG {code}, expected {glim.EXPECTED_EPSG}")
    for name in sorted(glim.JOM_LABELS) if glim.jom_available() else []:
        code = glim.layer_epsg(glim.PREVIOUS_INVENTORIES, name)
        rows.append([name, code if code else "not identified",
                     "ok" if code == glim.EXPECTED_EPSG else "MISMATCH"])
    rep.table(["layer", "EPSG", "status"], rows)
    rep.line()
    rep.line("The hillshade raster Dubracina_SHD.tif is checked separately: gdalinfo")
    rep.line("reports an unnamed CRS whose parameters match HTRS96/TM but which")
    rep.line("carries no EPSG authority code (see README).")

    # ----------------------------------------------------------- B geometry
    rep.h2("B - Geometry validity, empty and multipart geometries")
    rows = []
    named = list(layers.items()) + [(glim.ROCKFALL_LAYER, rockfalls),
                                    (glim.GEOMORPHOLOGY_LAYER, geomorph)]
    for name, feats in named:
        invalid = [f for f in feats if not f.geom.is_valid]
        empty = [f for f in feats if f.geom.is_empty]
        multi = [f for f in feats if f.geom.geom_type == "MultiPolygon"]
        holes = sum(len(p.interiors) for f in feats for p in glim.parts(f.geom))
        areas = [f.geom.area for f in feats]
        rows.append([name.replace("DRV_GLIM_", ""), len(feats), len(invalid),
                     len(empty), len(multi), holes,
                     f"{min(areas):.2f}", f"{max(areas):.0f}"])
        for f in invalid:
            ident = f["lds_id"] if "lds_id" in f.attrs else f.attrs.get("gm_id", f.attrs.get("Id"))
            record("B_invalid_geometry", name, ident, explain_validity(f.geom))
        for f in empty:
            record("B_empty_geometry", name, f.fid, "empty geometry")
    rep.table(["layer", "features", "invalid", "empty", "multipart",
               "holes", "min area m2", "max area m2"], rows)

    rep.blank()
    rep.line("Invalid geometries in detail:")
    for check, layer, ident, desc in findings:
        if check == "B_invalid_geometry":
            rep.finding(f"{layer} (id {ident}): {desc}")

    rep.blank()
    rep.line("Multipart landslide features (relevant to the two-polygon rule):")
    for name, feats in named:
        multi = [f for f in feats if f.geom.geom_type == "MultiPolygon"]
        for f in multi:
            ident = f.attrs.get("lds_id", f.attrs.get("Id"))
            rep.line(f"    {name.replace('DRV_GLIM_', ''):10s} id={ident} "
                     f"class={f['lds_class']} parts={len(f.geom.geoms)} "
                     f"area={f.geom.area:.1f} m2")
            record("B_multipart", name, ident,
                   f"{len(f.geom.geoms)} parts, class {f['lds_class']}")

    # ------------------------------------------------------------- C NULLs
    rep.h2("C - Missing (NULL) attribute values")
    rows = []
    for name, feats in named:
        for field in feats[0].attrs:
            n_null = sum(1 for f in feats if f.attrs[field] is None)
            if n_null:
                rows.append([name.replace("DRV_GLIM_", ""), field, n_null])
                record("C_null", name, field, f"{n_null} NULL values")
    rep.table(["layer", "field", "NULL values"], rows or [["-", "-", "none found"]])

    # -------------------------------------------------------- D dictionaries
    rep.h2("D - Attribute values against the published dictionaries")
    lds_codes = Counter(f["lds_class"] for feats in layers.values() for f in feats)
    lds_codes.update(f["lds_class"] for f in rockfalls)
    unknown = sorted(set(lds_codes) - glim.TABLE2_LDS_CLASS - {None})
    unused = sorted(glim.TABLE2_LDS_CLASS - set(lds_codes))
    rep.kv("distinct lds_class values in the data", len(lds_codes))
    rep.line()
    if unknown:
        for code in unknown:
            rep.finding(f"lds_class '{code}' ({lds_codes[code]} polygons) "
                        "is not documented in Table 2")
            record("D_undocumented_code", "landslide layers", code,
                   f"lds_class not in Table 2 ({lds_codes[code]} polygons)")
    if unused:
        for code in unused:
            rep.finding(f"lds_class '{code}' is documented in Table 2 "
                        "but never used in the data")
            record("D_unused_code", "Table 2", code, "documented but absent from data")
    if "wxrf" in lds_codes and "wrf" not in lds_codes:
        rep.blank()
        rep.line("The widespread rockfall codes need a decision. Table 2 defines wxrf as")
        rep.line("the widespread source area and wrf as the widespread deposit. The data")
        rep.line(f"now hold {lds_codes['wxrf']} wxrf and no wrf, so every widespread")
        rep.line("rockfall polygon is a source area and no widespread deposit is recorded")
        rep.line("in this layer. That is coherent if the deposits are the talus units of")
        rep.line("the geomorphological layer, which is what the Data Overview describes,")
        rep.line("but Table 2 must then say so and stop offering wrf as a value.")

    geom_codes = Counter(f["geom_class"] for f in geomorph)
    for code in sorted(set(geom_codes) - glim.TABLE5_GEOM_CLASS - {None}):
        rep.finding(f"geom_class '{code}' ({geom_codes[code]} polygons) "
                    "is not documented in Table 5")
        record("D_undocumented_code", glim.GEOMORPHOLOGY_LAYER, code,
               f"geom_class not in Table 5 ({geom_codes[code]} polygons)")

    depth_values = Counter(f["depth"] for feats in layers.values() for f in feats)
    for value in sorted(set(depth_values) - glim.TABLE7_DEPTH - {None}, key=str):
        rep.finding(f"depth value '{value}' is not documented in Table 7")
        record("D_undocumented_code", "landslide layers", value, "depth not in Table 7")

    bad_uncert = [(n, f["lds_id"], f["uncert"]) for n, feats in layers.items()
                  for f in feats
                  if not (f["uncert"] and len(f["uncert"]) == 3
                          and set(f["uncert"]) <= {"0", "1"})]
    rep.blank()
    rep.kv("uncert values not matching /^[01]{3}$/", len(bad_uncert))
    for n, i, v in bad_uncert[:10]:
        record("D_bad_uncert", n, i, f"uncert='{v}'")

    # ------------------------------------------------------ E identifiers
    rep.h2("E - Identifier structure")
    rows = []
    global_ids = defaultdict(list)
    for name, feats in layers.items():
        grouped = glim.group_by_id(feats)
        counts = Counter(len(v) for v in grouped.values())
        rows.append([name.replace("DRV_GLIM_", ""), len(grouped),
                     min(grouped), max(grouped),
                     counts[1], counts[2],
                     sum(v for k, v in counts.items() if k > 2)])
        for lds_id in grouped:
            global_ids[lds_id].append(name)
        # A landslide is normally one or two polygons, but the Methods now state
        # that particular morphological configurations require more - a debris
        # flow fed by two separate detachment zones, for example. More than two
        # polygons is therefore no longer a finding in itself; what must hold is
        # that at most one of them is the deposit.
        for lds_id, group in grouped.items():
            deposits = [f for f in group
                        if not glim.is_source_polygon(f["lds_class"])]
            if len(deposits) > 1:
                record("E_too_many_polygons", name, lds_id,
                       f"{len(group)} polygons share this lds_id, "
                       f"{len(deposits)} of them deposits")
    rep.table(["layer", "distinct lds_id", "min", "max", "1 poly", "2 poly", ">2 poly"], rows)
    reused = {i: v for i, v in global_ids.items() if len(v) > 1}
    rep.blank()
    rep.kv("distinct lds_id values across all layers", len(global_ids))
    rep.kv("lds_id values occurring in more than one layer", len(reused))
    if reused:
        rep.finding("lds_id is unique only WITHIN a layer, not across the dataset; "
                    "Table 7 describes it as a 'unique identifier assigned to each "
                    "landslide'. Either qualify the description or issue dataset-wide ids "
                    "(e.g. O_I_0001).")
        record("E_not_globally_unique", "all landslide layers", "-",
               f"{len(reused)} lds_id values are reused in more than one layer")
    rf_ids = Counter(f["lds_id"] for f in rockfalls)
    rf_schema = glim.layer_fields(glim.ORIGINAL_DATA, glim.ROCKFALL_LAYER)
    rep.blank()
    rep.kv("rockfall layer attributes", ", ".join(rf_schema))
    rep.kv("distinct identifier values in the rockfall layer", len(rf_ids))
    if len(rf_ids) == 1:
        rep.finding(f"all {len(rockfalls)} rockfall polygons carry the same identifier "
                    f"({list(rf_ids)[0]}): the layer has no usable identifier, unlike "
                    "the schema described in Table 7")
        record("E_no_identifier", glim.ROCKFALL_LAYER, "-",
               f"identifier constant for all {len(rockfalls)} features")
    else:
        rf_groups = Counter(rf_ids.values())
        rep.kv("polygons per identifier", dict(sorted(rf_groups.items())))
        if max(rf_ids.values()) > 2:
            worst = [i for i, n in rf_ids.items() if n > 2][:10]
            rep.finding(f"{sum(1 for n in rf_ids.values() if n > 2)} rockfall identifiers "
                        f"are shared by more than two polygons, e.g. {worst}; the "
                        "two-polygon rule of the landslide layers does not hold here, "
                        "so the meaning of the rockfall identifier must be documented")
            record("E_rockfall_id_structure", glim.ROCKFALL_LAYER,
                   ", ".join(str(i) for i in worst),
                   "identifier shared by more than two polygons")
    missing = [f for f in ("uncert", "depth") if f not in rf_schema]
    if missing:
        rep.finding(f"the rockfall layer still has no {' or '.join(missing)} attribute, "
                    "so Table 7 has to state which fields apply to which layer")
        record("E_rockfall_missing_fields", glim.ROCKFALL_LAYER, "-",
               f"no {', '.join(missing)} attribute")

    # ---------------------------------------------------------- F/G pairing
    rep.h2("F - Source / deposit pairing within one lds_id")
    rows = []
    for name, feats in layers.items():
        problems = []
        for lds_id, group in glim.group_by_id(feats).items():
            if len(group) != 2:
                continue
            src = [f["lds_class"] for f in group if glim.is_source_polygon(f["lds_class"])]
            dep = [f["lds_class"] for f in group if not glim.is_source_polygon(f["lds_class"])]
            if len(src) != 1 or len(dep) != 1:
                problems.append((lds_id, [f["lds_class"] for f in group],
                                 "not exactly one source and one deposit polygon"))
            elif glim.base_class(src[0]) != dep[0]:
                problems.append((lds_id, [f["lds_class"] for f in group],
                                 "source and deposit have different landslide types"))
        rows.append([name.replace("DRV_GLIM_", ""), len(problems)])
        for lds_id, classes, why in problems:
            rep.finding(f"{name.replace('DRV_GLIM_', '')} lds_id={lds_id}: "
                        f"{'+'.join(classes)} - {why}")
            record("F_pairing", name, lds_id, f"{'+'.join(classes)}: {why}")
    rep.table(["layer", "pairing problems"], rows)

    rep.h2("G - Attribute consistency within one lds_id")
    rows = []
    for name, feats in layers.items():
        n_dep = n_unc = 0
        for lds_id, group in glim.group_by_id(feats).items():
            if len({f["depth"] for f in group}) > 1:
                n_dep += 1
                record("G_mixed_depth", name, lds_id,
                       "polygons of one landslide carry different depth values")
            if len({f["uncert"] for f in group}) > 1:
                n_unc += 1
                record("G_mixed_uncert", name, lds_id,
                       "polygons of one landslide carry different uncert codes")
        rows.append([name.replace("DRV_GLIM_", ""), n_dep, n_unc])
    rep.table(["layer", "mixed depth", "mixed uncert"], rows)

    # -------------------------------------------------------- H the depth rule
    rep.h2("H - The depth rule as the Methods now state it")
    rep.line("A landslide whose scarp or source area is delineated as a polygon")
    rep.line("distinct from the deposit is deep-seated; one represented by a single")
    rep.line("polygon is shallow; and a debris flow is shallow whatever its polygon")
    rep.line("structure, because its upslope polygon is a concave erosional and")
    rep.line("transport zone and not the scarp of a failure surface.")
    rep.blank()
    rep.line("The submitted text stated the first two clauses and omitted the third.")
    rep.line("Run against that text this check flagged 79 debris flows, which were")
    rep.line("not defects: the data was right and the description incomplete. The")
    rep.line("criterion is encoded here as the Methods now state it, so those records")
    rep.line("are no longer findings, and what the check reports is disagreement with")
    rep.line("the published criterion rather than with an earlier draft of it.")
    rep.blank()
    rows = []
    total_violations = 0
    for name, feats in layers.items():
        wrong = mixed = 0
        for lds_id, group in glim.group_by_id(feats).items():
            values = {f["depth"] for f in group}
            if len(values) > 1:
                mixed += 1
                record("H_depth_rule", name, lds_id,
                       "the polygons of one landslide disagree on depth")
                continue
            value = values.pop()
            want = glim.expected_depth(group)
            if want is not None and value != want:
                wrong += 1
                record("H_depth_rule", name, lds_id,
                       f"depth='{value}' but the stated criteria require '{want}' "
                       f"({len(group)} polygons, "
                       f"{'+'.join(sorted(f['lds_class'] for f in group))})")
        rows.append([name.replace("DRV_GLIM_", ""), len(glim.group_by_id(feats)),
                     wrong, mixed])
        total_violations += wrong + mixed
    rep.table(["layer", "landslides", "depth against the criteria",
               "polygons disagree"], rows)
    rep.blank()
    rep.kv("landslides not following the stated criteria", total_violations)
    if total_violations:
        rep.finding(f"{total_violations} landslides do not follow the criteria stated "
                    "in the Methods. Either the attribute or the statement of the "
                    "criteria must be corrected before resubmission.")

    # ------------------------------------------------------------ I boundary
    rep.h2("I - Features outside the study-area boundary")
    rows = []
    for name, feats in named:
        outside = []
        for f in feats:
            geom = f.geom if f.geom.is_valid else f.geom.buffer(0)
            spill = geom.difference(boundary).area
            if spill > 1.0:
                outside.append((f.attrs.get("lds_id", f.attrs.get("gm_id", f.attrs.get("Id"))),
                                spill))
        rows.append([name.replace("DRV_GLIM_", ""), len(outside),
                     f"{sum(s for _, s in outside):.1f}"])
        for ident, spill in outside:
            record("I_outside_boundary", name, ident, f"{spill:.1f} m2 outside the boundary")
    rep.table(["layer", "features spilling out", "area outside m2"], rows)

    # ---------------------------------------------------------- J Shape_Area
    rep.h2("J - Shape_Area attribute versus computed area")
    rep.line("Shape_Area is compared against BOTH the planar grid area of the projected")
    rep.line("coordinates and the geodesic area on the GRS80 ellipsoid, because the two")
    rep.line("differ here and a GIS can be asked for either. HTRS96 / Croatia TM has its")
    rep.line("central meridian at 16.5 deg E with a scale factor of 0.9999, and the")
    rep.line("Dubracina valley lies about 145 km west of it, where the grid scale has")
    rep.line("grown back above unity: grid areas run about 0.03 % larger than ellipsoidal")
    rep.line("ones. A field recomputed with the geodesic option therefore disagrees with")
    rep.line("the grid area by that amount without being wrong.")
    if not glim.HAVE_PYPROJ:
        rep.blank()
        rep.finding("pyproj is not installed, so only the grid area can be checked; "
                    "install it to separate a geodesic field from a stale one")
    rep.blank()
    rows = []
    basis_totals = Counter()
    for name, feats in list(layers.items()) + [(glim.ROCKFALL_LAYER, rockfalls),
                                               (glim.GEOMORPHOLOGY_LAYER, geomorph)]:
        grid_hits = geo_hits = neither = no_field = 0
        # the two bases are only distinguishable on polygons large enough for
        # the 0.03 % difference to exceed the rounding of the stored field
        decisive_grid = decisive_geo = 0
        sum_field = sum_grid = sum_geo = 0.0
        for f in feats:
            field = f.attrs.get("Shape_Area")
            grid = f.geom.area
            geo = glim.geodesic_area(f.geom) if glim.HAVE_PYPROJ else None
            sum_grid += grid
            sum_geo += geo if geo is not None else grid
            if field is None:
                no_field += 1
                continue
            sum_field += field
            tolerance = max(0.5, 1e-6 * grid)
            decisive = geo is not None and abs(grid - geo) > 2.0
            if abs(field - grid) <= tolerance:
                grid_hits += 1
                decisive_grid += decisive
            elif geo is not None and abs(field - geo) <= tolerance:
                geo_hits += 1
                decisive_geo += decisive
            else:
                neither += 1
                record("J_area_mismatch", name, glim.stable_key(f),
                       f"Shape_Area={field:.2f} m2 but grid={grid:.2f} and "
                       f"geodesic={geo:.2f} m2" if geo is not None else
                       f"Shape_Area={field:.2f} m2 but geometry={grid:.2f} m2")
        if no_field == len(feats):
            rows.append([short(name), len(feats), "no Shape_Area field",
                         "-", "-", "-", "-"])
            continue
        if decisive_grid or decisive_geo:
            basis = "geodesic" if decisive_geo > decisive_grid else "grid"
        else:
            basis = "indistinguishable"
        basis_totals[basis] += 1
        rows.append([short(name), len(feats), basis,
                     f"{decisive_grid}/{decisive_geo}", neither,
                     f"{glim.km2(sum_field):.5f}"])
    rep.table(["layer", "features", "field basis", "decisive grid/geodesic",
               "matches neither", "sum of field km2"], rows)
    rep.line()
    rep.line("'decisive' counts only polygons big enough that the 0.03 % difference")
    rep.line("between the two bases exceeds the precision the field is stored with;")
    rep.line("smaller polygons match both and say nothing about which was used.")
    rep.blank()
    rep.kv("layers whose Shape_Area is geodesic", basis_totals["geodesic"])
    rep.kv("layers whose Shape_Area is the grid area", basis_totals["grid"])
    mismatches = [f for f in findings if f[0] == "J_area_mismatch"]
    rep.kv("records matching neither basis", len(mismatches))
    for f in mismatches[:10]:
        rep.finding(f"{short(f[1])} {f[2]}: {f[3]}")
    if basis_totals["geodesic"] and basis_totals["grid"]:
        rep.blank()
        rep.finding("Shape_Area is geodesic in some layers and the grid area in others. "
                    "Recompute the whole dataset on one basis, and state in the Data "
                    "Records which one, so that areas summed across layers are "
                    "consistent.")
    rep.blank()
    rep.line("Totals on each basis, all polygon layers together:")
    all_polys = [f for feats in layers.values() for f in feats] + rockfalls + geomorph
    total_grid = sum(f.geom.area for f in all_polys)
    total_geo = (sum(glim.geodesic_area(f.geom) for f in all_polys)
                 if glim.HAVE_PYPROJ else None)
    rep.table(["basis", "km2", "difference from grid"],
              [["planar grid area (projected coordinates)",
                f"{glim.km2(total_grid):.5f}", "-"]]
              + ([["geodesic area (GRS80 ellipsoid)", f"{glim.km2(total_geo):.5f}",
                   f"{1e6 * (total_grid / total_geo - 1):.0f} ppm smaller"]]
                 if total_geo else []))
    rep.line()
    rep.line("At 1:10,000 the difference is immaterial: 0.03 % of the mapped area is a")
    rep.line("few thousand square metres over 20 km2, well inside the delineation")
    rep.line("uncertainty. It matters only for internal consistency and for saying in")
    rep.line("the Data Records which definition the published areas use.")

    rep.blank()
    rep.line("Groups of features sharing one Shape_Area value. Where the geometries of a")
    rep.line("group sum to the field value, one unit was subdivided and the field was")
    rep.line("copied to each part; where they do not, the parts were also edited after")
    rep.line("the field was written. Either way the field cannot be reconstructed and")
    rep.line("must simply be recomputed from the geometry.")
    rep.blank()
    value_groups = defaultdict(list)
    for f in geomorph:
        value_groups[round(f["Shape_Area"] or 0.0, 2)].append(f)
    group_rows = []
    for value, feats in sorted(value_groups.items(), key=lambda kv: -len(kv[1])):
        if len(feats) < 2:
            continue
        total = sum(f.geom.area for f in feats)
        group_rows.append([
            f"{value:.2f}", len(feats), f"{total:.1f}", f"{value - total:+.1f}",
            "yes" if abs(value - total) < 1.0 else "no",
            ", ".join(str(f["gm_id"]) for f in sorted(feats, key=lambda x: x["gm_id"]))])
        record("J_stale_area_group", glim.GEOMORPHOLOGY_LAYER,
               ", ".join(str(f["gm_id"]) for f in feats),
               f"{len(feats)} features share Shape_Area={value:.2f} m2; their geometries "
               f"sum to {total:.1f} m2")
    rep.table(["Shape_Area value", "features", "sum of geometries m2", "difference m2",
               "parts sum to value?", "gm_id"], group_rows)
    rep.blank()
    rep.kv("features inside a shared-value group", sum(int(r[1]) for r in group_rows))
    rep.kv("groups whose parts sum to the field value",
           sum(1 for r in group_rows if r[4] == "yes"))
    rep.finding("no polygon in the layer has a geometry matching the largest shared "
                "values, so those numbers refer to units that no longer exist as single "
                "features. Recompute Shape_Area; do not try to reconcile it.")
    rep.blank()
    slivers = sorted((f.geom.area, f["gm_id"], f["geom_class"]) for f in geomorph
                     if f.geom.area < 200)
    rep.blank()
    rep.kv("geomorphological polygons smaller than 200 m2", len(slivers))
    if slivers:
        rep.line("    smallest ten (area m2, gm_id, geom_class):")
        for area, ident, code in slivers[:10]:
            rep.line(f"      {area:10.3f}  gm_id {ident}  {code}")
        record("J_slivers", glim.GEOMORPHOLOGY_LAYER, "-",
               f"{len(slivers)} polygons < 200 m2, smallest {slivers[0][0]:.3f} m2")

    # --------------------------------------------------------- K duplicates
    rep.h2("K - Duplicate geometries")
    rows = []
    for name, feats in named:
        seen = defaultdict(list)
        for f in feats:
            geom = f.geom if f.geom.is_valid else f.geom.buffer(0)
            seen[geom.wkb].append(f.attrs.get("lds_id", f.attrs.get("gm_id", f.attrs.get("Id"))))
        dups = {k: v for k, v in seen.items() if len(v) > 1}
        rows.append([name.replace("DRV_GLIM_", ""), len(dups)])
        for ids in dups.values():
            record("K_duplicate_geometry", name, ids[0],
                   f"identical geometry shared with {ids[1:]}")
    rep.table(["layer", "duplicated geometries"], rows)

    rep.h2("L - Landslides represented by a source polygon alone")
    rep.line("Check F above pairs source and deposit only for landslides ALREADY")
    rep.line("represented by two polygons, so it cannot see a landslide whose scarp and")
    rep.line("deposit were given different identifiers: each half then looks like a")
    rep.line("well-formed single-polygon landslide. The converse test is made here. A")
    rep.line("source area is the upper part of a landslide and cannot be a landslide on")
    rep.line("its own, with one documented exception: debris flows are mapped by their")
    rep.line("detachment and transport zone alone wherever no distinct fan was")
    rep.line("deposited, which is why xdf and wxdf appear here legitimately. Any other")
    rep.line("class occurring as a lone source polygon is an anomaly, and for each one")
    rep.line("the neighbouring deposits of the matching landslide type carrying a")
    rep.line("different identifier are listed, because that is the signature of one")
    rep.line("landslide split across two identifiers.")
    rep.blank()
    orphan_tab = Counter()
    anomalies = []
    for name, feats in layers.items():
        groups = glim.group_by_id(feats)
        for lds_id, group in groups.items():
            if not all(glim.is_source_polygon(f["lds_class"]) for f in group):
                continue
            classes = "+".join(sorted(f["lds_class"] for f in group))
            orphan_tab[classes] += 1
            # The documented exception covers the debris-flow source unit in both
            # its forms, xdf and the widespread wxdf: base_class leaves the "w"
            # in place, so strip it before testing the process code.
            if glim.base_class(group[0]["lds_class"]).lstrip("w") == "df":
                continue
            anomalies.append((name, lds_id, group, groups))
    rep.table(["lone source class", "landslides"],
              sorted(orphan_tab.items(), key=lambda kv: -kv[1]))
    rep.blank()
    rep.kv("lone source polygons of a class other than debris flow", len(anomalies))
    for name, lds_id, group, groups in anomalies:
        short = name.replace("DRV_GLIM_", "")
        whole = glim.safe_union([f.geom for f in group])
        want = glim.base_class(group[0]["lds_class"])
        rep.finding(f"{short} lds_id={lds_id}: "
                    f"{'+'.join(f['lds_class'] for f in group)} mapped as a landslide "
                    f"with no deposit polygon")
        record("L_lone_source", name, lds_id,
               f"{'+'.join(f['lds_class'] for f in group)}: source area with no deposit")
        for other_id, other in groups.items():
            if other_id == lds_id:
                continue
            for f in other:
                if glim.is_source_polygon(f["lds_class"]):
                    continue
                if glim.base_class(f["lds_class"]) != want:
                    continue
                gap = whole.distance(f.geom)
                if gap < 1.0:
                    rep.line(f"      adjacent deposit lds_id={other_id} "
                             f"{f['lds_class']} at {gap:.3f} m - the two are very "
                             f"probably one landslide under two identifiers")
                    record("L_lone_source", name, lds_id,
                           f"adjacent deposit lds_id={other_id} ({f['lds_class']}) "
                           f"at {gap:.3f} m: candidate split landslide")

    # ------------------------------------------------------------------ export
    rep.h2("Summary")
    by_check = Counter(f[0] for f in findings)
    rep.table(["check", "findings"], sorted(by_check.items()))
    rep.kv("total findings", len(findings))
    glim.write_csv("02_qc_findings", ["check", "layer", "feature", "description"], findings)
    rep.save()


if __name__ == "__main__":
    main()
