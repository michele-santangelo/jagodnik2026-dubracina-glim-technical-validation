#!/usr/bin/env python3
"""
01 - Descriptive statistics of the G-LIM, layer by layer.

How many polygons, how many *landslides* (unique `lds_id`), how many are
represented by one and by two polygons, and how the counts break down by
relative age, by type, by uncertainty code and by `depth`. The counts quoted in
the Data Records section of the paper are traceable to this report.

Output: outputs/01_inventory_stats.txt, outputs/01_landslide_counts.csv
"""

from collections import Counter

import glim
from glim import Report


def main():
    rep = Report("01_inventory_stats", "01 - G-LIM descriptive statistics")
    layers = glim.read_all_landslides(fix=True)

    # ---------------------------------------------------------------- per layer
    rep.h2("Polygons and landslides per layer")
    rows = []
    totals = Counter()
    per_layer_ids = {}
    for name, feats in layers.items():
        grouped = glim.group_by_id(feats)
        per_layer_ids[name] = set(grouped)
        counts = Counter(len(v) for v in grouped.values())
        area_field = sum(f["Shape_Area"] or 0.0 for f in feats)
        area_geom = sum(f.geom.area for f in feats)
        rows.append([
            name.replace("DRV_GLIM_", ""),
            glim.RELATIVE_AGE[name],
            glim.GENERATION[name] if glim.GENERATION[name] else "-",
            len(feats),
            len(grouped),
            counts[1],
            counts[2],
            sum(v for k, v in counts.items() if k > 2),
            f"{glim.km2(area_field):.4f}",
            f"{glim.km2(area_geom):.4f}",
        ])
        totals["polygons"] += len(feats)
        totals["landslides"] += len(grouped)
        totals["one_polygon"] += counts[1]
        totals["two_polygons"] += counts[2]
        totals["more_than_two"] += sum(v for k, v in counts.items() if k > 2)
        totals["area_field"] += area_field
        totals["area_geom"] += area_geom
    rep.table(
        ["layer", "age", "gen", "polygons", "landslides", "1-poly", "2-poly",
         ">2-poly", "sum Shape_Area km2", "sum geometry km2"],
        rows + [["TOTAL", "", "", totals["polygons"], totals["landslides"],
                 totals["one_polygon"], totals["two_polygons"],
                 totals["more_than_two"],
                 f"{glim.km2(totals['area_field']):.4f}",
                 f"{glim.km2(totals['area_geom']):.4f}"]],
    )
    rep.line()
    rep.line("'landslides' counts distinct lds_id values *within* each layer.")
    rep.line("Note that lds_id is only unique within a layer (see 02_qc_database),")
    rep.line("so the total is the sum of the per-layer counts, not a global distinct count.")

    # ------------------------------------------------------ widespread subtotal
    rep.h2("Widespread (clustered / coalescent) process units")
    widespread = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            if any(glim.is_widespread(f["lds_class"]) for f in group):
                widespread.append([name.replace("DRV_GLIM_", ""), lds_id,
                                   "+".join(f["lds_class"] for f in group),
                                   f"{sum(f.geom.area for f in group):.1f}"])
    rep.table(["layer", "lds_id", "lds_class", "area m2"], widespread)
    n_total = totals["landslides"]
    n_wide = len(widespread)
    rep.line()
    rep.kv("landslides, all", n_total)
    rep.kv("landslides, widespread units only", n_wide)
    rep.kv("landslides, excluding widespread units", n_total - n_wide)

    # ---------------------------------------------------- counts by age class
    rep.h2("Counts by relative-age class")
    by_age = Counter()
    by_age_nowide = Counter()
    for name, feats in layers.items():
        age = glim.RELATIVE_AGE[name]
        for lds_id, group in glim.group_by_id(feats).items():
            by_age[age] += 1
            if not any(glim.is_widespread(f["lds_class"]) for f in group):
                by_age_nowide[age] += 1
    rep.table(["relative age", "landslides", "excl. widespread"],
              [[a, by_age[a], by_age_nowide[a]] for a in ("O", "VO", "VOR")]
              + [["TOTAL", sum(by_age.values()), sum(by_age_nowide.values())]])
    rep.line()
    rep.line("The manuscript reports 1,071 O + 20 VO + 3 VOR = 1,094 landslides.")
    if sum(by_age.values()) != 1094:
        rep.finding(
            f"the data give {sum(by_age.values())} landslides including widespread units "
            f"and {sum(by_age_nowide.values())} excluding them; neither equals 1,094, "
            "so the counts in the text predate the current shapefiles and must be recomputed")

    # ------------------------------------------------------------ type, uncert
    rep.h2("lds_class frequency (polygons)")
    codes = Counter(f["lds_class"] for feats in layers.values() for f in feats)
    rep.table(["lds_class", "polygons", "source area?", "widespread?"],
              [[c, n, "yes" if glim.is_source_polygon(c) else "no",
                "yes" if glim.is_widespread(c) else "no"]
               for c, n in sorted(codes.items())])

    rep.h2("Landslide type per landslide (from the deposit polygon)")
    types = Counter()
    for feats in layers.values():
        for group in glim.group_by_id(feats).values():
            deposit = [f for f in group if not glim.is_source_polygon(f["lds_class"])]
            code = deposit[0]["lds_class"] if deposit else \
                glim.base_class(group[0]["lds_class"])
            types[code] += 1
    rep.table(["type", "landslides", "%"],
              [[c, n, f"{100 * n / sum(types.values()):.1f}"]
               for c, n in types.most_common()])

    rep.h2("uncert code frequency (polygons)")
    unc = Counter(f["uncert"] for feats in layers.values() for f in feats)
    rep.table(["uncert", "detection", "delineation", "classification", "polygons"],
              [[u, u[0], u[1], u[2], n] for u, n in sorted(unc.items())])

    rep.h2("depth frequency")
    dep_poly = Counter(f["depth"] for feats in layers.values() for f in feats)
    dep_lds = Counter()
    for feats in layers.values():
        for group in glim.group_by_id(feats).values():
            dep_lds[group[0]["depth"]] += 1
    rep.table(["depth", "polygons", "landslides"],
              [[d, dep_poly.get(d, 0), dep_lds.get(d, 0)]
               for d in sorted(set(dep_poly) | set(dep_lds), key=str)])

    # -------------------------------------------------------- ancillary layers
    rep.h2("Rockfall and geomorphological layers")
    rf = glim.read_rockfalls(fix=True)
    rf_codes = Counter(f["lds_class"] for f in rf)
    rep.kv("rockfall polygons", len(rf))
    rep.kv("rockfall cumulative area (km2)", f"{glim.km2(sum(f.geom.area for f in rf)):.4f}")
    rep.table(["lds_class", "polygons"], sorted(rf_codes.items()))

    gm = glim.read_layer(glim.ORIGINAL_DATA, glim.GEOMORPHOLOGY_LAYER,
                         ["gm_id", "geom_class", "Shape_Area"], fix=True)
    gm_codes = Counter(f["geom_class"] for f in gm)
    rep.blank()
    rep.kv("geomorphological polygons", len(gm))
    rep.kv("geomorphology cumulative area (km2)",
           f"{glim.km2(sum(f.geom.area for f in gm)):.4f}")
    rep.table(["geom_class", "polygons", "area km2"],
              [[c, n, f"{glim.km2(sum(f.geom.area for f in gm if f['geom_class'] == c)):.4f}"]
               for c, n in sorted(gm_codes.items())])

    # ------------------------------------------------------------------ export
    csv_rows = []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            deposit = [f for f in group if not glim.is_source_polygon(f["lds_class"])]
            source = [f for f in group if glim.is_source_polygon(f["lds_class"])]
            csv_rows.append([
                name, glim.RELATIVE_AGE[name], glim.GENERATION[name] or "", lds_id,
                len(group),
                deposit[0]["lds_class"] if deposit else "",
                source[0]["lds_class"] if source else "",
                group[0]["depth"], group[0]["uncert"],
                round(sum(f.geom.area for f in source), 2),
                round(sum(f.geom.area for f in deposit), 2),
                round(sum(f.geom.area for f in group), 2),
            ])
    glim.write_csv(
        "01_landslide_counts",
        ["layer", "relative_age", "generation", "lds_id", "n_polygons",
         "deposit_class", "source_class", "depth", "uncert",
         "source_area_m2", "deposit_area_m2", "total_area_m2"],
        csv_rows)
    rep.blank()
    rep.line("outputs/01_landslide_counts.csv holds one record per landslide with")
    rep.line("source, deposit and total area, so that the area of a landslide does not")
    rep.line("have to be computed by summing Shape_Area over its polygons.")
    rep.save()


if __name__ == "__main__":
    main()
