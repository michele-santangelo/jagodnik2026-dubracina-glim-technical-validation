#!/usr/bin/env python3
"""
04 - Comparison with the previous LiDAR-based inventory (Jagodnik et al. 2020, J. Maps).

Independent validation against a pre-existing inventory of the same area.

The inventory of Jagodnik et al. (2020) is not public. It was made available to
the authors of this work for the comparison, and it is not distributed with
this code; without it this script reports that it was skipped (see README).

Three levels of agreement are quantified, following Galli et al. (2008)
"Comparing landslide inventory maps":

  1. DETECTION - is each previously mapped landform also mapped as a landslide
     in the new G-LIM? Measured as the fraction of each previous polygon that
     the new inventory covers.
  2. AREAL AGREEMENT - the error index / confusion matrix over a common
     analysis domain: area mapped by both, by only one, by neither.
  3. BOUNDARY AND CLASS AGREEMENT - for every previous polygon, the best
     matching new polygon by intersection over union, with its landslide type
     and the layer it belongs to, giving a type cross-tabulation.
  4. CONFIRMATION AND INCREMENT - the converse question. Of the landslides the
     present inventory maps over the shared domain, which had already been
     mapped by the previous one and which are new, and what the new ones are
     made of. Agreement is not by itself the objective of the comparison: an
     inventory that reproduced the previous one exactly would add no knowledge.

The previous inventory classifies its features by type of movement, by material
and by estimated depth (Hungr et al. 2014). It assigns NO relative age and no
state of activity - not in the published paper, whose Table 2 is a typology, and
not in the delivered shapefiles, whose only descriptive fields are the
topographic location (TopoLok), whether the feature is whole or partial
(Cjelovitos) and the flank aspect (Bok_VD). No comparison of relative age
between the two inventories is therefore possible, and none is attempted here.

The mapping domain of the previous inventory was supplied by its authors as
`JOM_C_Pilot area_htrs`, so the areal comparison is run over the real domain and "not mapped" can be told apart from "outside the mapped
area". The convex hull that stood in for it before is still computed, for the
record, because it was 15 % larger and therefore diluted the agreement.

Output: outputs/04_compare_jom.txt, outputs/04_jom_feature_match.csv,
        outputs/04_jom_type_crosstab.csv
"""

from collections import Counter, defaultdict

from shapely.strtree import STRtree

import glim
from glim import Report


def main():
    rep = Report("04_compare_jom",
                 "04 - Comparison with Jagodnik et al. (2020, J. Maps) inventory")
    glim.require_jom(rep)
    boundary = glim.read_boundary()
    layers = glim.read_all_landslides(fix=True)
    rockfalls = glim.read_rockfalls(fix=True)
    jom = glim.read_jom_polygons(fix=True)

    # -------------------------------------------------------------- inventory
    rep.h2("Content of the previous inventory")
    rows = []
    jom_geoms = []
    for name, feats in sorted(jom.items()):
        n_parts = sum(len(glim.parts(f.geom)) for f in feats)
        area = sum(f.geom.area for f in feats)
        rows.append([name.replace("JOM_", "").replace("_htrs", ""),
                     glim.JOM_LABELS.get(name, "?"), len(feats), n_parts,
                     f"{glim.km2(area):.4f}"])
        jom_geoms.extend(f.geom for f in feats)
    rep.table(["layer", "landform (READ_ME.docx)", "features", "single parts",
               "cumulative km2"], rows)
    jom_union = glim.safe_union(jom_geoms)
    rep.blank()
    rep.kv("polygon features, total", len(jom_geoms))
    rep.kv("cumulative area (km2)", f"{glim.km2(sum(g.area for g in jom_geoms)):.4f}")
    rep.kv("union area (km2)", f"{glim.km2(jom_union.area):.4f}")
    rep.kv("internal overlap (km2)",
           f"{glim.km2(sum(g.area for g in jom_geoms) - jom_union.area):.4f}")
    rep.kv("union inside the G-LIM boundary (km2)",
           f"{glim.km2(jom_union.intersection(boundary).area):.4f}")
    rep.kv("union outside the G-LIM boundary (km2)",
           f"{glim.km2(jom_union.difference(boundary).area):.4f}")
    for name in glim.JOM_LINE_LAYERS:
        rep.line()
        rep.line(f"{name} is a line layer (debris-flow thalwegs) and is excluded "
                 "from all areal comparisons.")

    # -------------------------------------------------------- analysis domain
    rep.h2("Analysis domain")
    pilot = glim.read_jom_domain()
    domain = pilot.intersection(boundary)
    hull = jom_union.convex_hull.intersection(boundary)
    rep.table(
        ["domain", "km2", "% of the study area"],
        [["pilot area of the previous inventory, as supplied",
          f"{glim.km2(pilot.area):.4f}", "-"],
         ["pilot area clipped to the study boundary = analysis domain",
          f"{glim.km2(domain.area):.4f}",
          f"{100 * domain.area / boundary.area:.1f}"],
         ["pilot area falling outside the study boundary",
          f"{glim.km2(pilot.difference(boundary).area):.4f}", "-"],
         ["convex-hull proxy used before the domain was available",
          f"{glim.km2(hull.area):.4f}",
          f"{100 * hull.area / boundary.area:.1f}"]])
    rep.line()
    rep.line(f"The proxy was {hull.area / domain.area:.2f} times the real domain, so every")
    rep.line("agreement index computed against it was diluted. The figures below use the")
    rep.line("real domain.")
    outside = sum(1 for feats in jom.values() for f in feats
                  if not pilot.buffer(1).contains(f.geom))
    rep.blank()
    rep.kv("previously mapped polygons not fully inside the pilot area", outside)
    rep.kv("previous inventory area inside the pilot area",
           f"{100 * jom_union.intersection(pilot).area / jom_union.area:.1f} %")

    # ------------------------------------------------------- 1. detection
    rep.h2("1 - Detection agreement")
    rep.line("For each previously mapped polygon, the fraction of its area that the")
    rep.line("new G-LIM (landslide layers plus rockfall layer) covers.")
    rep.blank()
    new_geoms = [f.geom for feats in layers.values() for f in feats] + \
                [f.geom for f in rockfalls]
    new_union = glim.safe_union(new_geoms)
    rows = []
    detection = {}
    for name, feats in sorted(jom.items()):
        fractions = [f.geom.intersection(new_union).area / f.geom.area for f in feats]
        detection[name] = fractions
        rows.append([name.replace("JOM_", "").replace("_htrs", ""), len(feats),
                     f"{100 * sum(fractions) / len(fractions):.1f}",
                     sum(1 for x in fractions if x >= 0.5),
                     sum(1 for x in fractions if x > 0.0),
                     sum(1 for x in fractions if x == 0.0)])
    all_fractions = [x for v in detection.values() for x in v]
    rows.append(["ALL", len(all_fractions),
                 f"{100 * sum(all_fractions) / len(all_fractions):.1f}",
                 sum(1 for x in all_fractions if x >= 0.5),
                 sum(1 for x in all_fractions if x > 0.0),
                 sum(1 for x in all_fractions if x == 0.0)])
    rep.table(["layer", "features", "mean coverage %", ">=50 % covered",
               "any overlap", "no overlap"], rows)
    n_hit = sum(1 for x in all_fractions if x >= 0.5)
    rep.blank()
    rep.kv("previously mapped polygons recovered (>=50 % covered)",
           f"{n_hit} / {len(all_fractions)} ({100 * n_hit / len(all_fractions):.1f} %)")

    # Which generations carry the recovery? The very old and relict polygons are
    # large enough to overlap much of the valley, so a recovery figure computed
    # over every layer at once could be an artefact of that overlap. The same
    # figures are therefore recomputed over subsets of the present inventory.
    rep.h2("1b - Which generations of the present inventory carry the recovery")
    rep.line("The previous inventory targeted currently active, gully-related")
    rep.line("processes, which belong to the old class here. Restricting the present")
    rep.line("inventory to that class tests whether the recovery above depends on the")
    rep.line("large very old and relict polygons.")
    rep.blank()
    old_class = [n for n in layers if glim.RELATIVE_AGE[n] == "O"]
    jom_feats = [f for feats in jom.values() for f in feats]
    theirs_dom = jom_union.intersection(domain)
    subsets = [("all generations + rockfalls", sorted(layers), True),
               ("all generations, no rockfalls", sorted(layers), False),
               ("old class O1-O4 + rockfalls", old_class, True),
               ("old class O1-O4 only", old_class, False)]
    rows = []
    for label, names, with_rf in subsets:
        geoms_sub = [f.geom for n in names for f in layers[n]]
        if with_rf:
            geoms_sub += [f.geom for f in rockfalls]
        union_sub = glim.safe_union(geoms_sub)
        fr = [f.geom.intersection(union_sub).area / f.geom.area for f in jom_feats]
        hit = sum(1 for x in fr if x >= 0.5)
        mine_sub = union_sub.intersection(domain)
        inter = mine_sub.intersection(theirs_dom).area
        rows.append([label, f"{100 * sum(fr) / len(fr):.1f}",
                     f"{hit} / {len(fr)}", f"{100 * hit / len(fr):.1f}",
                     f"{100 * inter / theirs_dom.area:.1f}",
                     f"{inter / mine_sub.union(theirs_dom).area:.3f}"])
    rep.table(["subset of the present inventory", "mean coverage %",
               ">=50 % covered", "%", "their area recovered %", "IoU"], rows)
    rep.line()
    rep.line("Both figures are needed. The recovery falls when the older")
    rep.line("generations are excluded, which is expected, but the index of agreement")
    rep.line("does not: the older generations add mapped area on our side without")
    rep.line("adding agreement, which is precisely the increment quantified in 4.")

    # -------------------------------------------------- 2. areal agreement
    rep.h2("2 - Areal agreement over the analysis domain")
    mine = new_union.intersection(domain)
    theirs = jom_union.intersection(domain)
    both = mine.intersection(theirs)
    only_mine = mine.difference(theirs)
    only_theirs = theirs.difference(mine)
    neither = domain.difference(glim.safe_union([mine, theirs]))
    rep.table(
        ["class", "km2", "% of domain"],
        [["landslide in both inventories", f"{glim.km2(both.area):.4f}",
          f"{100 * both.area / domain.area:.1f}"],
         ["landslide only in the new G-LIM", f"{glim.km2(only_mine.area):.4f}",
          f"{100 * only_mine.area / domain.area:.1f}"],
         ["landslide only in the previous inventory", f"{glim.km2(only_theirs.area):.4f}",
          f"{100 * only_theirs.area / domain.area:.1f}"],
         ["landslide in neither", f"{glim.km2(neither.area):.4f}",
          f"{100 * neither.area / domain.area:.1f}"],
         ["DOMAIN", f"{glim.km2(domain.area):.4f}", "100.0"]])
    union_area = glim.safe_union([mine, theirs]).area
    rep.blank()
    rep.kv("new G-LIM landslide area in the domain (km2)", f"{glim.km2(mine.area):.4f}")
    rep.kv("previous inventory area in the domain (km2)", f"{glim.km2(theirs.area):.4f}")
    rep.kv("error index E = 1 - intersection/union", f"{1 - both.area / union_area:.3f}")
    rep.kv("index of agreement (Jaccard) = intersection/union",
           f"{both.area / union_area:.3f}")
    rep.kv("previous inventory area recovered by the new G-LIM",
           f"{100 * both.area / theirs.area:.1f} %")
    rep.kv("new G-LIM area also present in the previous inventory",
           f"{100 * both.area / mine.area:.1f} %")
    rep.line()
    rep.line("The asymmetry is expected and is the substance of the comparison: the")
    rep.line("previous inventory mapped active, mostly gully-related instabilities,")
    rep.line("whereas the new G-LIM also maps old and very old deep-seated landslides.")

    # ---------------------------------------- 3. best match, type and age
    rep.h2("3 - Best matching new polygon per previous polygon (by IoU)")
    records = []
    for name, feats in layers.items():
        for f in feats:
            records.append((glim.RELATIVE_AGE[name], glim.GENERATION[name],
                            name.replace("DRV_GLIM_", ""),
                            f["lds_id"], f["lds_class"], f["depth"], f.geom))
    for f in rockfalls:
        records.append(("-", None, "rockfalls", f["lds_id"], f["lds_class"], "-", f.geom))
    geoms = [r[6] for r in records]
    tree = STRtree(geoms)

    match_rows = []
    crosstab = defaultdict(Counter)
    age_tab = defaultdict(Counter)
    iou_by_layer = defaultdict(list)
    for name, feats in sorted(jom.items()):
        for f in feats:
            best, best_iou = None, 0.0
            for idx in tree.query(f.geom):
                value = glim.iou(f.geom, geoms[idx])
                if value > best_iou:
                    best_iou, best = value, records[idx]
            label = name.replace("JOM_", "").replace("_htrs", "")
            iou_by_layer[label].append(best_iou)
            coverage = f.geom.intersection(new_union).area / f.geom.area
            if best is None:
                crosstab[label]["<no overlap>"] += 1
                age_tab[label]["<no overlap>"] += 1
                match_rows.append([name, f.fid, round(f.geom.area, 1),
                                   round(coverage, 3), "", "", "", "", 0.0])
            else:
                crosstab[label][best[4]] += 1
                age_tab[label][best[2]] += 1
                match_rows.append([name, f.fid, round(f.geom.area, 1),
                                   round(coverage, 3), best[2], best[3], best[4],
                                   best[5], round(best_iou, 3)])

    rows = []
    for label, values in iou_by_layer.items():
        values_sorted = sorted(values)
        rows.append([label, len(values),
                     f"{sum(values) / len(values):.2f}",
                     f"{values_sorted[len(values) // 2]:.2f}",
                     sum(1 for x in values if x >= 0.5),
                     sum(1 for x in values if x >= 0.3)])
    rep.table(["previous layer", "features", "mean IoU", "median IoU",
               "IoU>=0.5", "IoU>=0.3"], rows)
    rep.line()
    rep.line("A caveat on interpretation: the very-old-relict layer contains a few")
    rep.line("polygons of several km2 which overlap almost everything, so a criterion")
    rep.line("based on the largest overlapping area would be biased towards them.")
    rep.line("Intersection over union is used instead, which selects the polygon of")
    rep.line("comparable size and position.")

    # ------------------------------------- the error index of Carrara et al.
    rep.h2("3b - Error index E of Carrara et al. (1992), pairwise")
    rep.line("E = (A u B - A n B) / (A u B) for two corresponding polygons, 0 for a")
    rep.line("perfect match and 1 for none. It is the complement of the intersection")
    rep.line("over union, E = 1 - IoU, so the two carry the same information; what")
    rep.line("differs is where it is applied.")
    rep.blank()
    rep.line("Applied PAIRWISE, as Carrara et al. intend and as Santangelo et al.")
    rep.line("(2015) apply it, E measures the discrepancy between a previously mapped")
    rep.line("polygon and its counterpart here. Applied GLOBALLY to the two unions")
    rep.line("over the domain it measures something else - how much of the combined")
    rep.line("mapped area is not shared - and it is then dominated by whatever one")
    rep.line("inventory maps and the other does not. The global value is reported in")
    rep.line("section 2 and in script 12; the pairwise distribution is reported here.")
    rep.line("A polygon with no counterpart has E = 1 by definition and is kept in,")
    rep.line("since dropping it would flatter the result.")
    rep.blank()

    def quantile(values, q):
        v = sorted(values)
        if not v:
            return float("nan")
        i = q * (len(v) - 1)
        lo, hi = int(i), min(int(i) + 1, len(v) - 1)
        return v[lo] + (v[hi] - v[lo]) * (i - lo)

    rows = []
    for label, values in iou_by_layer.items():
        e = [1 - x for x in values]
        rows.append([label, len(e), f"{sum(e) / len(e):.2f}",
                     f"{quantile(e, 0.25):.2f}", f"{quantile(e, 0.5):.2f}",
                     f"{quantile(e, 0.75):.2f}",
                     sum(1 for x in e if x <= 0.5)])
    all_e = [1 - x for v in iou_by_layer.values() for x in v]
    rows.append(["ALL", len(all_e), f"{sum(all_e) / len(all_e):.2f}",
                 f"{quantile(all_e, 0.25):.2f}", f"{quantile(all_e, 0.5):.2f}",
                 f"{quantile(all_e, 0.75):.2f}",
                 sum(1 for x in all_e if x <= 0.5)])
    rep.table(["previous layer", "features", "mean E", "Q1", "median E", "Q3",
               "E<=0.5"], rows)

    # Carrara's index is known to depend on feature size: a small polygon is
    # penalised by a positional shift that a large one absorbs. Reporting E
    # without this is reporting half the result.
    rep.blank()
    rep.line("E against the size of the previously mapped polygon. The index is known")
    rep.line("to penalise small features, whose E a positional shift of a few metres")
    rep.line("moves far more than that of a large one, so the distribution is given by")
    rep.line("size class rather than as one number.")
    rep.blank()
    sized = []
    for name, feats in sorted(jom.items()):
        label = name.replace("JOM_", "").replace("_htrs", "")
        for f, value in zip(feats, iou_by_layer[label]):
            sized.append((f.geom.area, 1 - value))
    bands = [(0, 500), (500, 2000), (2000, 10000), (10000, float("inf"))]
    names = ["< 500", "500 - 2,000", "2,000 - 10,000", "> 10,000"]
    rows = []
    for (lo, hi), band in zip(bands, names):
        sub = [e for a, e in sized if lo <= a < hi]
        if not sub:
            continue
        rows.append([band, len(sub), f"{sum(sub) / len(sub):.2f}",
                     f"{quantile(sub, 0.5):.2f}",
                     sum(1 for x in sub if x <= 0.5)])
    rep.table(["previous polygon area, m2", "features", "mean E", "median E",
               "E<=0.5"], rows)

    rep.h2("Landslide type of the best match, per previous layer")
    for label in crosstab:
        total = sum(crosstab[label].values())
        rep.line(f"  {label} (n={total}): " +
                 ", ".join(f"{k}={v}" for k, v in crosstab[label].most_common(6)))
    rep.h2("Layer of the best match, per previous layer")
    rep.line("This is NOT an agreement on relative age: the previous inventory")
    rep.line("assigns none. It says which generation of the present inventory the")
    rep.line("previously mapped features fall in, and is a coherence check internal")
    rep.line("to this inventory, not an external validation of its age attribution.")
    rep.blank()
    for label in age_tab:
        total = sum(age_tab[label].values())
        rep.line(f"  {label} (n={total}): " +
                 ", ".join(f"{k}={v}" for k, v in age_tab[label].most_common(6)))
    dated_layers = {n.replace("DRV_GLIM_", "") for n in glim.RELATIVE_AGE}
    counts = Counter()
    for label in age_tab:
        counts.update(age_tab[label])
    dated = sum(v for k, v in counts.items() if k in dated_layers)
    old_cls = sum(v for k, v in counts.items()
                  if k in dated_layers
                  and glim.RELATIVE_AGE["DRV_GLIM_" + k] == "O")
    undated = sum(counts.values()) - dated
    rep.blank()
    rep.kv("previous features whose best match carries a generation", dated)
    rep.kv("  of which matched to the old class O1-O4",
           f"{old_cls} / {dated} ({100 * old_cls / dated:.1f} %)")
    rep.kv("previous features with no datable match", undated)
    rep.line()
    rep.line("The features with no datable match are those whose best match lies in")
    rep.line("the rockfall layer, which is not subdivided by generation, plus those")
    rep.line("with no counterpart at all. They cannot enter the numerator and are")
    rep.line("excluded from the denominator; computing the share over all")
    rep.line(f"{sum(counts.values())} previous features instead would deflate it to "
             f"{100 * old_cls / sum(counts.values()):.1f} %")
    rep.line("for a reason that has nothing to do with agreement.")

    # ------------------------------------------ 4. confirmation and increment
    rep.h2("4 - Confirmation and increment over the shared domain")
    rep.line("The converse of section 1. For every landslide of the present inventory")
    rep.line("that reaches the shared domain, the fraction of its own area - clipped")
    rep.line("to that domain - which the previous inventory covers. A landslide")
    rep.line("covered over at least half of its area was already mapped; one below")
    rep.line("that threshold is an addition. The threshold is a convention, so the")
    rep.line("count with no overlap at all is reported as well, as a lower bound that")
    rep.line("does not depend on it.")
    rep.blank()
    confirmed, added = [], []
    for name, feats in layers.items():
        for lds_id, group in glim.group_by_id(feats).items():
            whole = glim.safe_union([f.geom for f in group])
            if not whole.intersects(domain):
                continue
            clipped = whole.intersection(domain)
            if clipped.area <= 0:
                continue
            covered = clipped.intersection(jom_union).area / clipped.area
            rec = {"layer": name, "age": glim.RELATIVE_AGE[name], "lds_id": lds_id,
                   "area": clipped.area, "covered": covered,
                   "depth": group[0]["depth"], "geom": clipped}
            (confirmed if covered >= 0.5 else added).append(rec)
    everything = confirmed + added

    def median(values):
        """Standard median: the mean of the two central values when n is even,
        so that the figure matches what a reader recomputes from the CSV."""
        values = sorted(values)
        if not values:
            return 0.0
        mid = len(values) // 2
        if len(values) % 2:
            return values[mid]
        return 0.5 * (values[mid - 1] + values[mid])

    union_conf = glim.safe_union([r["geom"] for r in confirmed])
    union_added = glim.safe_union([r["geom"] for r in added])
    rep.table(
        ["class", "landslides", "cumulative km2", "footprint km2", "median m2"],
        [["already mapped by the previous inventory (>=50 % covered)",
          len(confirmed), f"{glim.km2(sum(r['area'] for r in confirmed)):.4f}",
          f"{glim.km2(union_conf.area):.4f}",
          f"{median([r['area'] for r in confirmed]):.1f}"],
         ["added by the present inventory (<50 % covered)",
          len(added), f"{glim.km2(sum(r['area'] for r in added)):.4f}",
          f"{glim.km2(union_added.area):.4f}",
          f"{median([r['area'] for r in added]):.1f}"],
         ["  of which no overlap at all",
          sum(1 for r in added if r["covered"] == 0.0),
          f"{glim.km2(sum(r['area'] for r in added if r['covered'] == 0.0)):.4f}",
          "-",
          f"{median([r['area'] for r in added if r['covered'] == 0.0]):.1f}"]])
    rep.blank()
    rep.kv("footprint of the additions outside the confirmed landslides (km2)",
           f"{glim.km2(union_added.difference(union_conf).area):.4f}")
    rep.kv("footprint of the previous inventory in the domain (km2)",
           f"{glim.km2(theirs_dom.area):.4f}")
    rep.line()
    rep.line("Cumulative and footprint areas are both given because polygons of")
    rep.line("different generations overlap: the cumulative value is not a land area.")

    rep.h2("What the increment is made of - by size")
    bands = [(0, 1000), (1000, 5000), (5000, 20000),
             (20000, 100000), (100000, float("inf"))]
    labels = ["< 1,000", "1,000 - 5,000", "5,000 - 20,000",
              "20,000 - 100,000", "> 100,000"]
    rows = []
    for (lo, hi), label in zip(bands, labels):
        sub = [r for r in everything if lo <= r["area"] < hi]
        if not sub:
            continue
        n_add = sum(1 for r in sub if r["covered"] < 0.5)
        rows.append([label, len(sub), len(sub) - n_add, n_add,
                     f"{100 * n_add / len(sub):.1f}"])
    rep.table(["landslide area, m2", "landslides", "already mapped", "added",
               "% added"], rows)
    biggest = max(confirmed, key=lambda r: r["area"])
    rep.blank()
    rep.kv("largest landslide the previous inventory had already mapped (m2)",
           f"{biggest['area']:.0f}")
    rep.line("Jagodnik et al. (2020) report 49,462 m2 as the largest polygon of their")
    rep.line("inventory, so this is the same feature: it is the detection ceiling of")
    rep.line("the previous mapping, reached independently from both sides.")

    rep.h2("What the increment is made of - by generation and by depth")
    rows = []
    for age in ("O", "VO", "VOR"):
        sub = [r for r in everything if r["age"] == age]
        if not sub:
            continue
        n_add = sum(1 for r in sub if r["covered"] < 0.5)
        rows.append([age, len(sub), len(sub) - n_add, n_add,
                     f"{100 * n_add / len(sub):.1f}",
                     f"{glim.km2(sum(r['area'] for r in sub if r['covered'] < 0.5)):.4f}"])
    rep.table(["relative age", "landslides", "already mapped", "added", "% added",
               "cumulative km2 added"], rows)
    rep.blank()
    for label, recs in (("already mapped", confirmed), ("added", added)):
        deep = sum(1 for r in recs if r["depth"] == "deep-seated")
        rep.kv(f"{label}: deep-seated",
               f"{deep} / {len(recs)} ({100 * deep / len(recs):.1f} %)")
    rep.line()
    rep.line("Depth class does not separate the two groups; size and generation do.")
    rep.line("The increment is not a deeper reading of the same landslides, it is a")
    rep.line("different population: large, old, subdued features whose expression is")
    rep.line("an arrangement of slope facets over hundreds of metres rather than the")
    rep.line("sharp morphology of a recent failure. That is what the pseudo-")
    rep.line("stereoscopic interpretation is for, and it is why no external dataset")
    rep.line("can validate this part of the inventory: for this class of feature the")
    rep.line("method used here is more sensitive than the alternatives available to")
    rep.line("check it against.")

    glim.write_csv("04_increment_by_landslide",
                   ["layer", "relative_age", "lds_id", "area_in_domain_m2",
                    "fraction_covered_by_previous", "depth", "verdict"],
                   [[r["layer"], r["age"], r["lds_id"], round(r["area"], 1),
                     round(r["covered"], 3), r["depth"],
                     "already mapped" if r["covered"] >= 0.5 else "added"]
                    for r in sorted(everything,
                                    key=lambda r: -r["area"])])

    glim.write_csv("04_jom_feature_match",
                   ["previous_layer", "previous_fid", "previous_area_m2",
                    "fraction_covered_by_new", "matched_layer", "matched_lds_id",
                    "matched_lds_class", "matched_depth", "iou"], match_rows)
    ct_codes = sorted({c for v in crosstab.values() for c in v})
    glim.write_csv("04_jom_type_crosstab", ["previous_layer"] + ct_codes,
                   [[label] + [crosstab[label].get(c, 0) for c in ct_codes]
                    for label in crosstab])
    rep.save()


if __name__ == "__main__":
    main()
