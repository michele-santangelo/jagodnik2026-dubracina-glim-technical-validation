#!/usr/bin/env python3
"""
12 - The three-way overlap with the previous inventory, as a map and as areas.

Why this exists. The index of agreement between this inventory and Jagodnik et
al. (2020) is 0.156, which looks low against the 95.6 % of their polygons that
this inventory recovers. The two numbers measure different things and both are
correct: recovery asks how much of THEIR map is inside ours, while the index of
agreement is symmetric and is driven down by everything we map that they do not.
This script draws the partition so the asymmetry can be inspected rather than
argued about, and reports the area of each part.

Over the analysis domain - their pilot area clipped to the study boundary - the
ground falls into three classes plus the unmapped remainder:

    both          mapped as landslide by both inventories
    only_new      mapped only by the present G-LIM
    only_previous mapped only by Jagodnik et al. (2020)
    neither       mapped by neither

Two partitions are produced: one using every generation of the present
inventory, and one using the old class alone (O1-O4), which is the class the
previous inventory targeted.

Output: outputs/12_overlap_map.gpkg   the partition, for inspection in a GIS
        outputs/12_overlap_map.png    a two-panel map
        outputs/12_overlap_areas.csv  the areas
        outputs/12_overlap_map.txt    the report
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from osgeo import ogr, osr

import glim
from glim import Report

# Categorical slots 1-3 of the validated palette. Checked with
# validate_palette.js --pairs all (worst CVD dE 9.2, normal-vision 24.0): the
# classes are spatially interleaved, so every pair must separate, not only
# adjacent ones. The aqua warns on contrast against the surface, which the
# legend labels and the area table relieve.
COL = {"both": "#2a78d6", "only_new": "#eb6834", "only_previous": "#1baf7a"}
LABEL = {"both": "mapped by both",
         "only_new": "present inventory only",
         "only_previous": "Jagodnik et al. 2020 only"}
ORDER = ["both", "only_new", "only_previous"]


def patches(geom, **kw):
    """A matplotlib patch per polygon, holes included."""
    out = []
    if geom.geom_type == "Polygon":
        polys = [geom]
    else:
        # A difference can return a GeometryCollection holding slivers that
        # degenerated to lines or points; only the polygons are drawable.
        polys = [g for g in geom.geoms if g.geom_type == "Polygon"]
    for poly in polys:
        if poly.is_empty:
            continue
        verts, codes = [], []
        for ring in [poly.exterior] + list(poly.interiors):
            pts = list(ring.coords)
            verts.extend(pts)
            codes.extend([Path.MOVETO] + [Path.LINETO] * (len(pts) - 2) + [Path.CLOSEPOLY])
        out.append(PathPatch(Path(verts, codes), **kw))
    return out


def partition(mine, theirs, domain):
    return {"both": mine.intersection(theirs),
            "only_new": mine.difference(theirs),
            "only_previous": theirs.difference(mine),
            "neither": domain.difference(glim.safe_union([mine, theirs]))}


def main():
    rep = Report("12_overlap_map",
                 "12 - The three-way overlap with the previous inventory")
    glim.require_jom(rep)
    boundary = glim.read_boundary()
    layers = glim.read_all_landslides(fix=True)
    jom = glim.read_jom_polygons(fix=True)
    domain = glim.read_jom_domain().intersection(boundary)
    theirs = glim.safe_union([f.geom for v in jom.values() for f in v]).intersection(domain)

    old = [n for n in layers if glim.RELATIVE_AGE[n] == "O"]
    rockfalls = glim.read_rockfalls(fix=True)
    # The first variant must reproduce the areal matrix of script 04, which
    # includes the rockfall layer in the union; the second is the old class
    # alone, without rockfalls, as in section 1b of that script. Mixing the two
    # conventions would put two different numbers on the same quantity.
    variants = [
        ("all generations + rockfalls", sorted(layers), True),
        ("old class O1-O4 + rockfalls", old, True),
        ("old class O1-O4 only", old, False),
    ]
    parts, rows = {}, []
    for label, names, with_rf in variants:
        geoms = [f.geom for n in names for f in layers[n]]
        if with_rf:
            geoms += [f.geom for f in rockfalls]
        mine = glim.safe_union(geoms).intersection(domain)
        p = partition(mine, theirs, domain)
        parts[label] = p
        inter, union = p["both"].area, mine.union(theirs).area
        for k in ORDER + ["neither"]:
            rows.append([label, k, f"{glim.km2(p[k].area):.4f}",
                         f"{100 * p[k].area / domain.area:.1f}"])
        rows.append([label, "index of agreement", f"{inter / union:.3f}", ""])
        rows.append([label, "their area recovered %",
                     f"{100 * inter / theirs.area:.1f}", ""])
        rows.append([label, "our area also theirs %",
                     f"{100 * inter / mine.area:.1f}", ""])

    rep.h2("Areas of each part of the partition")
    rep.line(f"Analysis domain: {glim.km2(domain.area):.4f} km2, the pilot area of")
    rep.line("the previous inventory clipped to the study boundary. Areas are planar")
    rep.line("on EPSG:3765; ratios are unaffected by the basis (see script 03).")
    rep.blank()
    rep.table(["subset of the present inventory", "class", "km2 / value",
               "% of domain"], rows)
    rep.blank()
    rep.line("Why the index of agreement is low while 95.6 % of their polygons are")
    rep.line("recovered. The two measure different things. Recovery asks what share of")
    rep.line("THEIR polygons this inventory covers, and answers 95.6 % of the polygons")
    rep.line("and 63.0 % of their area. The index of agreement is symmetric: it divides")
    rep.line("the shared area by the union, and the union is dominated by the 8.20 km2")
    rep.line("this inventory maps and the previous one does not. The index therefore")
    rep.line("measures how much has been added, not how far the two disagree where both")
    rep.line("looked. The ground neither maps, 5.59 km2, does not enter it at all.")
    rep.blank()
    rep.line("A caution on the restricted comparison. Dropping to the old class alone")
    rep.line("removes the rockfall layer from OUR side while leaving the rock-fall")
    rep.line("chutes and cliff sectors on THEIRS, which is why the ground mapped only")
    rep.line("by the previous inventory grows from 1.00 to 1.77 km2 - visibly, along")
    rep.line("the north-eastern carbonate margin. The like-for-like restriction is the")
    rep.line("old class WITH the rockfall layer, and on that subset the index of")
    rep.line("agreement is the highest of the three, because removing the very old and")
    rep.line("relict polygons removes area that is ours alone without removing")
    rep.line("agreement. That is the same result the increment analysis reports from")
    rep.line("the other direction.")

    # ------------------------------------------------------------------ GPKG
    path = os.path.join(glim.OUTPUTS, "12_overlap_map.gpkg")
    if os.path.exists(path):
        os.remove(path)
    drv = ogr.GetDriverByName("GPKG")
    ds = drv.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(glim.EXPECTED_EPSG)
    for label, p in parts.items():
        name = (label.replace(" + ", "_plus_").replace("-", "_")
                     .replace(" ", "_").lower())
        lyr = ds.CreateLayer(name, srs, ogr.wkbMultiPolygon)
        if lyr is None:
            raise RuntimeError(f"could not create layer {name!r} in {path}")
        lyr.CreateField(ogr.FieldDefn("class", ogr.OFTString))
        lyr.CreateField(ogr.FieldDefn("area_km2", ogr.OFTReal))
        for k in ORDER + ["neither"]:
            if p[k].is_empty:
                continue
            f = ogr.Feature(lyr.GetLayerDefn())
            f.SetField("class", k)
            f.SetField("area_km2", round(glim.km2(p[k].area), 4))
            f.SetGeometry(ogr.CreateGeometryFromWkb(p[k].wkb))
            lyr.CreateFeature(f)
            f = None
    del ds
    rep.blank()
    rep.kv("map for inspection in a GIS", glim.shown_path(path))

    # ------------------------------------------------------------------- PNG
    fig, axes = plt.subplots(1, 3, figsize=(19.5, 7.6))
    for ax, (label, _, _) in zip(axes, variants):
        p = parts[label]
        for pt in patches(domain, facecolor="#f0efec", edgecolor="none", zorder=0):
            ax.add_patch(pt)
        for k in ORDER:
            for pt in patches(p[k], facecolor=COL[k], edgecolor="none",
                              alpha=0.95, zorder=2):
                ax.add_patch(pt)
        for pt in patches(domain, facecolor="none", edgecolor="#5b5750",
                          linewidth=1.1, zorder=3):
            ax.add_patch(pt)
        ax.set_title(label, fontsize=12, loc="left", color="#1b1a17", pad=10)
        sub = "  ·  ".join(f"{LABEL[k]} {glim.km2(p[k].area):.2f} km²" for k in ORDER)
        ax.set_xlabel(sub, fontsize=8, color="#45423c", labelpad=8)
        x0, y0, x1, y1 = domain.bounds
        m = 0.02 * max(x1 - x0, y1 - y0)
        ax.set_xlim(x0 - m, x1 + m)
        ax.set_ylim(y0 - m, y1 + m)
        ax.set_aspect("equal")
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
    handles = [Line2D([0], [0], marker="s", linestyle="none", markersize=11,
                      markerfacecolor=COL[k], markeredgecolor="none", label=LABEL[k])
               for k in ORDER]
    handles.append(Line2D([0], [0], marker="s", linestyle="none", markersize=11,
                          markerfacecolor="#f0efec", markeredgecolor="#5b5750",
                          label="analysis domain, mapped by neither"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Overlap with the inventory of Jagodnik et al. (2020) over the shared domain",
                 fontsize=13.5, x=0.055, ha="left", y=0.965, color="#1b1a17")
    fig.tight_layout(rect=[0, 0.055, 1, 0.94])
    png = os.path.join(glim.OUTPUTS, "12_overlap_map.png")
    fig.savefig(png, dpi=170, facecolor="#fcfcfb")
    plt.close(fig)
    rep.kv("map as an image", glim.shown_path(png))

    glim.write_csv("12_overlap_areas",
                   ["subset", "class", "value", "percent_of_domain"], rows)
    rep.save()


if __name__ == "__main__":
    main()
