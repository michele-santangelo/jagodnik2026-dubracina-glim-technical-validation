"""
Shared I/O, configuration and helpers for the Dubracina G-LIM technical validation.

All analyses in this folder import from here so that layer lists, attribute
dictionaries and geometry-handling conventions are defined exactly once.

Conventions
-----------
* All vector data are read with GDAL/OGR and converted to shapely geometries.
* Geometries are NEVER silently repaired when the purpose of the check is to
  find broken geometries (see `read_layer(fix=False)`, used by 02_qc_database).
  Where a valid geometry is a precondition for an areal computation (unions,
  intersections) `buffer(0)` is applied and the fact is reported.
* Areas are computed from the geometry in the layer CRS (EPSG:3765,
  HTRS96 / Croatia TM, metre) and are therefore in square metres.
  The `Shape_Area` attribute is treated as *data to be checked*, never as
  the source of truth.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import warnings
from dataclasses import dataclass

# The Debian GDAL python bindings are built against numpy 1.x; importing
# osgeo.gdal_array under numpy 2.x prints a spurious warning. We do not use
# gdal_array (no raster arrays are read), so the warning is suppressed here.
warnings.filterwarnings("ignore", category=FutureWarning, module="osgeo")

from osgeo import ogr, osr  # noqa: E402
from shapely import wkb as shapely_wkb  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
#: the dataset under analysis: the layers of the Zenodo data record, unpacked
#: into data/. Override with GLIM_DATA to point elsewhere.
ORIGINAL_DATA = os.environ.get("GLIM_DATA", os.path.join(HERE, "data"))
#: the inventory of Jagodnik et al. (2020), which is NOT public and is not
#: distributed here (see README). Override with GLIM_JOM to point elsewhere.
PREVIOUS_INVENTORIES = os.environ.get(
    "GLIM_JOM", os.path.join(ORIGINAL_DATA, "jagodnik2020"))
#: where the reports are written; override with GLIM_OUTPUTS so that a run over
#: another data state does not overwrite the results of the current one
OUTPUTS = os.environ.get("GLIM_OUTPUTS", os.path.join(HERE, "outputs"))


def shown_path(path):
    """A path as printed in a report: relative to this folder when inside it."""
    path = os.path.abspath(path)
    if os.path.commonpath([path, HERE]) == HERE:
        return os.path.relpath(path, HERE)
    return path

# --------------------------------------------------------------------------
# Dataset definition (must mirror Table 6 of the manuscript)
# --------------------------------------------------------------------------

#: landslide layers, in the relative-age / generation order of Table 6
LANDSLIDE_LAYERS = [
    "DRV_GLIM_O_I",
    "DRV_GLIM_O_II",
    "DRV_GLIM_O_III",
    "DRV_GLIM_O_IV",
    "DRV_GLIM_VO_I",
    "DRV_GLIM_VO_II",
    "DRV_GLIM_VO_R",
]
ROCKFALL_LAYER = "DRV_GLIM_rockfalls"
GEOMORPHOLOGY_LAYER = "DRV_GLIM_Geomorphology"
BOUNDARY_LAYER = "DRV_Boundary"
#: The contour layer is ancillary: it is not part of the landslide inventory,
#: it was not produced by the interpretation this paper describes, and it is
#: distributed as it stands. It is therefore no longer validated and no longer
#: enters the data-state digest, so that the digest tracks only the data the
#: results depend on.
CONTOUR_LAYER = "DRV_contour5"

#: relative-age class of each landslide layer
RELATIVE_AGE = {
    "DRV_GLIM_O_I": "O",
    "DRV_GLIM_O_II": "O",
    "DRV_GLIM_O_III": "O",
    "DRV_GLIM_O_IV": "O",
    "DRV_GLIM_VO_I": "VO",
    "DRV_GLIM_VO_II": "VO",
    "DRV_GLIM_VO_R": "VOR",
}
#: generation number within the relative-age class (VOR has none)
GENERATION = {
    "DRV_GLIM_O_I": 1,
    "DRV_GLIM_O_II": 2,
    "DRV_GLIM_O_III": 3,
    "DRV_GLIM_O_IV": 4,
    "DRV_GLIM_VO_I": 1,
    "DRV_GLIM_VO_II": 2,
    "DRV_GLIM_VO_R": None,
}

LANDSLIDE_FIELDS = ["lds_id", "lds_class", "uncert", "depth", "Shape_Area"]

#: Attributes added after submission, absent from the data as first delivered.
#: They are read when the layer carries them and set to None when it does not,
#: so that the same checks run over both states of the dataset.
LANDSLIDE_OPTIONAL = ["lds_area", "cluster_id"]

#: allowed values of `lds_class` as documented in Table 2 of the manuscript.
#: "(x)" in the table means the code exists with and without the `x` prefix
#: that marks the scarp/source-area polygon.
TABLE2_LDS_CLASS = {
    "ds", "xds",      # debris slide
    "s", "xs",        # earth slide
    "rsl", "xrsl",    # rock slide
    "sef", "xsef",    # earth slide-flow
    "df", "xdf",      # debris flow
    "f", "xf",        # earth flow
    "rav", "xrav",    # rock avalanche
    "rf", "xrf",      # rock fall
    "wxrf",           # widespread rock-fall source area
    "wxdf",           # widespread debris-flow source unit
    "ws",             # widespread earth slide
}
# Table 2 as revised. Two corrections the checks themselves produced: `wrf` was
# documented but occurs nowhere in the data - the widespread rock-fall polygons
# are all source areas on the cliffs, the matching deposits being the talus of
# the geomorphological layer - and `wxdf` occurs in the data but was missing
# from the table. Both are now stated as they are, so neither is a finding.

#: allowed values of `geom_class` as documented in Table 5 of the manuscript
#: `lake` was present in the data and missing from the published table; the
#: revised Table 5 lists it (the Tribalj reservoir, one polygon).
TABLE5_GEOM_CLASS = {"ad", "af", "oaf", "t", "ot", "gl", "lake"}

#: allowed values of `depth` (Table 7)
TABLE7_DEPTH = {"deep-seated", "shallow"}

#: Process types that are shallow whatever their polygon structure. A debris
#: flow is mapped with a separate upslope polygon wherever the detachment and
#: transport zone can be told from the fan, but that polygon is a concave
#: erosional zone and not the scarp of a failure surface, so it implies nothing
#: about the depth of one. The submitted text stated the two-polygon rule
#: without this exception; the check below encodes the rule as the Methods now
#: state it, which is why the 79 debris flows it used to flag are no longer
#: findings. They were never defects in the data.
ALWAYS_SHALLOW = {"df"}


def expected_depth(group):
    """The depth class the stated criteria require for one landslide.

    `group` is the list of polygons sharing an lds_id. Returns None when the
    criteria do not determine a class.
    """
    classes = {base_class(f["lds_class"]).lstrip("w") for f in group}
    if classes & ALWAYS_SHALLOW:
        return "shallow"
    if any(is_source_polygon(f["lds_class"]) for f in group):
        return "deep-seated"
    if len(group) == 1:
        return "shallow"
    return None

#: study-area figure quoted in the submitted manuscript, km2. It is the
#: hydrological basin area (43.55 km2 in Tosevski 2013), not the mapped domain.
#: Kept so that reports can show the old and the new figure side by side.
AREA_QUOTED_KM2 = 43.6

#: the study area the paper quotes: the measured area of DRV_Boundary, which is
#: the domain that was actually mapped. Every density and coverage percentage in the paper follows from
#: this value, so it is defined once here. Ellipsoidal, to match the basis the
#: dataset's own Shape_Area now uses (see AREA_BASIS below); the planar grid
#: value is 38.4159 km2, and both round to 38.4 km2.
AREA_STUDY_KM2 = 38.4043
AREA_STUDY_GRID_KM2 = 38.4159

#: How the areas stored in Shape_Area were measured, for the Data Records.
#: Established from the data and confirmed against the QGIS project settings:
#: QGIS `$area` with ellipsoidal measurement enabled, so the
#: coordinates are transformed to the chosen ellipsoid and the area is returned
#: in square metres on that ellipsoid rather than on the projection plane.
AREA_BASIS = {
    "software": "QGIS, field calculator, $area",
    "measurement": "ellipsoidal (project ellipsoidal measurement enabled)",
    "ellipsoid": "GRS 1980 (EPSG:7019)",
    "layer_crs": "HTRS96 / Croatia TM (EPSG:3765)",
    "units": "square metres",
    "note": ("grid areas computed from the projected coordinates are about "
             "300 ppm larger, because HTRS96/TM has k0 = 0.9999 on the central "
             "meridian at 16.5 deg E and the valley lies about 145 km west of it"),
}

#: cartographic completeness threshold quoted in the Usage Notes, m2
COMPLETENESS_THRESHOLD_M2 = 400.0

#: EPSG code the whole dataset is expected to use
EXPECTED_EPSG = 3765

# --------------------------------------------------------------------------
# Previous inventories
# --------------------------------------------------------------------------

#: Jagodnik et al. (2020, J. Maps) engineering-geological inventory, shared as
#: one shapefile per mapped landform. `JOM_Debris flow_thalweg_htrs` is a line
#: layer (thalwegs) and is excluded from all areal comparisons.
JOM_LINE_LAYERS = ["JOM_Debris flow_thalweg_htrs"]

#: The mapping domain of that inventory, supplied by its authors.
#: It shares the JOM_ prefix but is NOT an inventory layer: including it among
#: the landform layers would make the inventory's union equal to the whole
#: pilot area. Having the real domain replaces the convex-hull proxy the areal
#: comparison had to use before, so "not mapped" and "outside the mapped area"
#: can finally be told apart.
JOM_DOMAIN_LAYER = "JOM_C_Pilot area_htrs"

#: expansion of the abbreviations used in the JOM file names (from READ_ME.docx)
JOM_LABELS = {
    "JOM_Debris flow_thalweg_htrs": "debris-flow thalwegs (lines)",
    "JOM_Debris slide_Flow_htrs": "debris slide-flow",
    "JOM_Debris slide_htrs": "debris slide",
    "JOM_Fall_Chutes_htrs": "rock-fall related chutes",
    "JOM_RFall_Cliff_htrs": "rock-fall source area along cliffs",
    "JOM_RFall_Slope_htrs": "rock-fall source area along slopes",
    "JOM_RIrregSlide_htrs": "rock irregular slide",
    "JOM_RRotatSlide_htrs": "rock rotational slide",
    "JOM_RSDeformation_htrs": "rock slope deformation",
    "JOM_RotatSlide_htrs": "rotational slide",
}

#: CRS of the coordinates published in Table 1 of Tosevski (2013):
#: "MGI Balkans 5" = MGI 1901 / Balkans zone 5, Gauss-Krueger with central
#: meridian 15 deg E and false easting 5 500 000 m.
TOSEVSKI_EPSG = 31275
TOSEVSKI_TABLE = os.path.join(HERE, "tosevski2013_inventory.csv")

# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------


@dataclass
class Feature:
    """One vector feature: its attributes (as a dict) and its shapely geometry."""

    attrs: dict
    geom: object
    layer: str
    fid: int

    def __getitem__(self, key):
        return self.attrs[key]


def read_layer(directory: str, name: str, fields=(), fix: bool = False,
               optional=()):
    """Read a shapefile into a list of `Feature`.

    Parameters
    ----------
    directory : str
        folder holding the shapefile
    name : str
        layer name, with or without the ``.shp`` extension
    fields : iterable of str
        attributes to load; missing attributes raise, so this doubles as a
        schema check
    fix : bool
        if True, invalid geometries are repaired with ``buffer(0)``. Leave
        False whenever the analysis is meant to *detect* invalid geometries.
    optional : iterable of str
        attributes loaded when the layer carries them and set to None when it
        does not; unlike `fields` they are not a schema check.
    """
    path = os.path.join(directory, name if name.endswith(".shp") else name + ".shp")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    ds = ogr.Open(path)
    if ds is None:
        raise IOError(f"OGR could not open {path}")
    layer = ds.GetLayer()
    schema = {layer.GetLayerDefn().GetFieldDefn(i).GetName()
              for i in range(layer.GetLayerDefn().GetFieldCount())}
    missing = [f for f in fields if f not in schema]
    if missing:
        raise KeyError(f"{name}: missing attribute(s) {missing}; present: {sorted(schema)}")
    out = []
    for feature in layer:
        ref = feature.GetGeometryRef()
        geom = shapely_wkb.loads(bytes(ref.ExportToWkb()))
        if fix and not geom.is_valid:
            geom = geom.buffer(0)
        attrs = {f: feature.GetField(f) for f in fields}
        attrs.update({f: (feature.GetField(f) if f in schema else None)
                      for f in optional})
        out.append(Feature(attrs, geom, os.path.basename(path), feature.GetFID()))
    del ds  # release the datasource only after the geometries have been copied
    return out


def layer_epsg(directory: str, name: str):
    """Return the EPSG code declared by the layer's .prj, or None."""
    path = os.path.join(directory, name if name.endswith(".shp") else name + ".shp")
    ds = ogr.Open(path)
    srs = ds.GetLayer().GetSpatialRef()
    del ds
    if srs is None:
        return None
    srs.AutoIdentifyEPSG()
    code = srs.GetAuthorityCode(None)
    return int(code) if code else None


def layer_fields(directory: str, name: str):
    """Return the attribute names a layer actually carries."""
    path = os.path.join(directory, name if name.endswith(".shp") else name + ".shp")
    ds = ogr.Open(path)
    defn = ds.GetLayer().GetLayerDefn()
    names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
    del ds
    return names


def read_rockfalls(fix: bool = True):
    """Read the rockfall layer, whichever identifier field it carries.

    The layer originally had a single `Id` field holding 0 for all 486 records,
    which gave its polygons no usable identifier. It was later given a real
    `lds_id` (and a `Shape_Area`), so both schemas have to be readable: the
    reports must keep working on an archived copy of either version.
    Whatever the file holds, the returned features expose `lds_id`.
    """
    present = layer_fields(ORIGINAL_DATA, ROCKFALL_LAYER)
    wanted = [f for f in ("lds_id", "Id", "lds_class", "Shape_Area", "uncert", "depth")
              if f in present]
    feats = read_layer(ORIGINAL_DATA, ROCKFALL_LAYER, wanted, fix=fix)
    for feature in feats:
        if "lds_id" not in feature.attrs:
            feature.attrs["lds_id"] = feature.attrs.get("Id")
        feature.attrs.setdefault("Shape_Area", None)
    return feats


def read_boundary():
    """Return the study-area polygon as a single shapely geometry."""
    feats = read_layer(ORIGINAL_DATA, BOUNDARY_LAYER, ["Id", "area_km2"], fix=True)
    if len(feats) != 1:
        raise ValueError(f"{BOUNDARY_LAYER} holds {len(feats)} features, expected 1")
    return feats[0].geom


def read_all_landslides(fix: bool = True):
    """Read every landslide layer, tagging each feature with its layer name."""
    out = {}
    for name in LANDSLIDE_LAYERS:
        out[name] = read_layer(ORIGINAL_DATA, name, LANDSLIDE_FIELDS,
                               fix=fix, optional=LANDSLIDE_OPTIONAL)
    return out


def jom_available():
    """True when the Jagodnik et al. (2020) inventory is where GLIM_JOM points."""
    return os.path.exists(os.path.join(PREVIOUS_INVENTORIES, JOM_DOMAIN_LAYER + ".shp"))


def require_jom(rep):
    """End a script that cannot run without the Jagodnik et al. (2020) data.

    That inventory is not public. Without it the script writes a report saying
    so and exits normally, so that `run_all.sh` carries on with the others.
    """
    if jom_available():
        return
    rep.line("SKIPPED: the inventory of Jagodnik et al. (2020) was not found in")
    rep.line(f"  {shown_path(PREVIOUS_INVENTORIES)}")
    rep.line("Those data are not public and are not distributed with this code;")
    rep.line("they must be requested from the authors of that work (see README).")
    rep.save()
    sys.exit(0)


def read_jom_polygons(fix: bool = True):
    """Read the polygon layers of the JOM (2020) inventory, keyed by layer name."""
    out = {}
    for filename in sorted(os.listdir(PREVIOUS_INVENTORIES)):
        if not filename.endswith(".shp"):
            continue
        name = filename[:-4]
        if name in JOM_LINE_LAYERS or name == JOM_DOMAIN_LAYER:
            continue
        out[name] = read_layer(PREVIOUS_INVENTORIES, filename, fix=fix)
    return out


def read_jom_domain():
    """Return the mapping domain of the Jagodnik et al. (2020) inventory."""
    feats = read_layer(PREVIOUS_INVENTORIES, JOM_DOMAIN_LAYER + ".shp", fix=True)
    return safe_union([f.geom for f in feats])


def transformer(from_epsg: int, to_epsg: int = EXPECTED_EPSG):
    """Return an (x, y) -> (x, y) coordinate transformation.

    Axis order is forced to traditional GIS order (easting, northing) so that
    the call signature does not depend on the CRS definition.
    """
    src = osr.SpatialReference()
    src.ImportFromEPSG(from_epsg)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(to_epsg)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(src, dst)

    def _tr(x, y):
        px, py, _ = ct.TransformPoint(x, y)
        return px, py

    _tr.source_wkt = src.ExportToWkt()
    return _tr


# --------------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------------


#: True when pyproj is available, so geodesic areas can be computed
try:
    from pyproj import Geod, Transformer  # noqa: E402
    _TO_GEOGRAPHIC = Transformer.from_crs(f"EPSG:{EXPECTED_EPSG}", "EPSG:4326",
                                          always_xy=True)
    _GEOD = Geod(ellps="GRS80")
    HAVE_PYPROJ = True
except Exception:  # pragma: no cover - the grid area is still reported
    HAVE_PYPROJ = False


def geodesic_area(geom):
    """Area on the GRS80 ellipsoid, in square metres.

    This is not the same number as `geom.area`, which is the planar area of the
    projected coordinates. HTRS96 / Croatia TM has a scale factor of 0.9999 on
    the central meridian at 16.5 deg E, and the Dubracina valley sits about
    145 km west of it, where the grid scale has grown back above unity. Areas
    measured on the grid are therefore about 0.03 % larger than areas on the
    ellipsoid, and which of the two a GIS reports depends on the tool: ArcGIS
    "Calculate Geometry" offers both, and its geodesic option is what produced
    the Shape_Area values in this dataset.

    The difference is far below the mapping uncertainty at 1:10,000, but it has
    to be stated, because otherwise a recomputed field looks like a corrupted
    one.
    """
    if not HAVE_PYPROJ:
        return None
    from shapely.ops import transform as _transform
    geographic = _transform(lambda x, y: _TO_GEOGRAPHIC.transform(x, y), geom)
    total = 0.0
    for part in parts(geographic):
        area, _ = _GEOD.geometry_area_perimeter(part)
        total += abs(area)
    return total


def safe_union(geoms):
    """Dissolve geometries, repairing invalid ones first (required by GEOS)."""
    return unary_union([g if g.is_valid else g.buffer(0) for g in geoms])


def iou(a, b):
    """Intersection over union (Jaccard index) of two polygons."""
    inter = a.intersection(b).area
    if inter <= 0:
        return 0.0
    return inter / (a.area + b.area - inter)


def parts(geom):
    """Return the list of single-part polygons composing `geom`."""
    return list(geom.geoms) if geom.geom_type.startswith("Multi") else [geom]


def is_source_polygon(lds_class):
    """True if the code marks a scarp/source-area polygon (the `x` prefix)."""
    if not lds_class:
        return False
    return lds_class.startswith("x") or lds_class.startswith("wx")


def base_class(lds_class):
    """Strip the source-area prefix: 'xs' -> 's', 'wxdf' -> 'wdf'."""
    if not lds_class:
        return lds_class
    if lds_class.startswith("wx"):
        return "w" + lds_class[2:]
    if lds_class.startswith("x"):
        return lds_class[1:]
    return lds_class


def is_widespread(lds_class):
    """True for the 'widespread' (clustered/coalescent) process units."""
    return bool(lds_class) and lds_class.startswith("w")


def stable_key(feature):
    """Return a delete-proof identifier for one polygon.

    Shapefile feature ids are POSITIONAL: deleting a record renumbers every
    record after it, so an fid noted today can point at a different polygon
    tomorrow. Never use fid as the key of a correction worklist. These
    attribute keys survive editing:

      landslide layers    (lds_id, lds_class) - verified unique within each layer
      geomorphological    gm_id               - unique, 0 to 294
      rockfall layer      none, because Id is 0 for all 486 records; the
                          centroid is the only locator

    Returns a short string, or "-" when the layer offers no attribute key.
    """
    attrs = feature.attrs
    if "lds_id" in attrs and "lds_class" in attrs:
        return f"{attrs['lds_id']}/{attrs['lds_class']}"
    if "gm_id" in attrs:
        return str(attrs["gm_id"])
    if "ELEV" in attrs:
        return f"ELEV {attrs['ELEV']:.0f}"
    return "-"


def centroid_xy(feature):
    """Return the centroid as (easting, northing) rounded to the metre.

    A position is the one locator that no editing operation invalidates, so it
    is carried in every findings table alongside the attribute key.
    """
    geom = feature.geom if feature.geom.is_valid else feature.geom.buffer(0)
    point = geom.centroid
    if point.is_empty:
        point = geom.representative_point()
    return round(point.x), round(point.y)


def group_by_id(features):
    """Group landslide features by `lds_id`, preserving insertion order."""
    grouped = {}
    for f in features:
        grouped.setdefault(f["lds_id"], []).append(f)
    return grouped


# --------------------------------------------------------------------------
# Data provenance
# --------------------------------------------------------------------------

#: the files of the deposited dataset whose content the results depend on.
#: The Jagodnik et al. (2020) layers are deliberately left out, so that anyone
#: who downloads the data record obtains the same digest, with or without them.
def _input_files():
    names = (LANDSLIDE_LAYERS + [ROCKFALL_LAYER, GEOMORPHOLOGY_LAYER,
                                 BOUNDARY_LAYER])
    for name in names:
        for extension in (".shp", ".dbf"):
            path = os.path.join(ORIGINAL_DATA, name + extension)
            if os.path.exists(path):
                yield name + extension, path


def data_state():
    """Fingerprint the input data so every report says which version it describes.

    Results are comparable only when they were computed from the same data, so
    every report records which data it describes.

    The digest is a short hash over the SHA-256 content hashes of the input
    files, so it depends on their content only, not on where or when they were
    unpacked, and it can be compared with the one in the reference outputs.
    """
    entries = []
    for label, path in _input_files():
        stat = os.stat(path)
        entries.append({"file": label, "sha256": sha256(path),
                        "mtime": time.strftime("%Y-%m-%d %H:%M:%S",
                                               time.localtime(stat.st_mtime))})
    blob = json.dumps([[e["file"], e["sha256"]] for e in entries]).encode()
    return {
        "digest": hashlib.sha256(blob).hexdigest()[:12],
        "files": entries,
        "newest": max((e["mtime"] for e in entries), default="-"),
    }


def sha256(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


class Report:
    """Collects the lines of an analysis and writes them to stdout and a file.

    Every script produces one plain-text report under `outputs/` so that the
    numbers quoted in the manuscript can be traced back to a dated run.
    """

    def __init__(self, name: str, title: str):
        self.name = name
        self.lines = []
        os.makedirs(OUTPUTS, exist_ok=True)
        self.h1(title)
        self.kv("python", sys.version.split()[0])
        from osgeo import gdal
        self.kv("GDAL", gdal.__version__)
        import shapely
        self.kv("shapely", shapely.__version__)
        self.kv("data", shown_path(ORIGINAL_DATA))
        self.kv("Jagodnik et al. (2020) data",
                "found" if jom_available() else "not available")
        state = data_state()
        self.kv("data state digest", state["digest"])
        self.line()
        self.line("The data state digest identifies the exact version of the input")
        self.line("shapefiles these numbers describe. If it differs between two reports,")
        self.line("the data changed in between and the numbers are not comparable.")
        self.line("Run 00_data_fingerprint.py for the per-file content hashes.")
        self.blank()

    def _emit(self, text=""):
        print(text)
        self.lines.append(text)

    def h1(self, text):
        self._emit("=" * 78)
        self._emit(text)
        self._emit("=" * 78)

    def h2(self, text):
        self.blank()
        self._emit(text)
        self._emit("-" * len(text))

    def line(self, text=""):
        self._emit(text)

    def blank(self):
        self._emit("")

    def kv(self, key, value):
        self._emit(f"{key:<38s} {value}")

    def table(self, header, rows):
        widths = [len(str(h)) for h in header]
        rows = [[str(c) for c in r] for r in rows]
        for r in rows:
            widths = [max(w, len(c)) for w, c in zip(widths, r)]
        fmt = "  ".join("{:<%d}" % w for w in widths)
        self._emit(fmt.format(*[str(h) for h in header]))
        self._emit(fmt.format(*["-" * w for w in widths]))
        for r in rows:
            self._emit(fmt.format(*r))

    def finding(self, text):
        """Flag a discrepancy between the data and their description."""
        self._emit(f"  [!] {text}")

    def save(self):
        path = os.path.join(OUTPUTS, self.name + ".txt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(self.lines) + "\n")
        print(f"\n-> written {shown_path(path)}")
        return path


def write_csv(name, header, rows):
    """Write a result table to `outputs/<name>.csv`."""
    import csv

    os.makedirs(OUTPUTS, exist_ok=True)
    path = os.path.join(OUTPUTS, name + ".csv")
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"-> written {shown_path(path)}")
    return path


def km2(area_m2):
    return area_m2 / 1e6
