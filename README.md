# Dubračina Valley G-LIM — technical validation code

Code for the quality control and technical validation of the geomorphological
landslide inventory map (G-LIM) of the Dubračina Valley, Croatia, described in

> Jagodnik, P., Fiorucci, F., Bernat Gazibara, S. & Santangelo, M.
> Geomorphological landslide inventory of the Dubračina Valley using
> pseudo-stereoscopic LiDAR hillshade maps. *Scientific Data* (2026).
> https://doi.org/XXXXXXX

The inventory itself is published separately, under CC BY 4.0:

> Jagodnik, P., Fiorucci, F., Bernat Gazibara, S. & Santangelo, M.
> Geomorphological landslide inventory of the Dubračina Valley using pseudo-stereoscopic LiDAR hillshade maps. Version 2. Zenodo. https://doi.org/10.5281/zenodo.22814494

Every number in the Technical Validation section of the paper is produced by
these scripts from the published data, with a single command.

## What each script does

| Script | Purpose |
|---|---|
| `glim.py` | Shared I/O, layer lists, attribute dictionaries (Tables 2, 5 and 7 of the paper), geometry helpers, ellipsoidal areas, data fingerprinting, report writer |
| `00_data_fingerprint.py` | SHA-256 content hash, size and feature count of every input file, plus one digest for the whole data state |
| `01_inventory_stats.py` | Polygon and landslide counts per layer, relative age, type, uncertainty code and depth class; one record per landslide with source, deposit and total area |
| `02_qc_database.py` | CRS consistency, geometry validity, empty and multipart geometries, missing values, values outside the published dictionaries, identifier structure, source/deposit pairing, attribute consistency, the depth rule, features outside the boundary, `Shape_Area` against the geometry, duplicate geometries |
| `03_area_statistics.py` | Cumulative polygon area against the union footprint, overlap between generations, stacking depth, size distribution against the 400 m² threshold, density |
| `04_compare_jom.py` | Comparison with Jagodnik et al. (2020): detection agreement, areal agreement and error index (Galli et al. 2008), best match by intersection over union, increment. **Needs non-public data, see below** |
| `05_compare_tosevski.py` | Comparison with Toševski (2013): transcription checks, reprojection of the published coordinates, point-in-polygon test, metadata comparison |
| `06_depth_rule_records.py` | Record-level listing of the landslides that break the shallow / deep-seated rule, with the depth-by-type cross-tabulation |
| `07_qc_topology.py` | Record-level geometric and topological QC: duplicates and near-duplicates, degenerate and sliver polygons, vertex defects, overlaps within a generation, source/deposit topology, holes, overlaps in the geomorphological layer |
| `09_validation_tables.py` | The tables of the Technical Validation section: integrity checks (flagged / corrected / remaining), the two comparisons, the uncertainty attribute |
| `10_specification_findings.py` | The places where the submitted description of the data, rather than the data, was wrong |
| `11_landslide_areas.py` | One row per landslide with its identifier and area; checks the stored `lds_area` against the polygons |
| `12_overlap_map.py` | The three-way partition (both / only this inventory / only the previous one) over the shared domain, as areas and as a map. **Needs non-public data, see below** |

Each script writes a plain-text report and one or more CSV tables to
`outputs/`. Lines marked `[!]` are discrepancies between the data and their
description.

## Requirements

Python 3.10 or later, with the GDAL/OGR Python bindings, shapely, pyproj and
matplotlib:

```
apt install python3-gdal             # or: pip install gdal==$(gdal-config --version)
pip install -r requirements.txt
```

Tested with Python 3.10.12, GDAL 3.8.4, shapely 2.0.6, pyproj 3.7.1 and
matplotlib 3.10.8. No network access and no proprietary software are required.

## Running

1. Download version 2 of the data record (https://doi.org/10.5281/zenodo.22814494) and
   place its layers in `data/`, so that `data/DRV_GLIM_O_I.shp` and the other
   layers sit directly in that folder. To keep them elsewhere, set
   `GLIM_DATA=/path/to/layers`.
2. Run

   ```
   ./run_all.sh
   ```

   or one script at a time, for example `python3 03_area_statistics.py`.

`outputs/` in this repository holds the reference results of the published
run. A run over the published data overwrites them, so `git diff outputs/`
shows at once whether anything differs.

### The data state digest

Every report starts with a **data state digest**, a short hash of the content
of the input layers. The reference outputs were produced with digest
**`2547949e5806`**. A run over the unmodified data record gives the same
digest; a different digest means the input data are not those of the paper,
and the numbers are not comparable. `00_data_fingerprint.py` lists the
SHA-256 hash of each file, to find which layer differs.

## Data not distributed here

### The inventory of Jagodnik et al. (2020)

Scripts `04_compare_jom.py` and `12_overlap_map.py`, and block C of
`09_validation_tables.py`, compare the G-LIM with the engineering-geological
inventory of Jagodnik et al. (2020). **Those data are not public.** They were
made available to the authors of the present work for this comparison only;
that permission does not extend to third parties, and the data are not
included in this repository. To reproduce these results, request the data from
the authors of Jagodnik et al. (2020), place the shapefiles in
`data/jagodnik2020/` (or set `GLIM_JOM`), and run the scripts again.

Without those data, scripts 04 and 12 write a report saying they were skipped,
block C of script 09 is skipped, and every other result is produced in full.
The reference outputs in `outputs/` include the aggregate results of that
comparison, as reported in the paper, but no feature-level table and no
geometry derived from the previous inventory.

### The dataset as originally submitted

The integrity table of the paper reports, for each check, how many records
were flagged in the dataset as originally submitted, how many were corrected
and how many remain. The "flagged" column comes from a run of scripts 02 and 07
over the dataset as first submitted, before the correction round. That version
of the dataset is not distributed. The findings of that run are, in
`reference/as_submitted/`, and script 09 reads them from there. The dataset as
originally submitted is available on request from the author of this
repository.

## Conventions

* All vector layers are EPSG:3765 (HTRS96 / Croatia TM). Areas quoted in the
  paper are ellipsoidal, on GRS 1980, matching the `Shape_Area` field of the
  data; planar grid areas, about 300 ppm larger, are reported alongside.
* `Shape_Area` is treated as **data to be checked**, never as the source of
  truth: script 02 compares it with the geometry.
* Geometries are **not** repaired in script 02, whose purpose is to find broken
  ones. Elsewhere `buffer(0)` is applied where a valid geometry is a
  precondition for an areal operation, and the reports say so.
* A *landslide* is a distinct `lds_id` **within one layer**; a *polygon* is one
  record. A landslide may be represented by more than one polygon, so the two
  counts differ by design.
* Findings are keyed by stable attribute identifiers (`lds_id/lds_class` for
  the landslide layers, `gm_id` for the geomorphological layer), never by the
  positional shapefile feature id. Every finding also carries the coordinates
  of its centroid.

## Transcribed data

`tosevski2013_inventory.csv` holds Tables 1 and 2 of Toševski (2013),
transcribed from the published paper, with the provenance and the published
control totals in its header. Script 05 checks the transcription against those
totals before using it. If you use these values, cite Toševski (2013).

## How to cite

Cite the paper above and the archived version of this code (the DOI is on the
Zenodo badge of the release; see also `CITATION.cff`).

## Licence

The code is released under the MIT licence (`LICENSE`). The inventory is
released separately under CC BY 4.0.

## References

* Galli, M., Ardizzone, F., Cardinali, M., Guzzetti, F. & Reichenbach, P.
  Comparing landslide inventory maps. *Geomorphology* **94**, 268–289 (2008).
* Jagodnik, P., Bernat Gazibara, S., Arbanas, Ž. & Mihalić Arbanas, S.
  Engineering geological mapping using airborne LiDAR datasets – an example
  from the Vinodol Valley, Croatia. *J. Maps* **16**, 855–866 (2020).
* Toševski, A. Inventar klizišta u slivu rijeke Dubračine [Landslide inventory
  of the Dubračina river basin]. *Rudarsko-geološko-naftni zbornik* **27**,
  1–16 (2013).

## Contact

Michele Santangelo, CNR-GEO (National Research Council of Italy, Institute of Geosciences) — michele.santangelo@cnr.it
