"""Step 4. Baton Rouge fault trace and historical withdrawals by sand, both taken
from the public USGS SEAWAT model archive for SIR 2019-5102 (Heywood, 2019,
doi:10.5066/P96YE03Q).

Fault: the model represents the Baton Rouge fault as horizontal-flow-barrier (HFB)
cell pairs. The shared edge of each cell pair in the "2,000-foot" sand layer
(model layer 24) is extracted and merged into a line. This is a model-grid
approximation (cells are 500 ft wide near Baton Rouge), not a surveyed fault trace.
The Denham Springs-Scotlandville fault to the north is not in this model and is not mapped here.

Withdrawals: the SEAWAT Well Package for scenario 1 lists pumping (ft3/d) for each
model cell and stress period. Stress periods 2 to 78 are calendar years 1940 to 2016.
Rates are summed by sand and converted to million gallons per day (Mgal/d).
"""
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.ops import linemerge, unary_union

from common import RAW, PROC, GISDATA, TAB, MODEL_LAYER_SAND, COUNTY_FIPS

SB = RAW / "sciencebase_P96YE03Q"
FT3D_TO_MGD = 7.48052 / 1e6


def unzip(name):
    d = SB / name.replace(".zip", "")
    if not d.exists():
        with zipfile.ZipFile(SB / name) as z:
            z.extractall(SB)
    return d


def parish():
    c = gpd.read_file(f"zip://{RAW / 'cb_2023_us_county_500k.zip'}")
    return c[c.GEOID == COUNTY_FIPS].to_crs(4326)


def read_hfb(path, group="hfb_2000ft"):
    rows, on = [], False
    for line in open(path, encoding="latin-1"):
        if line.startswith("hfb_"):
            on = line.split()[0] == group
            continue
        p = line.split()
        if on and len(p) == 6:
            rows.append([int(v) for v in p[:5]])
    return pd.DataFrame(rows, columns=["lay", "r1", "c1", "r2", "c2"])


def fault_line(grid):
    hfb = read_hfb(unzip("ancillary.zip") / "hfb.tpl")
    cell = grid.set_index(["ROW", "COL"]).geometry
    edges = []
    for _, h in hfb.iterrows():
        a, b = cell.get((h.r1, h.c1)), cell.get((h.r2, h.c2))
        if a is not None and b is not None:
            e = a.boundary.intersection(b.boundary)
            if e.length > 0:
                edges.append(e)
    line = linemerge(unary_union(edges))
    g = gpd.GeoDataFrame({"name": ["Baton Rouge fault (approximate, USGS model HFB cells, 2,000-foot sand layer)"],
                          "source": ["Heywood 2019, doi:10.5066/P96YE03Q, hfb.tpl group hfb_2000ft"],
                          "n_cell_edges": [len(edges)]}, geometry=[line], crs=grid.crs)
    return g.to_crs(4326), len(hfb)


def read_wel(path):
    lines = open(path, encoding="latin-1").read().replace("\r", "").split("\n")
    i, sp, recs, last = 1, 0, [], None
    while i < len(lines) and lines[i].strip():
        n = int(lines[i].split()[0])
        sp += 1
        i += 1
        if n < 0:
            block = last.copy()
            block["sp"] = sp
        else:
            arr = np.array([l.split()[:4] for l in lines[i:i + n]], dtype=float)
            block = pd.DataFrame(arr, columns=["lay", "row", "col", "q"])
            block["sp"] = sp
            i += n
        last = block
        recs.append(block)
    return pd.concat(recs, ignore_index=True)


def main():
    grid = gpd.read_file(unzip("georef.zip") / "sir20195102.shp")
    ebr = parish()

    # fault
    fault, n_hfb = fault_line(grid)
    fault.to_file(GISDATA / "baton_rouge_fault_approx.geojson", driver="GeoJSON")
    print(f"fault: {n_hfb} HFB cell pairs; length {fault.to_crs(26915).length.iloc[0] / 1000:.1f} km over model domain")

    # withdrawals
    unzip("model.zip")
    w = read_wel(SB / "model" / "SEAWAT-scen1" / "scenario1.wel")
    w = w[(w.sp >= 2) & (w.sp <= 78)].copy()
    w["year"] = w.sp + 1938
    w["withdrawal_mgd"] = np.where(w.q < 0, -w.q, 0.0) * FT3D_TO_MGD
    w["sand"] = w.lay.astype(int).map(MODEL_LAYER_SAND).fillna("confining-unit cells")
    cent = grid.copy()
    cent["geometry"] = grid.centroid
    inside = gpd.sjoin(cent.to_crs(4326), ebr[["geometry"]], predicate="within")
    ebr_cells = set(zip(inside.ROW, inside.COL))
    w["in_ebr"] = [(int(r), int(c)) in ebr_cells for r, c in zip(w.row, w.col)]

    out = []
    for area, sub in (("model domain", w), ("East Baton Rouge Parish", w[w.in_ebr])):
        t = sub.groupby(["year", "sand"]).withdrawal_mgd.sum().unstack(fill_value=0.0)
        t["total"] = t.sum(axis=1)
        t = t.reset_index().melt(id_vars="year", var_name="sand", value_name="withdrawal_mgd")
        t["area"] = area
        out.append(t)
    wd = pd.concat(out)
    wd.to_csv(PROC / "withdrawals_by_sand_1940_2016.csv", index=False, float_format="%.3f")
    tot = wd[wd.sand == "total"].pivot(index="year", columns="area", values="withdrawal_mgd")
    print(tot.loc[[1940, 1960, 1970, 1975, 1980, 1990, 2000, 2010, 2016]].round(1))
    peak = tot.idxmax()
    summary = []
    for area in tot.columns:
        s = wd[(wd.area == area) & (wd.year.isin([1975, 2016]))].pivot(index="sand", columns="year", values="withdrawal_mgd")
        s["peak_year"] = [wd[(wd.area == area) & (wd.sand == k)].set_index("year").withdrawal_mgd.idxmax() for k in s.index]
        s["peak_mgd"] = [wd[(wd.area == area) & (wd.sand == k)].withdrawal_mgd.max() for k in s.index]
        s["area"] = area
        summary.append(s.reset_index())
    summary = pd.concat(summary)
    summary.columns = [str(c) for c in summary.columns]
    summary.to_csv(TAB / "withdrawals_summary_by_sand.csv", index=False, float_format="%.2f")
    print(summary.round(2).to_string())
    print("peak total year:", peak.to_dict())


if __name__ == "__main__":
    main()
