"""Step 3. Clean Water Quality Portal chloride and specific conductance for wells in
East Baton Rouge Parish and attach aquifer, depth and position relative to the
Baton Rouge fault (fault line comes from step 4, so run step 4 first).

Rules
- Routine samples only (field replicates dropped).
- Chloride: USGS parameter 00940 (dissolved, mg/L). WQX 3.0 rows (2023 onward) are
  added where the legacy profile has no result for the same well and date.
- Specific conductance: field value (00095) preferred over laboratory value (90095)
  for the same sample.
- Several results for the same well and date are averaged.
- Aquifer and well depth come from NWIS site metadata (step 2); for the few non-USGS
  wells, WQP station well depth is used and the aquifer is left unknown.
- Distance to the fault is signed: positive north of the fault, negative south.

Output: data/processed/chloride_spcond_samples.csv
"""
import geopandas as gpd
import numpy as np
import pandas as pd

from common import RAW, PROC, GISDATA, UTM15


def legacy(name, char):
    d = pd.read_csv(RAW / f"wqp_{name}.csv", dtype=str, low_memory=False)
    d = d[d.ActivityTypeCode.isin(["Sample-Routine", "Sample", "Field Msr/Obs"])]
    out = pd.DataFrame({
        "wqp_id": d.MonitoringLocationIdentifier,
        "date": pd.to_datetime(d.ActivityStartDate),
        "value": pd.to_numeric(d.ResultMeasureValue, errors="coerce"),
        "pcode": d.USGSPCode.fillna(""),
        "censored": d.ResultDetectionConditionText.notna(),
    })
    return out


def wqx3(name):
    d = pd.read_csv(RAW / f"wqp3_{name}_2023on.csv", dtype=str, low_memory=False)
    d = d[d.Activity_TypeCode.isin(["Sample-Routine", "Sample", "Field Msr/Obs"])]
    lab = d.Result_CharacteristicUserSupplied.fillna("").str.contains("laboratory")
    return pd.DataFrame({
        "wqp_id": d.Location_Identifier,
        "date": pd.to_datetime(d.Activity_StartDate),
        "value": pd.to_numeric(d.Result_Measure, errors="coerce"),
        "pcode": np.where(name == "spcond", np.where(lab, "90095", "00095"), "00940"),
        "censored": d.Result_ResultDetectionCondition.notna() if "Result_ResultDetectionCondition" in d else False,
    })


def combine(name, char):
    a, b = legacy(name, char), wqx3(name)
    key = set(zip(a.wqp_id, a.date))
    b = b[[k not in key for k in zip(b.wqp_id, b.date)]]
    d = pd.concat([a, b], ignore_index=True).dropna(subset=["value"])
    if name == "spcond":
        d["pri"] = (d.pcode != "00095").astype(int)
        d = d.sort_values("pri")
        best = d.groupby(["wqp_id", "date"]).pri.transform("min")
        d = d[d.pri == best]
    g = d.groupby(["wqp_id", "date"]).agg(value=("value", "mean"), censored=("censored", "max")).reset_index()
    return g, len(a), len(b)


def fault_side(df):
    f = gpd.read_file(GISDATA / "baton_rouge_fault_approx.geojson").to_crs(UTM15)
    line = f.geometry.iloc[0]
    pts = gpd.GeoSeries(gpd.points_from_xy(df.lon, df.lat), crs=4326).to_crs(UTM15)
    dist = pts.distance(line).values
    # fault northing at each well easting (fault strikes roughly east-west)
    coords = np.vstack([np.asarray(g.coords) for g in getattr(line, "geoms", [line])])
    o = np.argsort(coords[:, 0])
    fy = np.interp(pts.x.values, coords[o, 0], coords[o, 1])
    north = pts.y.values > fy
    df["x_utm"], df["y_utm"] = pts.x.values, pts.y.values
    df["dist_fault_km"] = np.where(north, 1, -1) * dist / 1000.0
    df["side"] = np.where(north, "north", "south")
    return df


def main():
    wells = pd.read_csv(PROC / "wells.csv", dtype={"site_no": str})
    st = pd.read_csv(RAW / "wqp_stations.csv", dtype=str)
    cl, na, nb = combine("chloride", "Chloride")
    sc, sa, sb = combine("spcond", "Specific conductance")
    print(f"chloride rows: legacy {na}, WQX3 extra {nb}; well-dates {len(cl)}")
    print(f"spcond rows: legacy {sa}, WQX3 extra {sb}; well-dates {len(sc)}")
    d = cl.rename(columns={"value": "chloride_mgL", "censored": "cl_censored"}).merge(
        sc.rename(columns={"value": "spcond_uScm", "censored": "sc_censored"}), on=["wqp_id", "date"], how="outer")
    d["site_no"] = d.wqp_id.str.replace("USGS-", "", regex=False)
    meta = wells.set_index("site_no")
    usgs = d.wqp_id.str.startswith("USGS-") & d.site_no.isin(meta.index)
    for col in ("lat", "lon", "aqfr_cd", "sand", "well_depth_ft", "station_nm"):
        d[col] = np.where(usgs, d.site_no.map(meta[col]), np.nan)
    s = st.set_index("MonitoringLocationIdentifier")
    other = ~usgs
    d.loc[other, "lat"] = pd.to_numeric(d.loc[other, "wqp_id"].map(s.LatitudeMeasure), errors="coerce")
    d.loc[other, "lon"] = pd.to_numeric(d.loc[other, "wqp_id"].map(s.LongitudeMeasure), errors="coerce")
    d.loc[other, "well_depth_ft"] = pd.to_numeric(d.loc[other, "wqp_id"].map(s["WellDepthMeasure/MeasureValue"]), errors="coerce")
    d.loc[other, "sand"] = "other or unknown"
    d["lat"], d["lon"] = d.lat.astype(float), d.lon.astype(float)
    d["well_depth_ft"] = pd.to_numeric(d.well_depth_ft, errors="coerce")
    d = d.dropna(subset=["lat", "lon"])
    d = fault_side(d)
    d["year"] = d.date.dt.year
    d.sort_values(["wqp_id", "date"]).to_csv(PROC / "chloride_spcond_samples.csv", index=False, float_format="%.5g")
    c = d.dropna(subset=["chloride_mgL"])
    print(f"chloride: {len(c)} well-dates at {c.wqp_id.nunique()} wells, {c.year.min()} to {c.year.max()}")
    print(f"spcond: {d.spcond_uScm.notna().sum()} well-dates at {d[d.spcond_uScm.notna()].wqp_id.nunique()} wells")
    print(c.groupby("side").wqp_id.nunique())
    print(c[c.year >= 2000].groupby("sand").wqp_id.nunique().sort_values(ascending=False))


if __name__ == "__main__":
    main()
