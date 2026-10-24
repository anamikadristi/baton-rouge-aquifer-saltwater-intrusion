"""Step 2. Clean NWIS well metadata and groundwater levels for East Baton Rouge Parish.

Outputs (data/processed):
  wells.csv              one row per NWIS groundwater site: location, aquifer code, sand name, well depth
  gwlevels_annual.csv    annual mean water-level altitude (ft, NGVD29) per well and year
"""
import json

import numpy as np
import pandas as pd

from common import RAW, PROC, SANDS

BAD_QUAL = ("Pumping", "NoMeasurement", "ForeignSubstance", "Dry", "Obstruct", "Recently")


def read_rdb(path):
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    return df.iloc[1:].reset_index(drop=True)


def wells():
    s = read_rdb(RAW / "nwis_sites_ebr_gw.rdb")
    s = s[(s.state_cd == "22") & (s.county_cd == "033")]
    w = pd.DataFrame({
        "site_no": s.site_no,
        "station_nm": s.station_nm.str.strip(),
        "lat": pd.to_numeric(s.dec_lat_va),
        "lon": pd.to_numeric(s.dec_long_va),
        "alt_ft": pd.to_numeric(s.alt_va, errors="coerce"),
        "alt_datum": s.alt_datum_cd,
        "aqfr_cd": s.aqfr_cd.fillna("").str.strip(),
        "well_depth_ft": pd.to_numeric(s.well_depth_va, errors="coerce"),
        "hole_depth_ft": pd.to_numeric(s.hole_depth_va, errors="coerce"),
        "construction_dt": s.construction_dt,
    })
    w["sand"] = w.aqfr_cd.map(SANDS).fillna("other or unknown")
    return w.reset_index(drop=True)


def field_levels(site_set):
    d = pd.read_json(RAW / "gwlevels_62610.jsonl", lines=True, dtype={"value": str})
    d["site_no"] = d.monitoring_location_id.str.replace("USGS-", "", regex=False)
    d = d[d.site_no.isin(site_set)]
    n0 = len(d)
    q = d.qualifier.apply(lambda v: " ".join(v) if isinstance(v, list) else str(v))
    d = d[~q.str.contains("|".join(BAD_QUAL), na=False)]
    d["value"] = pd.to_numeric(d.value, errors="coerce")
    d = d.dropna(subset=["value"])
    d["date"] = pd.to_datetime(d.time, utc=True).dt.tz_localize(None)
    print(f"field measurements (62610) in parish: {n0}; after qualifier screen: {len(d)}")
    return d[["site_no", "date", "value"]].assign(source="field")


def daily_levels(w):
    p = RAW / "gwlevels_daily.jsonl"
    if not p.exists():
        print("no daily file; skipping")
        return pd.DataFrame(columns=["site_no", "date", "value", "source"])
    d = pd.read_json(p, lines=True, dtype={"value": str})
    d["site_no"] = d.monitoring_location_id.str.replace("USGS-", "", regex=False)
    d = d[d.site_no.isin(set(w.site_no))]
    d["value"] = pd.to_numeric(d.value, errors="coerce")
    d = d.dropna(subset=["value"])
    d["date"] = pd.to_datetime(d.time)
    alt = w.set_index("site_no").alt_ft
    # depth below land surface (72019) converted to altitude with the NWIS land-surface altitude
    dep = d.parameter_code.astype(str) == "72019"
    d.loc[dep, "value"] = d.loc[dep, "site_no"].map(alt) - d.loc[dep, "value"]
    d = d.dropna(subset=["value"])
    print(f"daily values: {len(d)} days at {d.site_no.nunique()} wells")
    return d[["site_no", "date", "value"]].assign(source="daily")


def main():
    w = wells()
    f = field_levels(set(w.site_no))
    dl = daily_levels(w)
    lev = pd.concat([f, dl], ignore_index=True)
    lev["year"] = lev.date.dt.year
    ann = (lev.groupby(["site_no", "year", "source"]).value.agg(["mean", "size"]).reset_index())
    # one value per well-year: daily mean if a recorder was running, else mean of field visits
    ann["pri"] = (ann.source == "daily").astype(int)
    ann = ann.sort_values("pri").groupby(["site_no", "year"]).tail(1).drop(columns="pri")
    ann = ann.rename(columns={"mean": "wl_alt_ft", "size": "n_obs"}).sort_values(["site_no", "year"])
    ann.to_csv(PROC / "gwlevels_annual.csv", index=False, float_format="%.2f")

    st = ann.groupby("site_no").agg(n_years=("year", "size"), first_year=("year", "min"), last_year=("year", "max"))
    w = w.merge(st, on="site_no", how="left")
    w[["n_years"]] = w[["n_years"]].fillna(0).astype(int)
    w.to_csv(PROC / "wells.csv", index=False)
    print(f"wells in parish: {len(w)}; with any water level: {(w.n_years > 0).sum()}; "
          f"well-years: {len(ann)}; years {ann.year.min()} to {ann.year.max()}")
    print(w[w.n_years > 0].groupby("sand").agg(wells=("site_no", "size"), wells_20yr=("n_years", lambda x: (x >= 20).sum())).sort_values("wells", ascending=False))


if __name__ == "__main__":
    main()
