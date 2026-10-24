"""Step 1. Download all public data used in this project (keyless endpoints only).

Sources
1. USGS NWIS site service: groundwater sites in East Baton Rouge Parish with
   aquifer code and well depth (expanded site output, RDB).
2. USGS Water Data OGC API (api.waterdata.usgs.gov), collection
   "field-measurements": discrete groundwater levels. The legacy NWIS gwlevels
   service was retired in fall 2025, so the new API is used. Parameter codes:
   62610 (water-level altitude, NGVD29), 62611 (NAVD88), 72019 (depth below land surface).
   Daily means from continuous recorders come from the "daily" collection.
3. Water Quality Portal (WQX 3.0 profile): chloride and specific conductance in
   wells in East Baton Rouge Parish, plus station metadata.
4. USGS ScienceBase data release P96YE03Q (Heywood, 2019), the SEAWAT model archive
   for SIR 2019-5102: model grid shapefile (georef.zip), calibration templates with the
   Baton Rouge fault horizontal-flow-barrier cells (ancillary.zip), and model input
   (model.zip, which holds the 1940 to 2016 well withdrawals by model layer).
5. U.S. Census Bureau 2023 cartographic boundary file of counties (parish outline).

Run: python scripts/01_download_data.py   (about 10 to 20 minutes; files are skipped if present)
"""
import json
import time
import urllib.parse
from pathlib import Path

import requests

from common import RAW, BBOX, COUNTY_FIPS

S = requests.Session()
S.headers["User-Agent"] = "public-data-research-script (baton-rouge-aquifer project)"


def get(url, params=None, out=None, timeout=600, tries=4):
    for k in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if r.status_code == 429:
                wait = 60 * (k + 1)
                print(f"  rate limited, waiting {wait} s")
                time.sleep(wait)
                continue
            r.raise_for_status()
            if out:
                Path(out).write_bytes(r.content)
            return r
        except requests.RequestException as e:
            print(f"  attempt {k + 1} failed: {e}")
            time.sleep(10 * (k + 1))
    raise RuntimeError(f"failed: {url}")


def nwis_sites():
    out = RAW / "nwis_sites_ebr_gw.rdb"
    if out.exists():
        return
    print("NWIS site metadata")
    get("https://waterservices.usgs.gov/nwis/site/",
        params={"format": "rdb", "countyCd": COUNTY_FIPS, "siteType": "GW",
                "siteOutput": "expanded"}, out=out)


def field_measurements():
    base = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements/items"
    for pc in ("62610", "62611", "72019"):
        out = RAW / f"gwlevels_{pc}.jsonl"
        if out.exists():
            continue
        print(f"Water Data API field measurements, parameter {pc}")
        params = {"bbox": ",".join(map(str, BBOX)), "parameter_code": pc,
                  "f": "json", "limit": 10000, "skipGeometry": "true",
                  "properties": "monitoring_location_id,parameter_code,value,unit_of_measure,"
                                "time,qualifier,vertical_datum,approval_status,observing_procedure"}
        url = base + "?" + urllib.parse.urlencode(params)
        n = 0
        tmp = out.with_suffix(".part")
        with open(tmp, "w") as fh:
            while url:
                d = get(url).json()
                for f in d.get("features", []):
                    fh.write(json.dumps(f["properties"]) + "\n")
                    n += 1
                url = next((l["href"] for l in d.get("links", []) if l.get("rel") == "next"), None)
                print(f"  {n} records")
        tmp.rename(out)


def daily_levels():
    """Daily mean groundwater levels from continuous recorders (statistic 00003)."""
    out = RAW / "gwlevels_daily.jsonl"
    if out.exists():
        return
    print("Water Data API daily groundwater-level series")
    meta = get("https://api.waterdata.usgs.gov/ogcapi/v0/collections/time-series-metadata/items",
               params={"bbox": ",".join(map(str, BBOX)), "f": "json", "limit": 5000,
                       "skipGeometry": "true"}).json()
    series = [f["properties"] for f in meta["features"]
              if f["properties"]["parameter_code"] in ("72019", "62610", "62611")
              and f["properties"]["computation_period_identifier"] == "Daily"
              and f["properties"]["statistic_id"] == "00003"]
    (RAW / "daily_series_metadata.json").write_text(json.dumps(series, indent=1))
    base = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items"
    tmp = out.with_suffix(".part")
    with open(tmp, "w") as fh:
        for s in series:
            params = {"monitoring_location_id": s["monitoring_location_id"],
                      "parameter_code": s["parameter_code"], "statistic_id": "00003",
                      "f": "json", "limit": 10000, "skipGeometry": "true",
                      "properties": "monitoring_location_id,parameter_code,statistic_id,time,value,"
                                    "unit_of_measure,approval_status,qualifier"}
            url = base + "?" + urllib.parse.urlencode(params)
            n = 0
            while url:
                d = get(url).json()
                for f in d.get("features", []):
                    fh.write(json.dumps(f["properties"]) + "\n")
                    n += 1
                url = next((l["href"] for l in d.get("links", []) if l.get("rel") == "next"), None)
            print(f"  {s['monitoring_location_id']} {s['parameter_code']}: {n} days")
    tmp.rename(out)


def wqp_get(url, params, out):
    """WQP writes an INCOMPLETE DATA line at the end when a request times out; retry then."""
    for k in range(4):
        get(url, params=params, out=out, timeout=1800)
        tail = Path(out).read_bytes()[-400:].decode("latin-1")
        if "INCOMPLETE DATA" not in tail:
            return
        print(f"  incomplete response, retry {k + 1}")
        time.sleep(30)
    raise RuntimeError(f"WQP response still incomplete: {out}")


def wqp():
    """Legacy WQP profile for the full record (for this query its USGS results end in March 2023),
    plus the WQX 3.0 profile from 2023 onward for newer USGS results. Duplicates are
    removed in step 3."""
    common = {"countycode": "US:22:033", "siteType": "Well", "mimeType": "csv", "zip": "no"}
    for name, char in (("chloride", "Chloride"), ("spcond", "Specific conductance")):
        out = RAW / f"wqp_{name}.csv"
        if not out.exists():
            print(f"WQP results (legacy profile): {char}")
            wqp_get("https://www.waterqualitydata.us/data/Result/search",
                    {**common, "characteristicName": char}, out)
        out = RAW / f"wqp3_{name}_2023on.csv"
        if not out.exists():
            print(f"WQP results (WQX 3.0, 2023 on): {char}")
            wqp_get("https://www.waterqualitydata.us/wqx3/Result/search",
                    {**common, "characteristicName": char, "dataProfile": "basicPhysChem",
                     "startDateLo": "01-01-2023"}, out)
    out = RAW / "wqp_stations.csv"
    if not out.exists():
        print("WQP stations")
        wqp_get("https://www.waterqualitydata.us/data/Station/search",
                {**common, "characteristicName": ["Chloride", "Specific conductance"]}, out)


def sciencebase():
    item = "631405c6d34e36012efa31db"  # DOI 10.5066/P96YE03Q
    d = get(f"https://www.sciencebase.gov/catalog/item/{item}", params={"format": "json"}).json()
    sb = RAW / "sciencebase_P96YE03Q"
    sb.mkdir(exist_ok=True)
    (sb / "item.json").write_text(json.dumps(d, indent=1))
    for f in d["files"]:
        if f["name"] in ("georef.zip", "ancillary.zip", "model.zip", "readme.txt", "modelgeoref.txt"):
            out = sb / f["name"]
            if not out.exists():
                print(f"ScienceBase {f['name']} ({f['size']} bytes)")
                get(f["url"], out=out)


def census():
    out = RAW / "cb_2023_us_county_500k.zip"
    if not out.exists():
        print("Census county boundaries")
        get("https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_county_500k.zip", out=out)


if __name__ == "__main__":
    print(f"Raw data folder: {RAW}")
    nwis_sites()
    field_measurements()
    daily_levels()
    wqp()
    sciencebase()
    census()
    (RAW / "ACCESS_DATE.txt").write_text(time.strftime("%Y-%m-%d") + "\n")
    print("done")
