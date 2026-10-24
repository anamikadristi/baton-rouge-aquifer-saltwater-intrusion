"""Shared paths and constants for the Baton Rouge aquifer project."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Raw downloads go to data/raw by default. Set BRAQ_RAW_DIR to keep them elsewhere.
RAW = Path(os.environ.get("BRAQ_RAW_DIR", ROOT / "data" / "raw"))
PROC = ROOT / "data" / "processed"
GISDATA = ROOT / "data" / "gis"
FIG = ROOT / "outputs" / "figures"
TAB = ROOT / "outputs" / "tables"
GISOUT = ROOT / "outputs" / "gis"
for p in (RAW, PROC, GISDATA, FIG, TAB, GISOUT):
    p.mkdir(parents=True, exist_ok=True)

COUNTY_FIPS = "22033"  # East Baton Rouge Parish, Louisiana
# Search box slightly larger than the parish; wells are filtered to the parish afterwards.
BBOX = (-91.36, 30.27, -90.84, 30.74)

# NWIS local aquifer codes for the Baton Rouge area (names from the USGS
# Water Data API aquifer-codes collection). Order is shallow to deep.
SANDS = {
    "112SLBR": "Shallow sands",
    "112MRVA": "Mississippi River alluvial",
    "11204BR": "400-foot sand",
    "11205BR": "400- and 600-foot sands",
    "11206BR": "600-foot sand",
    "11207BR": "600- and 800-foot sands",
    "12108BR": "800-foot sand",
    "12109BR": "800- and 1,200-foot sands",
    "12110BR": "1,000-foot sand",
    "12111BR": "1,000- and 1,200-foot sands",
    "12112BR": "1,200-foot sand",
    "12113BR": "1,200- and 1,500-foot sands",
    "12115BR": "1,500-foot sand",
    "12116BR": "1,500- and 1,700-foot sands",
    "12117BR": "1,700-foot sand",
    "12220BR": "2,000-foot sand",
    "12223BR": "2,000- and 2,400-foot sands",
    "12224BR": "2,400-foot sand",
    "12226BR": "2,400- and 2,800-foot sands",
    "12228BR": "2,800-foot sand",
}
# Single named sands used for trend analysis, shallow to deep.
MAIN_SANDS = ["11204BR", "11206BR", "12108BR", "12110BR", "12112BR", "12115BR",
              "12117BR", "12220BR", "12224BR", "12228BR"]

# Layers of the USGS SEAWAT model (Heywood and others, 2019, SIR 2019-5102)
# that represent each sand. Even layers 2 to 10 and 21 to 27 are confining units.
MODEL_LAYER_SAND = {1: "400-foot sand", 3: "600-foot sand", 5: "800-foot sand",
                    7: "1,000-foot sand", 9: "1,200-foot sand",
                    **{k: "1,500-foot sand" for k in range(11, 21)},
                    22: "1,700-foot sand", 24: "2,000-foot sand",
                    26: "2,400-foot sand", 28: "2,800-foot sand"}

UTM15 = "EPSG:26915"  # NAD83 / UTM zone 15N, metres
