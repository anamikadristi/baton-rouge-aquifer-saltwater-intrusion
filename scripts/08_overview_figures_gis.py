"""Step 8. Overview map, withdrawals figure, and GIS exports for ArcGIS Pro, QGIS and Google Earth.

GIS outputs (outputs/gis):
  wells_water_level.shp        NWIS wells with water levels, aquifer code, sand, depth, record length
  baton_rouge_fault_approx.shp fault line from the USGS model grid (approximate)
  baton_rouge_fault_approx.kml same line for Google Earth
  chloride_latest_since2000.shp and kriged GeoTIFFs are written by step 6
"""
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import RAW, PROC, GISDATA, GISOUT, FIG, TAB, COUNTY_FIPS, MAIN_SANDS, SANDS, UTM15


def main():
    wells = pd.read_csv(PROC / "wells.csv", dtype={"site_no": str})
    cl = pd.read_csv(PROC / "chloride_spcond_samples.csv", dtype={"aqfr_cd": str})
    fault = gpd.read_file(GISDATA / "baton_rouge_fault_approx.geojson")
    c = gpd.read_file(f"zip://{RAW / 'cb_2023_us_county_500k.zip'}")
    ebr = c[c.GEOID == COUNTY_FIPS].to_crs(4326)
    near = c[c.STATEFP == "22"].to_crs(4326)

    gw = wells[wells.n_years > 0]
    g = gpd.GeoDataFrame(gw, geometry=gpd.points_from_xy(gw.lon, gw.lat), crs=4326)
    g.rename(columns={"station_nm": "station", "well_depth_ft": "depth_ft", "first_year": "first_yr",
                      "last_year": "last_yr", "construction_dt": "constr_dt", "hole_depth_ft": "hole_ft"})[
        ["site_no", "station", "aqfr_cd", "sand", "depth_ft", "alt_ft", "n_years", "first_yr", "last_yr", "geometry"]
    ].to_file(GISOUT / "wells_water_level.shp")
    fault.rename(columns={"n_cell_edges": "n_edges"}).to_file(GISOUT / "baton_rouge_fault_approx.shp")
    try:
        fault[["name", "geometry"]].rename(columns={"name": "Name"}).to_file(GISOUT / "baton_rouge_fault_approx.kml", driver="KML")
    except Exception as e:
        print("KML export skipped:", e)

    # fig01 overview map
    fig, ax = plt.subplots(figsize=(7.5, 7.0))
    near.plot(ax=ax, color="#f2f2f2", edgecolor="0.75", lw=0.5)
    ebr.plot(ax=ax, color="white", edgecolor="0.2", lw=1.0)
    depth = g.well_depth_ft.clip(upper=3000)
    sc = ax.scatter(g.lon, g.lat, c=depth, cmap="viridis_r", s=6 + 2 * np.sqrt(g.n_years), alpha=0.8,
                    edgecolor="none", vmin=0, vmax=3000)
    clw = cl.dropna(subset=["chloride_mgL"]).drop_duplicates("wqp_id")
    ax.scatter(clw.lon, clw.lat, marker="+", s=12, color="#c0392b", lw=0.6, label=f"Chloride wells (WQP), n = {len(clw)}")
    fault.plot(ax=ax, color="k", lw=1.8, ls="--", label="Baton Rouge fault (approximate, USGS model HFB cells)")
    xmin, ymin, xmax, ymax = ebr.total_bounds
    ax.set_xlim(xmin - 0.03, xmax + 0.03)
    ax.set_ylim(ymin - 0.03, ymax + 0.03)
    ax.set_aspect(1 / np.cos(np.radians(30.5)))
    cb = fig.colorbar(sc, ax=ax, shrink=0.7)
    cb.set_label("Well depth (ft); marker size grows with years of water-level record")
    ax.legend(loc="lower left", fontsize=7.5, frameon=True)
    ax.set_title(f"East Baton Rouge Parish: {len(g)} NWIS wells with water levels, {len(clw)} wells with chloride", fontsize=9.5)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    fig.tight_layout()
    fig.savefig(FIG / "fig01_overview_map.png", dpi=150)
    plt.close(fig)

    # fig03 withdrawals and water levels
    wd = pd.read_csv(PROC / "withdrawals_by_sand_1940_2016.csv")
    e = wd[(wd.area == "East Baton Rouge Parish") & (wd.sand != "total")].pivot(index="year", columns="sand", values="withdrawal_mgd")
    order = [SANDS[k] for k in MAIN_SANDS if SANDS[k] in e.columns]
    order = [s for s in order if e[s].max() > 0.5]
    comp = pd.read_csv(PROC / "composite_hydrograph_by_sand.csv")
    seg = pd.read_csv(TAB / "piecewise_regression_by_sand.csv")
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 9), sharex=True)
    axes[0].stackplot(e.index, [e[s] for s in order], labels=order, colors=plt.get_cmap("tab10").colors[:len(order)])
    axes[0].set_ylabel("Withdrawals (Mgal/d)")
    axes[0].set_title("Groundwater withdrawals in East Baton Rouge Parish by sand, SEAWAT model input (Heywood, 2019)", fontsize=9)
    axes[0].legend(fontsize=7, ncol=2, loc="upper left", frameon=False)
    for ax, s in zip(axes[1:], ["2,000-foot sand", "1,500-foot sand"]):
        cc = comp[comp.sand == s]
        ax.plot(cc.year, cc.index_ft, ".", color="0.3", ms=3, label="Composite water-level index (NWIS)")
        r = seg[seg.sand == s].iloc[0]
        ax.axvline(r.bp1_year, color="#c0392b", ls="--", lw=1, label=f"Single breakpoint {int(r.bp1_year)} (90% CI {int(r.bp1_ci90_lo)} to {int(r.bp1_ci90_hi)})")
        ax.set_ylabel("Index (ft)")
        a2 = ax.twinx()
        a2.plot(e.index, e[s], color="#1f5fa6", lw=1.4, label="Withdrawals")
        a2.set_ylabel("Withdrawals (Mgal/d)", color="#1f5fa6")
        ax.set_title(s, fontsize=9)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = a2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="lower left", frameon=False)
    axes[-1].set_xlim(1935, 2027)
    fig.tight_layout()
    fig.savefig(FIG / "fig03_withdrawals_and_levels.png", dpi=150)
    plt.close(fig)
    print("figures and GIS exports written")


if __name__ == "__main__":
    main()
