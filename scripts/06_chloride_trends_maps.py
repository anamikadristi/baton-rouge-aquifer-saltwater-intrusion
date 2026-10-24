"""Step 6. Chloride trends at wells near the Baton Rouge fault, and chloride maps by
depth interval with interpolation north of the fault.

Trends: wells with at least 8 sample dates over at least 10 years. Mann-Kendall test
(original test; samples are irregular, so no serial-correlation correction) and Sen
slope in mg/L per year, computed on actual decimal dates. "Near fault" means within
8 km of the approximate fault line on either side.

Maps: most recent chloride sample since 2000 at each well, in three depth intervals
defined by NWIS aquifer code:
  shallow       "400-foot" to "800-foot" sands
  intermediate  "1,000-foot" to "1,700-foot" sands
  deep          "2,000-foot" to "2,800-foot" sands
Ordinary kriging of log10(chloride) is run only for wells north of the fault and only
inside the parish north of the fault (the fault is a hydraulic barrier, so values are
not interpolated across it). Leave-one-out cross-validation (LOOCV) is reported, with
skill = 1 - MSE(kriging) / MSE(mean). An interval is kriged only if it has at least 20
wells north of the fault, and the surface is shown and written to GeoTIFF only if the
LOOCV skill is above zero.
Kriging pools several sands in each interval, so the surfaces are screening maps, not
plume delineations.
"""
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import pymannkendall as mk
import rasterio
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from pykrige.ok import OrdinaryKriging
from shapely.geometry import box
from shapely.ops import split

from common import RAW, PROC, GISDATA, GISOUT, FIG, TAB, UTM15, COUNTY_FIPS

INTERVALS = {
    "shallow (400- to 800-foot sands)": ["11204BR", "11205BR", "11206BR", "11207BR", "12108BR"],
    "intermediate (1,000- to 1,700-foot sands)": ["12109BR", "12110BR", "12111BR", "12112BR", "12113BR",
                                                   "12115BR", "12116BR", "12117BR"],
    "deep (2,000- to 2,800-foot sands)": ["12220BR", "12223BR", "12224BR", "12226BR", "12228BR"],
}
SHORT = {k: k.split()[0] for k in INTERVALS}
NEAR_KM = 8.0
CELL = 250.0  # m


def decimal_year(d):
    d = pd.to_datetime(d)
    return d.dt.year + (d.dt.dayofyear - 0.5) / 365.25


def trends(c):
    rows = []
    for wid, g in c.groupby("wqp_id"):
        g = g.sort_values("date")
        t = decimal_year(g.date).values
        if len(g) < 8 or t.max() - t.min() < 10:
            continue
        y = g.chloride_mgL.values
        r = mk.original_test(y)
        i, j = np.triu_indices(len(t), 1)
        ok = t[j] > t[i]
        sen = np.median((y[j] - y[i])[ok] / (t[j] - t[i])[ok])
        rows.append(dict(wqp_id=wid, station_nm=g.station_nm.iloc[0], sand=g.sand.iloc[0],
                         well_depth_ft=g.well_depth_ft.iloc[0], side=g.side.iloc[0],
                         dist_fault_km=g.dist_fault_km.iloc[0], n=len(g),
                         first=int(g.year.min()), last=int(g.year.max()),
                         first_cl=y[0], last_cl=y[-1], max_cl=y.max(),
                         trend=r.trend, p=r.p, tau=r.Tau, sen_mgL_per_yr=sen))
    return pd.DataFrame(rows)


def north_domain():
    c = gpd.read_file(f"zip://{RAW / 'cb_2023_us_county_500k.zip'}")
    ebr = c[c.GEOID == COUNTY_FIPS].to_crs(UTM15)
    f = gpd.read_file(GISDATA / "baton_rouge_fault_approx.geojson").to_crs(UTM15)
    poly = ebr.geometry.iloc[0]
    parts = split(poly, f.geometry.iloc[0])
    north = max(parts.geoms, key=lambda p: p.centroid.y)
    return ebr, f, north


def krige(sub, north, name):
    x, y = sub.x_utm.values / 1000.0, sub.y_utm.values / 1000.0
    z = np.log10(sub.chloride_mgL.values)
    best = None
    for model in ("exponential", "spherical"):
        pred = np.empty_like(z)
        for k in range(len(z)):
            m = np.ones(len(z), bool)
            m[k] = False
            ok = OrdinaryKriging(x[m], y[m], z[m], variogram_model=model, nlags=8)
            pred[k] = ok.execute("points", x[k:k + 1], y[k:k + 1])[0][0]
        rmse = float(np.sqrt(np.mean((pred - z) ** 2)))
        if best is None or rmse < best[1]:
            best = (model, rmse, pred)
    model, rmse, pred = best
    null_rmse = float(np.sqrt(np.mean((z - z.mean()) ** 2)))
    ok = OrdinaryKriging(x, y, z, variogram_model=model, nlags=8)
    xmin, ymin, xmax, ymax = north.bounds
    gx = np.arange(xmin + CELL / 2, xmax, CELL)
    gy = np.arange(ymax - CELL / 2, ymin, -CELL)
    zz, ss = ok.execute("grid", gx / 1000.0, gy / 1000.0)
    tr = from_origin(xmin, ymax, CELL, CELL)
    mask = geometry_mask([north], out_shape=zz.shape, transform=tr)
    zz = np.where(mask, np.nan, 10 ** np.asarray(zz))
    ss = np.where(mask, np.nan, np.asarray(ss))
    skill = 1 - rmse ** 2 / null_rmse ** 2
    # rasters are written only when kriging beats the plain mean in cross-validation
    for arr, tag in (((zz, "chloride_mgL"), (ss, "log10_kriging_variance")) if skill > 0 else ()):
        with rasterio.open(GISOUT / f"kriged_{tag}_{name}_north_of_fault.tif", "w", driver="GTiff",
                           height=arr.shape[0], width=arr.shape[1], count=1, dtype="float32",
                           crs=UTM15, transform=tr, nodata=-9999, compress="deflate") as dst:
            dst.write(np.where(np.isnan(arr), -9999, arr).astype("float32"), 1)
    stats = dict(interval=name, n_north=len(z), variogram=model, loocv_rmse_log10=rmse,
                 null_rmse_log10=null_rmse, loocv_skill=skill,
                 variogram_params=";".join(f"{p:.3g}" for p in ok.variogram_model_parameters))
    return zz, (xmin, xmax, ymin, ymax), stats


def main():
    d = pd.read_csv(PROC / "chloride_spcond_samples.csv", dtype={"aqfr_cd": str, "site_no": str}, parse_dates=["date"])
    c = d.dropna(subset=["chloride_mgL"]).copy()
    c = c[c.chloride_mgL > 0]

    # ---- trends
    t = trends(c)
    t["near_fault"] = t.dist_fault_km.abs() <= NEAR_KM
    t.sort_values("dist_fault_km").to_csv(TAB / "chloride_trends_by_well.csv", index=False, float_format="%.4g")
    s = (t.groupby(["near_fault", "side"]).agg(wells=("wqp_id", "size"),
                                               increasing=("trend", lambda v: (v == "increasing").sum()),
                                               decreasing=("trend", lambda v: (v == "decreasing").sum()),
                                               median_sen=("sen_mgL_per_yr", "median")).reset_index())
    s.to_csv(TAB / "chloride_trend_summary.csv", index=False, float_format="%.3f")
    print(s.to_string())
    near = t[t.near_fault & (t.side == "north")].sort_values("sen_mgL_per_yr", ascending=False)
    print(near[["station_nm", "sand", "dist_fault_km", "n", "first", "last", "first_cl", "last_cl",
                "trend", "p", "sen_mgL_per_yr"]].round(3).to_string())

    # figure: chloride time series, wells north of fault with significant increase
    inc = near[near.trend == "increasing"].head(8)
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    cmap = plt.get_cmap("tab10")
    for k, (_, r) in enumerate(inc.iterrows()):
        g = c[c.wqp_id == r.wqp_id].sort_values("date")
        ax.plot(g.date, g.chloride_mgL, "-o", ms=3, lw=1, color=cmap(k),
                label=f"{r.station_nm}, {r.sand}, {r.dist_fault_km:.1f} km N, {r.sen_mgL_per_yr:+.1f} mg/L/yr")
    ax.set_yscale("log")
    ax.axhline(250, color="0.4", ls=":", lw=1)
    ax.text(0.99, 250, "250 mg/L secondary drinking-water standard ", fontsize=7, color="0.3", va="bottom", ha="right",
            transform=ax.get_yaxis_transform())
    ax.set_ylabel("Chloride (mg/L, log scale)")
    ax.set_title(f"Wells within {NEAR_KM:.0f} km north of the Baton Rouge fault with a significant chloride increase", fontsize=9.5)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig04_chloride_trends_near_fault.png", dpi=150)
    plt.close(fig)

    # ---- maps by depth interval
    ebr, fault, north = north_domain()
    recent = c[c.year >= 2000].sort_values("date").groupby("wqp_id").tail(1)
    code2int = {code: k for k, codes in INTERVALS.items() for code in codes}
    recent["interval"] = recent.aqfr_cd.map(code2int)
    recent = recent.dropna(subset=["interval"])
    recent.to_csv(PROC / "chloride_latest_since2000_by_well.csv", index=False, float_format="%.5g")
    gpd.GeoDataFrame(recent[["wqp_id", "station_nm", "sand", "interval", "well_depth_ft", "date", "chloride_mgL",
                             "side", "dist_fault_km"]].assign(date=lambda x: x.date.dt.strftime("%Y-%m-%d")),
                     geometry=gpd.points_from_xy(recent.lon, recent.lat), crs=4326).to_file(
        GISOUT / "chloride_latest_since2000.shp")

    stats, fig, axes = [], *plt.subplots(1, 3, figsize=(13.5, 5.0))
    norm = LogNorm(vmin=1, vmax=5000)
    for ax, (name, codes) in zip(axes, INTERVALS.items()):
        sub = recent[recent.interval == name]
        nsub = sub[sub.side == "north"]
        ebr.boundary.plot(ax=ax, color="0.3", lw=0.8)
        if len(nsub) >= 20:
            zz, ext, st = krige(nsub, north, SHORT[name])
            stats.append(st)
            if st["loocv_skill"] > 0:
                ax.imshow(zz, extent=ext, origin="upper", cmap="YlOrRd", norm=norm, alpha=0.75)
                sk = f"kriged, LOOCV skill {st['loocv_skill']:.2f}"
            else:
                sk = f"not mapped: LOOCV skill {st['loocv_skill']:.2f}"
        else:
            sk = "too few wells to interpolate"
        fault.plot(ax=ax, color="k", lw=1.6, ls="--")
        sc = ax.scatter(sub.x_utm, sub.y_utm, c=sub.chloride_mgL, cmap="YlOrRd", norm=norm, s=22,
                        edgecolor="k", linewidth=0.4, zorder=3)
        xmin, ymin, xmax, ymax = ebr.total_bounds
        ax.set_xlim(xmin - 2000, xmax + 2000)
        ax.set_ylim(ymin - 2000, ymax + 2000)
        ax.set_title(f"{name}\n{len(sub)} wells ({len(nsub)} north of fault)\n{sk}", fontsize=8.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_aspect("equal")
    cb = fig.colorbar(sc, ax=axes, shrink=0.75, pad=0.01)
    cb.set_label("Chloride (mg/L), latest sample since 2000")
    fig.suptitle("Chloride by depth interval, East Baton Rouge Parish (WQP). Dashed line: Baton Rouge fault, approximate (USGS model grid).\n"
                 "Shading: ordinary kriging of log10 chloride north of the fault only, shown where cross-validation skill > 0.", fontsize=9)
    fig.savefig(FIG / "fig05_chloride_maps_by_depth_interval.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    st = pd.DataFrame(stats)
    st.to_csv(TAB / "kriging_cross_validation.csv", index=False, float_format="%.3f")
    print(st.to_string())
    cnt = recent.groupby(["interval", "side"]).agg(wells=("wqp_id", "size"), median_cl=("chloride_mgL", "median"),
                                                  over_250=("chloride_mgL", lambda v: (v > 250).sum()))
    cnt.to_csv(TAB / "chloride_latest_by_interval_side.csv", float_format="%.2f")
    print(cnt)


if __name__ == "__main__":
    main()
