"""Step 7. ILLUSTRATIVE saltwater-intrusion calculations for the "2,000-foot" sand.

These are back-of-envelope calculations meant to show orders of magnitude. They are not
a calibrated model and should not be used for management decisions. The USGS SEAWAT
and MODFLOW models (SIR 2013-5227, SIR 2019-5102) are the appropriate tools for that.

1. Two-end-member mixing. Saline fraction f = (Cl - Cl_fresh) / (Cl_saline - Cl_fresh).
   Cl_fresh: median latest chloride at 2,000-foot sand wells more than 8 km north of the fault.
   Cl_saline is poorly known, so two values are used: (a) the median of all samples since
   1990 at 2,000-foot sand wells south of the fault, and (b) the highest chloride recorded
   in the sand on either side of the fault. Fractions are capped at 1.
   Specific conductance vs chloride is checked for the linear relation that simple
   mixing of two waters implies.
2. Hydraulic gradient north of the fault: a plane fitted to mean 2015 to 2025 water-level
   altitudes at 2,000-foot sand wells north of the fault.
3. Advective front speed v = K i / n_e with ASSUMED ranges of hydraulic conductivity K and
   effective porosity n_e (stated below; not taken from the cited reports), and the travel
   time to cover 3 miles (the distance between the fault and industrial wells given in SIR 2013-5227).
4. Apparent front speed from data: the first year chloride exceeded 50 mg/L at 2,000-foot sand
   wells north of the fault that started below 50 mg/L, regressed on distance from the fault.
   Only two wells qualify in the current data, so this is a two-point estimate.
5. One-dimensional advection-dispersion (Ogata-Banks) breakthrough curves using the apparent
   speed and a longitudinal dispersivity of 10 percent of travel distance.
"""
import json

import numpy as np
import pandas as pd
from scipy.special import erfc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import PROC, TAB, FIG

SAND = "12220BR"
K_FT_D = (30.0, 100.0, 200.0)        # assumed hydraulic conductivity range, ft/d
NE = (0.35, 0.25, 0.20)              # assumed effective porosity (paired low to high speed)
THRESH = 50.0                        # mg/L, front arrival threshold
DIST_3MI_M = 3 * 1609.344


def main():
    d = pd.read_csv(PROC / "chloride_spcond_samples.csv", dtype={"aqfr_cd": str}, parse_dates=["date"])
    w = d[d.aqfr_cd == SAND].copy()
    cl = w.dropna(subset=["chloride_mgL"])
    res = {}

    # 1. mixing
    latest = cl.sort_values("date").groupby("wqp_id").tail(1)
    fresh = latest[(latest.side == "north") & (latest.dist_fault_km > 8)]
    saline = cl[(cl.side == "south") & (cl.year >= 1990)]
    c_f, c_s = float(fresh.chloride_mgL.median()), float(saline.chloride_mgL.median())
    c_max = float(cl.chloride_mgL.max())
    res.update(cl_fresh_mgL=c_f, n_fresh_wells=len(fresh), cl_saline_median_south_mgL=c_s,
               n_saline_samples=len(saline), n_saline_wells=saline.wqp_id.nunique(), cl_saline_max_in_sand_mgL=c_max)
    near = latest[(latest.side == "north") & (latest.dist_fault_km <= 8) & (latest.year >= 2000)].copy()
    near["saline_fraction_median_south"] = ((near.chloride_mgL - c_f) / (c_s - c_f)).clip(0, 1)
    near["saline_fraction_max_in_sand"] = ((near.chloride_mgL - c_f) / (c_max - c_f)).clip(0, 1)
    near[["wqp_id", "station_nm", "dist_fault_km", "date", "chloride_mgL", "saline_fraction_median_south",
          "saline_fraction_max_in_sand"]].sort_values(
        "dist_fault_km").to_csv(TAB / "mixing_fraction_2000ft_sand_near_fault.csv", index=False, float_format="%.4g")
    res["near_fault_wells_since2000"] = len(near)
    res["near_fault_wells_over_1pct_saline_maxEM"] = int((near.saline_fraction_max_in_sand > 0.01).sum())

    both = d.dropna(subset=["chloride_mgL", "spcond_uScm"])
    both = both[both.chloride_mgL > 0]
    b = np.polyfit(both.spcond_uScm, both.chloride_mgL, 1)
    r2 = np.corrcoef(both.spcond_uScm, both.chloride_mgL)[0, 1] ** 2
    res.update(sc_cl_n=len(both), sc_cl_slope=float(b[0]), sc_cl_intercept=float(b[1]), sc_cl_r2=float(r2))

    # 2. gradient north of fault
    ann = pd.read_csv(PROC / "gwlevels_annual.csv", dtype={"site_no": str})
    wells = pd.read_csv(PROC / "wells.csv", dtype={"site_no": str})
    side = d.drop_duplicates("site_no").set_index("site_no")[["x_utm", "y_utm", "side", "dist_fault_km"]]
    lv = ann[(ann.year >= 2015)].merge(wells[["site_no", "aqfr_cd"]], on="site_no")
    lv = lv[lv.aqfr_cd == SAND].groupby("site_no").wl_alt_ft.mean().to_frame().join(side, how="inner")
    lv = lv[lv.side == "north"]
    A = np.column_stack([np.ones(len(lv)), lv.x_utm, lv.y_utm])
    coef, *_ = np.linalg.lstsq(A, lv.wl_alt_ft.values * 0.3048, rcond=None)
    grad = float(np.hypot(coef[1], coef[2]))
    res.update(gradient_wells=len(lv), gradient_m_per_m=grad,
               head_range_ft=[float(lv.wl_alt_ft.min()), float(lv.wl_alt_ft.max())])

    # 3. Darcy front speed with assumed parameters
    rows = []
    for k, n in zip(K_FT_D, NE):
        v = (k * 0.3048) * grad / n  # m/d
        rows.append(dict(K_ft_per_day=k, eff_porosity=n, v_m_per_yr=v * 365.25,
                         years_to_travel_3_miles=DIST_3MI_M / (v * 365.25)))
    darcy = pd.DataFrame(rows)

    # 4. apparent front speed
    arr = []
    for wid, g in cl[cl.side == "north"].groupby("wqp_id"):
        g = g.sort_values("date")
        hit = g[g.chloride_mgL >= THRESH]
        if len(hit) and g.chloride_mgL.iloc[0] < THRESH:
            arr.append(dict(wqp_id=wid, station_nm=g.station_nm.iloc[0], dist_km=g.dist_fault_km.iloc[0],
                            first_sample=int(g.year.min()), first_year_over_50=int(hit.year.iloc[0])))
    arr = pd.DataFrame(arr).sort_values("dist_km")
    arr.to_csv(TAB / "front_arrival_2000ft_sand.csv", index=False, float_format="%.3f")
    if len(arr) >= 2:
        s = np.polyfit(arr.first_year_over_50, arr.dist_km, 1)
        res["apparent_front_speed_m_per_yr"] = float(s[0] * 1000)
        res["apparent_front_wells"] = len(arr)
    darcy.to_csv(TAB / "darcy_front_speed_illustrative.csv", index=False, float_format="%.3g")
    print(darcy.round(2).to_string())
    print(arr.to_string())

    # 5. Ogata-Banks curves
    v = max(res.get("apparent_front_speed_m_per_yr", darcy.v_m_per_yr.iloc[1]), 1.0)
    t = np.linspace(0.1, 150, 600)
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    for x_m in (500, 1000, 2000, DIST_3MI_M):
        D = 0.1 * x_m * v
        c = 0.5 * (erfc((x_m - v * t) / (2 * np.sqrt(D * t))) + np.exp(v * x_m / D) * erfc((x_m + v * t) / (2 * np.sqrt(D * t))))
        c = np.nan_to_num(c, nan=1.0)
        ax.plot(t, np.clip(c, 0, 1), label=f"{x_m / 1000:.1f} km north of fault")
    ax.set_xlabel("Years after saline water starts crossing the fault")
    ax.set_ylabel("Relative concentration C/C0")
    ax.set_title(f"ILLUSTRATIVE 1-D advection-dispersion (Ogata-Banks), v = {v:.0f} m/yr, dispersivity = 0.1 x distance", fontsize=8.5)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig06_illustrative_intrusion_front.png", dpi=150)
    plt.close(fig)
    res["ogata_banks_v_m_per_yr"] = v
    (TAB / "intrusion_illustration_summary.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
