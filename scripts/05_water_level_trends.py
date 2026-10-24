"""Step 5. Long-term water-level trends per sand and piecewise regression for recovery.

1. Per-well Mann-Kendall test (Hamed and Rao variance correction for serial
   correlation) and Sen slope on annual mean water-level altitude, full record and 2000 onward.
2. Per-sand composite hydrograph: a two-way fixed-effects fit
   wl[well, year] = a[well] + b[year] on all wells in the sand, so the year
   effects b[year] are not biased by wells entering and leaving the network.
   The year effects are expressed relative to their own maximum. Years with
   only one measured well are dropped.
3. Piecewise (segmented) linear regression on each composite with 0, 1 or 2
   breakpoints chosen by BIC. Breakpoints are found by grid search with at least
   6 years in each segment. A 90% interval for the single-breakpoint year comes from
   a moving-block bootstrap of residuals (block length 5 years, 500 replicates).
"""
import numpy as np
import pandas as pd
import pymannkendall as mk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import PROC, TAB, FIG, SANDS, MAIN_SANDS

RNG = np.random.default_rng(42)
MIN_YEARS, MIN_SPAN = 10, 15


def mk_sen(y, x):
    if len(y) < MIN_YEARS:
        return None
    try:
        r = mk.hamed_rao_modification_test(y)
    except Exception:
        r = mk.original_test(y)
    slope = mk.sens_slope(y).slope  # per index step
    # Sen slope on the actual years (gaps allowed)
    xs, ys = np.asarray(x, float), np.asarray(y, float)
    i, j = np.triu_indices(len(xs), 1)
    sl = np.median((ys[j] - ys[i]) / (xs[j] - xs[i]))
    return dict(trend=r.trend, p=r.p, tau=r.Tau, sen_ft_per_yr=sl)


def well_trends(ann, wells):
    rows = []
    for sid, g in ann.groupby("site_no"):
        g = g.sort_values("year")
        for label, sub in (("full", g), ("2000_on", g[g.year >= 2000])):
            if len(sub) < MIN_YEARS or sub.year.max() - sub.year.min() < MIN_SPAN:
                continue
            r = mk_sen(sub.wl_alt_ft.values, sub.year.values)
            rows.append(dict(site_no=sid, period=label, n_years=len(sub),
                             first=sub.year.min(), last=sub.year.max(), **r))
    t = pd.DataFrame(rows).merge(wells[["site_no", "station_nm", "aqfr_cd", "sand", "well_depth_ft", "lat", "lon"]], on="site_no")
    return t


def composite(sub):
    """Two-way fixed effects year index from an unbalanced well-year panel."""
    obs = sub[["year", "site_no", "wl_alt_ft"]].dropna().rename(columns={"wl_alt_ft": "y"})
    years, sites = np.sort(obs.year.unique()), np.sort(obs.site_no.unique())
    yi = obs.year.map({y: k for k, y in enumerate(years)}).values
    si = obs.site_no.map({s: k for k, s in enumerate(sites)}).values
    n, ny, ns = len(obs), len(years), len(sites)
    X = np.zeros((n, ny + ns - 1))
    X[np.arange(n), yi] = 1.0
    m = si > 0
    X[np.arange(n)[m], ny + si[m] - 1] = 1.0
    beta, *_ = np.linalg.lstsq(X, obs.y.values, rcond=None)
    nwell = obs.groupby("year").size().reindex(years).values
    return pd.DataFrame({"year": years, "index_ft": beta[:ny] - beta[:ny].max(), "n_wells": nwell})


def seg_design(x, bps):
    cols = [np.ones_like(x), x]
    cols += [np.maximum(x - b, 0) for b in bps]
    return np.column_stack(cols)


def fit_seg(x, y, k, min_seg=6):
    cand = x[(x >= x.min() + min_seg) & (x <= x.max() - min_seg)]
    best = (np.inf, None, None)
    if k == 0:
        X = seg_design(x, [])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(((y - X @ b) ** 2).sum()), [], b
    if k == 1:
        for b1 in cand:
            X = seg_design(x, [b1])
            b, *_ = np.linalg.lstsq(X, y, rcond=None)
            s = ((y - X @ b) ** 2).sum()
            if s < best[0]:
                best = (s, [b1], b)
    if k == 2:
        for b1 in cand:
            for b2 in cand[cand >= b1 + min_seg]:
                X = seg_design(x, [b1, b2])
                b, *_ = np.linalg.lstsq(X, y, rcond=None)
                s = ((y - X @ b) ** 2).sum()
                if s < best[0]:
                    best = (s, [b1, b2], b)
    return float(best[0]), best[1], best[2]


def bic(rss, n, k):
    p = 2 + 2 * k  # intercept, slope, plus slope change and breakpoint for each break
    return n * np.log(rss / n) + p * np.log(n)


def block_boot_bp(x, y, b_hat, coef, B=500, L=5):
    fit = seg_design(x, b_hat) @ coef
    res = y - fit
    n = len(y)
    out = []
    for _ in range(B):
        idx = np.concatenate([np.arange(s, s + L) for s in RNG.integers(0, n - L + 1, int(np.ceil(n / L)))])[:n]
        yb = fit + res[idx]
        out.append(fit_seg(x, yb, 1)[1][0])
    return np.percentile(out, [5, 95])


def main():
    ann = pd.read_csv(PROC / "gwlevels_annual.csv", dtype={"site_no": str})
    wells = pd.read_csv(PROC / "wells.csv", dtype={"site_no": str})
    ann = ann.merge(wells[["site_no", "aqfr_cd"]], on="site_no")

    wt = well_trends(ann, wells)
    wt.to_csv(TAB / "well_trends_mann_kendall_sen.csv", index=False, float_format="%.4g")
    summ = (wt.groupby(["sand", "period"]).agg(
        wells=("site_no", "size"),
        rising_sig=("trend", lambda s: (s == "increasing").sum()),
        falling_sig=("trend", lambda s: (s == "decreasing").sum()),
        median_sen_ft_per_yr=("sen_ft_per_yr", "median")).reset_index())
    order = {SANDS[c]: k for k, c in enumerate(SANDS)}
    summ["o"] = summ.sand.map(order)
    summ = summ.sort_values(["period", "o"], ascending=[False, True]).drop(columns="o")
    summ.to_csv(TAB / "trend_summary_by_sand.csv", index=False, float_format="%.3f")
    print(summ.to_string())

    comps, segs = [], []
    for code in MAIN_SANDS:
        sub = ann[ann.aqfr_cd == code]
        # keep wells with at least 5 years so each well ties into the panel
        keep = sub.groupby("site_no").year.transform("size") >= 5
        sub = sub[keep]
        if sub.site_no.nunique() < 2:
            continue
        c = composite(sub)
        c = c[c.n_wells >= 2]  # drop years measured at only one well
        c["sand"] = SANDS[code]
        comps.append(c)
        x, y = c.year.values.astype(float), c.index_ft.values
        if len(x) < 25:
            continue
        fits = {k: fit_seg(x, y, k) for k in (0, 1, 2)}
        bics = {k: bic(fits[k][0], len(x), k) for k in fits}
        kbest = min(bics, key=bics.get)
        rss1, bp1, coef1 = fits[1]
        ci = block_boot_bp(x, y, bp1, coef1)
        r_all = mk_sen(y, x)
        rec = dict(sand=SANDS[code], n_wells=sub.site_no.nunique(), first=int(x.min()), last=int(x.max()),
                   n_years=len(x), mk_trend=r_all["trend"], mk_p=r_all["p"], sen_ft_per_yr=r_all["sen_ft_per_yr"],
                   bic_k0=bics[0], bic_k1=bics[1], bic_k2=bics[2], best_k=kbest,
                   bp1_year=bp1[0], bp1_ci90_lo=ci[0], bp1_ci90_hi=ci[1],
                   slope_before_bp1=coef1[1], slope_after_bp1=coef1[1] + coef1[2],
                   bp2_years=";".join(str(int(b)) for b in fits[2][1]),
                   slopes_k2=";".join(f"{s:.2f}" for s in np.cumsum(fits[2][2][1:])))
        segs.append(rec)
    comp = pd.concat(comps)
    comp.to_csv(PROC / "composite_hydrograph_by_sand.csv", index=False, float_format="%.3f")
    seg = pd.DataFrame(segs)
    seg.to_csv(TAB / "piecewise_regression_by_sand.csv", index=False, float_format="%.3f")
    print(seg[["sand", "n_wells", "first", "last", "sen_ft_per_yr", "best_k", "bp1_year", "bp1_ci90_lo", "bp1_ci90_hi",
               "slope_before_bp1", "slope_after_bp1", "bp2_years", "slopes_k2"]].round(2).to_string())

    # figure: composite hydrographs with fitted segments
    sands = seg.sand.tolist()
    nc = 2
    nr = int(np.ceil(len(sands) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(10, 2.3 * nr), sharex=True)
    for ax, s in zip(axes.flat, sands):
        c = comp[comp.sand == s]
        r = seg[seg.sand == s].iloc[0]
        ax.plot(c.year, c.index_ft, ".", color="0.35", ms=3)
        k = int(r.best_k)
        x = c.year.values.astype(float)
        _, bp, coef = fit_seg(x, c.index_ft.values, k)
        ax.plot(x, seg_design(x, bp) @ coef, color="#1f5fa6", lw=1.6)
        for b in bp:
            ax.axvline(b, color="#c0392b", lw=0.8, ls="--")
        ax.set_title(f"{s} ({int(r.n_wells)} wells), breakpoints: {', '.join(str(int(b)) for b in bp) or 'none'}", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
    for ax in axes.flat[len(sands):]:
        ax.axis("off")
    fig.supylabel("Water-level index (ft, relative to series maximum)", fontsize=9)
    fig.suptitle("Composite water-level hydrographs by sand, East Baton Rouge Parish (USGS NWIS)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig02_composite_hydrographs_piecewise.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
