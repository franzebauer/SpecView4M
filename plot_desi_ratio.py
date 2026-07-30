#!/usr/bin/env python3
"""
plot_desi_ratio.py — 4MOST / DESI flux-ratio percentiles vs wavelength,
optionally split by a per-spectrum property (r-band SNR or r magnitude).

For every LR spectrum with a DESI DR1 match, interpolate both spectra onto a
common 1 Å grid, smooth each on a 10 Å scale (Gaussian, FWHM = 10 Å), and form
the ratio 4MOST / DESI.  The per-spectrum ratios are cached; the plotter bins
them and shows, per bin, the wavelength-resolved median and the 14/86, 5/95 and
1/99 percentile envelopes, plus a panel overlaying the per-bin medians.

Splits:
  --split snr    r-band SNR (SNR_R) bins 0-1, 1-2, 2-4, 4-8, >8
  --split rmag   r magnitude (CAT_MAG) bins 17-18 … 21-22
  --seeing-max X keep only spectra with DIMM seeing (FWHM_AMBI) < X arcsec

Outputs: PLOTS/desi_ratio_<split>.png
         PLOTS/desi_ratio_cache.npz   (per-spectrum ratios + filenames)

Usage:
  python plot_desi_ratio.py --split rmag --seeing-max 1.0
  python plot_desi_ratio.py --split snr
  python plot_desi_ratio.py --recompute      # rebuild the cache
"""

import os, sys, re, time, argparse
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
from viewer import (read_spectrum, read_desi_spectrum, build_desi_index,
                    build_alt_index, _coord_key, CSV_DEFAULT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── configuration ───────────────────────────────────────────────────────────
GRID        = np.arange(3700.0, 9401.0, 1.0)   # common wavelength grid (Å)
SMOOTH_FWHM = 10.0                             # smoothing scale (Å)
SIGMA_PIX   = SMOOTH_FWHM / 2.3548 / (GRID[1] - GRID[0])
MIN_DESI_FLUX = 0.0
MIN_COUNT     = 8          # need ≥ this many spectra at a wavelength to plot
PCTS = [1, 5, 14, 50, 86, 95, 99]
NORM_WAVE = 5500.0         # per-object normalisation wavelength (Å)
NORM_HALF = 25.0           # half-width of the normalisation window (Å)
OUTDIR = os.path.join(BASE, "PLOTS")
COND_CSV = os.path.join(BASE, "conditions.csv")


def _cache_path(recon):
    return os.path.join(OUTDIR, f"desi_ratio_cache"
                        f"{'_alt' if recon == 'alt' else ''}.npz")

SPLITS = {
    "none": (None,    [(-np.inf, np.inf, "all")], "all spectra"),
    "snr":  ("SNR_R", [(0, 1, "0–1"), (1, 2, "1–2"), (2, 4, "2–4"),
                       (4, 8, "4–8"), (8, np.inf, ">8")], "r-band SNR"),
    "rmag": ("rmag",  [(17, 18, "17–18"), (18, 19, "18–19"), (19, 20, "19–20"),
                       (20, 21, "20–21"), (21, 22, "21–22")], "r magnitude"),
}


def normalize_at(ratios, wave=NORM_WAVE, half=NORM_HALF):
    """Divide each spectrum's ratio by its median value in a window at `wave`
    (so ratio == 1 there); rows with no valid data in the window → all-NaN."""
    win = (GRID >= wave - half) & (GRID <= wave + half)
    with np.errstate(all="ignore"):
        norm = np.nanmedian(ratios[:, win], axis=1)
    out = ratios / norm[:, None]
    out[~np.isfinite(norm) | (norm <= 0)] = np.nan
    return out


def _specuid(fn):
    m = re.match(r"qmost_[^_]+_(\d+)_", os.path.basename(str(fn)))
    return m.group(1) if m else None


def _gauss_kernel(sigma_pix):
    hw = int(4 * sigma_pix) + 1
    x  = np.arange(-hw, hw + 1)
    k  = np.exp(-0.5 * (x / sigma_pix) ** 2)
    return k / k.sum()


_KERNEL = _gauss_kernel(SIGMA_PIX)


def _smooth(y):
    m   = np.isfinite(y).astype(float)
    yy  = np.where(np.isfinite(y), y, 0.0)
    num = np.convolve(yy, _KERNEL, mode="same")
    den = np.convolve(m,  _KERNEL, mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[den < 0.5] = np.nan
    return out


def _regrid(wave, flux, extra_mask=None):
    m = np.isfinite(wave) & np.isfinite(flux)
    if extra_mask is not None:
        m &= extra_mask
    if m.sum() < 10:
        return None
    return np.interp(GRID, wave[m], flux[m], left=np.nan, right=np.nan)


def compute_ratios(csv_path, recon="primary"):
    """Return (ratios[Nspec×Ngrid] float32, filenames[Nspec]).

    recon='primary' uses the LR/ files; recon='alt' uses the alternate
    reduction in LR4/ (matched by coord), keeping the same LR∩DESI sample."""
    df = pd.read_csv(csv_path, low_memory=False)
    csv_dir = os.path.dirname(os.path.abspath(csv_path))
    desi_index = build_desi_index()
    alt_lr = build_alt_index().get("LR", {}) if recon == "alt" else {}
    lr = df[df["category"] == "LR"].copy()
    lr["coord"] = lr["filename"].apply(_coord_key)
    lr = lr[lr["coord"].map(lambda k: k in desi_index)].reset_index(drop=True)
    print(f"{len(lr)} LR spectra with DESI matches "
          f"(reduction: {recon})")

    ratios = np.full((len(lr), GRID.size), np.nan, dtype=np.float32)
    names  = lr["filename"].to_numpy()
    t0, kept = time.time(), 0
    for i, row in lr.iterrows():
        if recon == "alt":
            fp = alt_lr.get(row["coord"])
            if not fp:
                continue
        else:
            fp = str(row["filepath"])
            if not os.path.isabs(fp):
                fp = os.path.join(csv_dir, fp)
        w, f, e, _sky = read_spectrum(fp)
        if w is None:
            continue
        wd, fd = read_desi_spectrum(desi_index[row["coord"]])
        if wd is None:
            continue
        f4 = _regrid(w, f, extra_mask=(np.isfinite(e) & (e > 0)))
        fD = _regrid(wd, fd)
        if f4 is None or fD is None:
            continue
        f4s, fDs = _smooth(f4), _smooth(fD)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = f4s / fDs
        r[~np.isfinite(fDs) | (fDs <= MIN_DESI_FLUX)] = np.nan
        ratios[i] = r
        kept += 1
        if (i + 1) % 500 == 0:
            print(f"  {i + 1}/{len(lr)}  ({(i + 1) / (time.time() - t0):.0f}/s)",
                  flush=True)
    print(f"Used {kept} pairs in {time.time() - t0:.0f}s")
    return ratios, names


def load_or_compute(csv_path, recompute, recon="primary"):
    cache = _cache_path(recon)
    if os.path.exists(cache) and not recompute:
        z = np.load(cache, allow_pickle=True)
        if ("filenames" in z.files and z["grid"].shape == GRID.shape
                and np.allclose(z["grid"], GRID)):
            print(f"Loaded cached ratios {z['ratios'].shape} from "
                  f"{os.path.basename(cache)} (use --recompute to rebuild)")
            return z["ratios"], z["filenames"]
    ratios, names = compute_ratios(csv_path, recon=recon)
    os.makedirs(OUTDIR, exist_ok=True)
    np.savez_compressed(cache, ratios=ratios,
                        filenames=np.array(names, dtype=object), grid=GRID)
    print(f"Cached ratios → {os.path.basename(cache)}")
    return ratios, names


def build_meta(filenames, csv_path):
    """Per-spectrum SNR_R, rmag (CAT_MAG), seeing (DIMM FWHM_AMBI)."""
    df = pd.read_csv(csv_path, low_memory=False)
    df = df.drop_duplicates("filename").set_index("filename")
    meta = pd.DataFrame(index=pd.Index(filenames, name="filename"))
    meta["SNR_R"] = pd.to_numeric(df.reindex(filenames)["SNR_R"].values,
                                  errors="coerce")
    meta["rmag"]  = pd.to_numeric(df.reindex(filenames)["CAT_MAG"].values,
                                  errors="coerce")
    seeing = np.full(len(filenames), np.nan)
    if os.path.exists(COND_CSV):
        c = pd.read_csv(COND_CSV, low_memory=False)
        c["specuid"] = c["filename"].apply(_specuid)
        smap = (c.dropna(subset=["specuid"]).drop_duplicates("specuid")
                 .set_index("specuid")["FWHM_AMBI"])
        seeing = pd.to_numeric(
            pd.Series([_specuid(fn) for fn in filenames]).map(smap),
            errors="coerce").to_numpy()
    else:
        print(f"  Note: {os.path.basename(COND_CSV)} not found — no seeing cut")
    meta["seeing"] = seeing
    return meta.reset_index(drop=True)


def percentiles(ratios):
    count = np.isfinite(ratios).sum(axis=0)
    with np.errstate(all="ignore"):
        Q = np.nanpercentile(ratios, PCTS, axis=0)
    Q[:, count < MIN_COUNT] = np.nan
    return Q


def _bands(ax, Q):
    q1, q5, q14, q50, q86, q95, q99 = Q
    ax.fill_between(GRID, q1,  q99, color="#4C72B0", alpha=0.15, lw=0,
                    label="1–99%")
    ax.fill_between(GRID, q5,  q95, color="#4C72B0", alpha=0.28, lw=0,
                    label="5–95%")
    ax.fill_between(GRID, q14, q86, color="#4C72B0", alpha=0.45, lw=0,
                    label="14–86%")
    ax.plot(GRID, q50, color="#C44E52", lw=1.5, label="median")


def plot_single(ratios, sel, title, out, ylim, ylabel, norm_wave=None):
    fig, ax = plt.subplots(figsize=(13, 6.5))
    fig.patch.set_facecolor("white")
    n = int(sel.sum())
    _bands(ax, percentiles(ratios[sel]))
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.7)
    for wsp in (5540, 6900):
        ax.axvline(wsp, color="grey", lw=1.0, ls=":", alpha=0.6)
    if norm_wave is not None:
        ax.axvline(norm_wave, color="#2ca02c", lw=1.0, ls="-", alpha=0.6)
    ax.set_yscale("log"); ax.set_ylim(*ylim); ax.set_xlim(GRID[0], GRID[-1])
    ax.set_xlabel("Observed wavelength (Å)"); ax.set_ylabel(ylabel)
    ax.set_title(f"{title}   (N = {n})")
    ax.legend(loc="upper right", framealpha=0.9); ax.grid(alpha=0.2, which="both")
    fig.tight_layout(); fig.savefig(out, dpi=140)
    print(f"Saved {out}")


def plot_binned(ratios, values, bins, split_label, out, base_mask=None,
                mask_note="", ylim=(0.1, 10.0),
                ylabel="Flux ratio  4MOST / DESI", suptitle=None):
    fig, axes = plt.subplots(3, 2, figsize=(15, 12), sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    axes = axes.ravel()
    cmap = plt.get_cmap("viridis")
    med_colors = [cmap(x) for x in np.linspace(0.05, 0.9, len(bins))]

    medians = []
    for k, (lo, hi, lab) in enumerate(bins):
        ax = axes[k]
        sel = np.isfinite(values) & (values >= lo) & (values < hi)
        if base_mask is not None:
            sel &= base_mask
        sel &= np.isfinite(ratios).any(axis=1)   # only spectra with data
        n = int(sel.sum())
        if n == 0:
            ax.text(0.5, 0.5, f"{lab}: no spectra", ha="center", va="center",
                    transform=ax.transAxes)
            medians.append(np.full(GRID.size, np.nan)); continue
        Q = percentiles(ratios[sel])
        medians.append(Q[3])
        _bands(ax, Q)
        ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.7)
        ax.set_title(f"{split_label} {lab}   (N = {n})")
        if k == 0:
            ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    axo = axes[len(bins)]
    for (lo, hi, lab), med, c in zip(bins, medians, med_colors):
        axo.plot(GRID, med, color=c, lw=1.6, label=f"{split_label} {lab}")
    axo.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.7)
    axo.set_title("Median ratio by bin")
    axo.legend(loc="lower right", fontsize=8, framealpha=0.9)

    for ax in axes:
        ax.set_yscale("log"); ax.set_ylim(*ylim)
        ax.set_xlim(GRID[0], GRID[-1])
        for wsp in (5540, 6900):
            ax.axvline(wsp, color="grey", lw=1.0, ls=":", alpha=0.6)
        ax.grid(alpha=0.2, which="both")
    for ax in axes[-2:]:
        ax.set_xlabel("Observed wavelength (Å)")
    for ax in (axes[0], axes[2], axes[4]):
        ax.set_ylabel(ylabel)

    if suptitle is None:
        suptitle = (f"4MOST / DESI DR1 flux ratio by {split_label}{mask_note} "
                    f"(smoothed on {SMOOTH_FWHM:g} Å scale)")
    fig.suptitle(suptitle, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out, dpi=140)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default=CSV_DEFAULT)
    ap.add_argument("--split", choices=list(SPLITS), default="rmag")
    ap.add_argument("--recon", choices=["primary", "alt"], default="primary",
                    help="which 4MOST reduction to use (alt = LR4/)")
    ap.add_argument("--seeing-max", type=float, default=None,
                    help="keep only spectra with DIMM seeing < this (arcsec)")
    ap.add_argument("--normalize", action="store_true",
                    help=f"normalise each ratio to 1 at the norm wavelength first")
    ap.add_argument("--norm-wave", type=float, default=NORM_WAVE,
                    help=f"normalisation wavelength in Å (default {NORM_WAVE:g})")
    ap.add_argument("--recompute", action="store_true")
    args = ap.parse_args()

    ratios, filenames = load_or_compute(args.csv, args.recompute, args.recon)
    meta = build_meta(filenames, args.csv)

    col, bins, label = SPLITS[args.split]
    values = meta[col].to_numpy() if col else np.zeros(len(filenames))

    recon_tag  = "_alt" if args.recon == "alt" else ""
    recon_note = "  [alt reduction]" if args.recon == "alt" else ""
    base_mask, note, tag = None, recon_note, args.split
    if args.seeing_max is not None:
        base_mask = np.isfinite(meta["seeing"].to_numpy()) & \
                    (meta["seeing"].to_numpy() < args.seeing_max)
        note = f", seeing < {args.seeing_max:g}\"" + recon_note
        tag  = f"{args.split}_see{args.seeing_max:g}"
        print(f"seeing cut < {args.seeing_max:g}\": "
              f"{int(base_mask.sum())}/{len(filenames)} spectra")

    # normalisation: divide each ratio by its value at the norm wavelength
    ylim, ylabel, norm_wave = (0.1, 10.0), "Flux ratio  4MOST / DESI", None
    if args.normalize:
        ratios = normalize_at(ratios, wave=args.norm_wave)
        ylim   = (0.3, 3.0)
        ylabel = f"Flux ratio  4MOST / DESI  (norm @ {args.norm_wave:g} Å)"
        norm_wave = args.norm_wave
        note += f", norm@{args.norm_wave:g}Å"
        tag  += f"_norm{args.norm_wave:g}"
    tag += recon_tag

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"desi_ratio_{tag}.png")

    if args.split == "none":
        sel = np.isfinite(ratios).any(axis=1)
        if base_mask is not None:
            sel &= base_mask
        title = (f"4MOST / DESI DR1 flux ratio{note} "
                 f"(smoothed on {SMOOTH_FWHM:g} Å scale)")
        plot_single(ratios, sel, title, out, ylim, ylabel, norm_wave=norm_wave)
    else:
        suptitle = (f"4MOST / DESI DR1 flux ratio by {label}{note} "
                    f"(smoothed on {SMOOTH_FWHM:g} Å scale)")
        plot_binned(ratios, values, bins, label, out, base_mask=base_mask,
                    ylim=ylim, ylabel=ylabel, suptitle=suptitle)


if __name__ == "__main__":
    main()
