#!/usr/bin/env python3
"""
viewer.py  –  Interactive 4MOST spectra viewer
=======================================================
Features
--------
  • Tabs           : arms (LR|HR|LRD|HIZ) + subsurveys (M_SED|W_SED|VAR_L|
                     VAR_Z|VAR_G|W_HIZ) + DESI-matched
  • HR view        : combined HJ14 file (all 3 arms in one FITS)
  • Batch view     : grid of spectra per page (1–16, default 9)
  • Smoothing      : Gaussian kernel slider (σ = 0 … 20 px)
  • Sky scale      : divide sky flux by 1–10 via slider
  • Navigation     : ◀ Prev / Next ▶ + page-jump text box + ← → keyboard
  • Statistics     : SNR vs mag, redshift histogram, sky map, etc.
  • Tagging         : per-page=1 shows a quality (ok/unclear/bad minor/bad major)
                     + issue + redshift-override sidebar; Save tags →
                     spectra_tags.csv, reloaded at startup
  • Reductions     : default (blue) = the L1 combination (LR4/HIZ4/… dirs; the
                     panel is blank where no …4 file exists); "alt" overlay
                     (orange) = my combination (LR)
  • DESI overlay   : matched DESI DR1 spectrum (pink) via SPV_DESI_match.fits
  • Spectral theme : gray=raw  blue=smooth  green=sky  red=err  pink=DESI
                     orange=alt (my LR combination)

Keyboard shortcuts
------------------
  ← / →       previous / next page
  1 … 9       switch tab (arm / subsurvey / DESI)
  0           toggle statistics dashboard
  s           cycle smoothing (+2 px each press)
  q / Esc     quit
"""

import os, sys, gzip, time
import numpy as np
import pandas as pd

# ── backend selection ──────────────────────────────────────────────────────
import matplotlib
_interactive = False
for _backend in ["TkAgg", "MacOSX", "Qt5Agg", "Qt6Agg", "GTK4Agg", "GTK3Agg"]:
    try:
        matplotlib.use(_backend, force=True)
        import matplotlib.pyplot as plt
        _fig_test = plt.figure()
        plt.close(_fig_test)
        _interactive = True
        break
    except Exception:
        pass
if not _interactive:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    print("INFO: No interactive display — saving PNGs instead.")

import matplotlib.gridspec as gridspec
import matplotlib.ticker
from matplotlib.widgets import Button, Slider, TextBox

# ── paths ──────────────────────────────────────────────────────────────────
BASE        = os.path.dirname(os.path.abspath(__file__))
CSV_DEFAULT = os.path.join(BASE, "spectra_headers.csv")

# ── LRS arm splice wavelengths (Å) ─────────────────────────────────────────
LRS_SPLICE_WAVES = [5540, 6900]
LRS_CATS         = {"LR", "LRD", "HIZ"}


# ── Filter effective wavelengths (Å) ───────────────────────────────────────
FILTER_WAVE_AA = {
    "SDSS_U_AB": 3551, "SDSS_G_AB": 4686, "SDSS_R_AB": 6166,
    "SDSS_I_AB": 7480, "SDSS_Z_AB": 8932,
    "SDSS_u_AB": 3551, "SDSS_g_AB": 4686, "SDSS_r_AB": 6166,
    "SDSS_i_AB": 7480, "SDSS_z_AB": 8932,
    "DECam_u_AB": 3567, "DECam_g_AB": 4811, "DECam_r_AB": 6417,
    "DECam_i_AB": 7838, "DECam_z_AB": 9259, "DECam_Y_AB": 9888,
    "decam_u_AB": 3567, "decam_g_AB": 4811, "decam_r_AB": 6417,
    "decam_i_AB": 7838, "decam_z_AB": 9259, "decam_Y_AB": 9888,
    "DECAM_G_AB_PSF": 4811, "DECAM_R_AB_PSF": 6417, "DECAM_I_AB_PSF": 7838,
    "DECAM_Z_AB_PSF": 9259, "DECAM_Y_AB_PSF": 9888,
    "LSST_u_AB": 3671, "LSST_g_AB": 4827, "LSST_r_AB": 6233,
    "LSST_i_AB": 7542, "LSST_z_AB": 8691, "LSST_y_AB": 9712,
    "GAIA_G_Vega": 6740, "GAIA_BP_Vega": 5320, "GAIA_RP_Vega": 7970,
    "Gaia_G_Vega": 6740, "Gaia_BP_Vega": 5320, "Gaia_RP_Vega": 7970,
}
_FILTER_PARTIALS = [
    ("sdss_u", 3551), ("sdss_g", 4686), ("sdss_r", 6166),
    ("sdss_i", 7480), ("sdss_z", 8932),
    ("decam_u", 3567), ("decam_g", 4811), ("decam_r", 6417),
    ("decam_i", 7838), ("decam_z", 9259),
    ("lsst_u", 3671), ("lsst_g", 4827), ("lsst_r", 6233),
    ("lsst_i", 7542), ("lsst_z", 8691), ("lsst_y", 9712),
    ("gaia_bp", 5320), ("gaia_rp", 7970), ("gaia_g", 6740),
    ("_u",  3551), ("_g",  4686), ("_r",  6166),
    ("_i",  7480), ("_z",  8932), ("_y",  9712),
]
_MAG_ID_WAVE = {
    "u": 3551, "g": 4686, "r": 6166, "i": 7480, "z": 8932, "y": 9712,
    "G": 6740, "B": 5320, "R": 7970, "Y": 9888,
}

def _filter_wave(mag_type_or_id):
    if not mag_type_or_id:
        return None
    mt = str(mag_type_or_id).strip()
    if mt in FILTER_WAVE_AA:
        return FILTER_WAVE_AA[mt]
    if len(mt) == 1 and mt in _MAG_ID_WAVE:
        return _MAG_ID_WAVE[mt]
    ml = mt.lower()
    for key, wave in _FILTER_PARTIALS:
        if key in ml:
            return wave
    return None

def _arm_color_for_wave(wave_aa):
    if wave_aa is None:
        return "#aaaaaa"
    if wave_aa < LRS_SPLICE_WAVES[0]:
        return "#7ecbff"
    if wave_aa < LRS_SPLICE_WAVES[1]:
        return "#7dffb3"
    return "#ffaa88"

# ── colour palette ─────────────────────────────────────────────────────────
CAT_COLORS = {
    "LR":  "#2196F3",
    "HR":  "#E91E63",
    "LRD": "#4CAF50",
    "HIZ": "#FF9800",
}
CAT_LABELS = {
    "LR":  "Low-Res (LR)",
    "HR":  "High-Res (HR)",
    "LRD": "Low-Res Diffuse (LRD)",
    "HIZ": "High-z (HIZ)",
}

# Spectral plot colours
C_RAW    = "#c8c8c8"   # very light gray  – raw spectrum
C_SMOOTH = "#5599ff"   # blue             – smoothed spectrum
C_SKY    = "#44cc55"   # green            – sky (offset below zero)
C_ERR    = "#ff5555"   # red              – error spectrum
C_DESI   = "#ff66cc"   # pink             – DESI DR1 overlay
C_ALT    = "#ff9800"   # orange           – "alt" overlay = my combination (LR/…)
C_EM     = "#00e5ff"   # cyan             – emission-line markers
C_ABS    = "#ffd54f"   # amber            – absorption-line markers

# ── Reductions ──────────────────────────────────────────────────────────────
#   default (blue trace) = the L1 combination in <category>4 dirs (LR4/HIZ4/…);
#     the panel is blank ("No LR4 spectrum") where a source has no <category>4 file;
#   "alt" overlay        = my combination in <category> dirs (LR/HIZ/…, i.e. the
#     spectra_headers.csv `filepath`), drawn in orange when the L1 default exists.
L1_SUFFIX = "4"

# ── Rest-frame spectral lines (Å) — major AGN + galaxy features ─────────────
EM_LINES = [
    (1215.67, "Lyα"),  (1240.1, "NV"),    (1399.8, "SiIV"),  (1549.1, "CIV"),
    (1640.4, "HeII"),  (1908.7, "CIII]"), (2798.8, "MgII"),  (3426.0, "[NeV]"),
    (3727.5, "[OII]"), (3869.0, "[NeIII]"), (4101.7, "Hδ"),  (4340.5, "Hγ"),
    (4861.3, "Hβ"),    (4958.9, "[OIII]"), (5006.8, "[OIII]"), (5875.6, "HeI"),
    (6562.8, "Hα"),    (6583.5, "[NII]"), (6716.4, "[SII]"),
]
ABS_LINES = [
    (3933.7, "CaK"), (3968.5, "CaH"), (4304.4, "G"),
    (5175.3, "Mgb"), (5892.9, "NaD"), (8542.1, "CaT"),
]

# ── Quality tagging ─────────────────────────────────────────────────────────
TAGS_CSV     = os.path.join(BASE, "spectra_tags.csv")
QUALITIES    = ["ok", "unclear", "bad minor", "bad major"]
QUAL_COLORS  = {"ok": "#4CAF50", "unclear": "#FF9800",
                "bad minor": "#FF7043", "bad major": "#f44336"}
# legacy quality values → new scheme
QUAL_MIGRATE = {"good": "ok", "bad": "bad minor"}
ISSUE_DEFS   = [   # (key, button label)
    ("blue_shape", "blue cont. shape"),
    ("neg_cont",   "negative cont"),
    ("bg_join",    "B/G join"),
    ("gr_join",    "G/R join"),
    ("bumpy",      "unusual cont (bumpy)"),
    ("norm",       "unusual cont (norm)"),
    ("telluric",   "telluric lines"),
    ("gap",        "gap"),
    ("incorrect_z", "incorrect redshift"),
]

# ── DESI overlay ───────────────────────────────────────────────────────────
import re as _re
DESI_DIR           = os.path.join(BASE, "DESI_spectra")
DESI_MATCH_FITS    = os.path.join(BASE, "SPV_DESI_match.fits")
DESI_RADIUS_ARCSEC = 3.0

# ── xpca redshift catalogues (one per arm; joined by filename == PROV) ──────
XPCA_FILES = ["LR_xpca.fits", "LRD_xpca.fits", "HIZ_xpca.fits", "HR_xpca.fits"]

# ── Subsurvey categories (via S16 target catalogue) ─────────────────────────
S16_CAT                = os.path.join(BASE, "S16_20250429T0456Z_target_catalogue.fits.gz")
SUBSURVEY_CACHE        = os.path.join(BASE, ".subsurvey_cache.csv")
SUBSURVEY_MATCH_ARCSEC = 1.0
# catalogue SUBSURVEY value → viewer tab label, in display order
SUBSURVEY_VIEWS = [
    ("M_SED", "M_SED"), ("W_SED", "W_SED"), ("W_VARL", "VAR_L"),
    ("W_VARZ", "VAR_Z"), ("W_VARG", "VAR_G"), ("W_HIZ", "W_HIZ"),
]
SUBSURVEY_COLORS = {
    "M_SED": "#9C27B0", "W_SED": "#00BCD4", "W_VARL": "#8BC34A",
    "W_VARZ": "#FFC107", "W_VARG": "#FF5722", "W_HIZ": "#795548",
}

# ── FITS spectrum reader ───────────────────────────────────────────────────

def _parse_cards(block):
    cards, end_found = {}, False
    for i in range(0, min(len(block), 2880), 80):
        card = block[i:i+80].decode("ascii", errors="replace")
        key  = card[:8].strip()
        rest = card[8:]
        if key == "END":
            end_found = True
            break
        if "=" in rest[:2]:
            vc = rest[2:]
            val = (vc.split("'")[1].strip() if "'" in vc
                   else vc.split("/")[0].strip())
            if key:
                cards[key] = val
    return cards, end_found


def read_spectrum(filepath):
    """
    Read WAVE, FLUX, ERR, SKY from a 4MOST L1 FITS file.
    SKY is derived as FLUX_NOSS − FLUX (sky contribution).
    Returns (wave, flux, err, sky) as float32 arrays, or (None,…) on error.
    """
    try:
        opener = gzip.open if filepath.endswith(".gz") else open
        with opener(filepath, "rb") as fh:
            in_primary = True
            nelem      = None
            while True:
                block = fh.read(2880)
                if not block:
                    return None, None, None, None
                cards, end_found = _parse_cards(block)
                if "NELEM" in cards:
                    try:
                        nelem = int(cards["NELEM"])
                    except ValueError:
                        pass
                if end_found:
                    if in_primary:
                        in_primary = False
                        continue
                    if nelem is None:
                        return None, None, None, None
                    ne  = nelem
                    # Column layout (big-endian, stored column-major):
                    #   WAVE     : ne × float32  (bytes 0   : 4ne)
                    #   FLUX     : ne × float32  (bytes 4ne : 8ne)
                    #   ERR      : ne × float32  (bytes 8ne : 12ne)
                    #   QUAL     : ne × int16    (bytes 12ne: 14ne)
                    #   FLUX_NOSS: ne × float32  (bytes 14ne: 18ne)
                    #   ERR_NOSS : ne × float32  (bytes 18ne: 22ne)
                    #   SENSFUNC : ne × float32  (HR only,   22ne: 26ne)
                    full_bytes = ne * 22
                    raw = fh.read(full_bytes)
                    if len(raw) < ne * 12:
                        return None, None, None, None
                    wave = np.frombuffer(raw[     : 4*ne], dtype=">f4").astype(np.float32)
                    flux = np.frombuffer(raw[4*ne : 8*ne], dtype=">f4").astype(np.float32)
                    err  = np.frombuffer(raw[8*ne :12*ne], dtype=">f4").astype(np.float32)
                    # Sky: FLUX_NOSS − FLUX
                    sky = None
                    if len(raw) >= 18 * ne:
                        fnoss = np.frombuffer(raw[14*ne:18*ne], dtype=">f4").astype(np.float32)
                        sky_raw = fnoss - flux
                        finite_sky = sky_raw[np.isfinite(sky_raw)]
                        if len(finite_sky) > 10 and np.any(finite_sky != 0):
                            sky = sky_raw
                    # Insert NaN breaks at large wavelength gaps (HR inter-arm gaps)
                    dw = np.diff(wave)
                    gap_idx = np.where(dw > 50.0)[0] + 1  # first pixel of each new arm
                    if len(gap_idx):
                        nan1 = np.full(1, np.nan, dtype=np.float32)
                        for gi in sorted(gap_idx, reverse=True):
                            mid_w = np.array([(wave[gi-1] + wave[gi]) / 2.0],
                                             dtype=np.float32)
                            wave = np.concatenate([wave[:gi], mid_w, wave[gi:]])
                            flux = np.concatenate([flux[:gi], nan1,  flux[gi:]])
                            err  = np.concatenate([err[:gi],  nan1,  err[gi:]])
                            if sky is not None:
                                sky = np.concatenate([sky[:gi], nan1, sky[gi:]])
                    return wave, flux, err, sky
    except Exception:
        return None, None, None, None


def gaussian_smooth(flux, sigma):
    """Gaussian smoothing using pure numpy."""
    if sigma <= 0:
        return flux.copy()
    hw     = int(4 * sigma) + 1
    x      = np.arange(-hw, hw + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel /= kernel.sum()
    n   = len(flux)
    pad = np.pad(flux, hw, mode="edge")
    out = np.convolve(pad, kernel, mode="valid")
    return out[:n]


# ── DESI overlay helpers ────────────────────────────────────────────────────

def _coord_key(filename):
    """First field of a qmost_<coord>_... filename, used to join across runs."""
    m = _re.match(r"qmost_([^_]+)_", str(filename))
    return m.group(1) if m else None


def read_desi_spectrum(json_path):
    """Return (wave, flux) float32 arrays from a DESI JSON file.
    DESI flux is in 1e-17 erg/s/cm²/Å; scale to match 4MOST CGS units."""
    import json
    try:
        with open(json_path) as f:
            data = json.load(f)
        wave = np.array([x if x is not None else np.nan for x in data["wave"]],
                        dtype=np.float32)
        flux = np.array([x if x is not None else np.nan for x in data["flux"]],
                        dtype=np.float32) * 1e-17
        return wave, flux
    except Exception:
        return None, None


def build_desi_index(match_fits=DESI_MATCH_FITS, desi_dir=DESI_DIR):
    """coord_key → best DESI json path, from the SPV_DESI_match catalogue."""
    index = {}
    if not (os.path.exists(match_fits) and os.path.isdir(desi_dir)):
        return index
    try:
        from astropy.table import Table
        t = Table.read(match_fits)
    except Exception as exc:
        print(f"  Note: could not read {os.path.basename(match_fits)} "
              f"({exc}) — no DESI overlays")
        return index

    best_sep = {}
    for prov, tid, sep in zip(t["PROV"], t["TARGETID"], t["Separation"]):
        coord = _coord_key(prov)
        if coord is None:
            continue
        try:
            sep = float(sep)
        except (TypeError, ValueError):
            sep = 999.0
        if sep > DESI_RADIUS_ARCSEC:
            continue
        p = os.path.join(desi_dir, f"targetid_{int(tid)}.json")
        if not os.path.exists(p):
            continue
        if coord not in best_sep or sep < best_sep[coord]:
            best_sep[coord] = sep
            index[coord]    = p

    if index:
        print(f"  {len(index)} DESI spectra indexed from "
              f"{os.path.basename(match_fits)}")
    return index


def build_l1_index(base=BASE, cats=("LR", "HIZ", "LRD", "HR")):
    """{category: {coord_key: path}} for the L1-combination <category>4 dirs.

    This is the *default* reduction (LR4/HIZ4/…); it shares each source's
    coordinate string but has a different SPECUID/date, so it is matched by
    coord (as with DESI)."""
    import glob
    index = {}
    for cat in cats:
        d = os.path.join(base, cat + L1_SUFFIX)
        if not os.path.isdir(d):
            continue
        cmap = {}
        for p in sorted(glob.glob(os.path.join(d, "qmost_*.fits"))):
            k = _coord_key(os.path.basename(p))
            if k:
                cmap.setdefault(k, p)
        if cmap:
            index[cat] = cmap
    if index:
        print("  L1 combination (default): " +
              ", ".join(f"{c}{L1_SUFFIX}={len(m)}" for c, m in index.items()))
    return index


def assign_subsurvey(df, cat_path=S16_CAT, cache_path=SUBSURVEY_CACHE):
    """Per-spectrum catalogue SUBSURVEY via 1″ match to the S16 target catalogue.

    Returns an object array (len == len(df)) of SUBSURVEY strings ("NONE" where
    unmatched), or None if the catalogue/astropy is unavailable. Result is cached
    to a filename-keyed sidecar CSV so only the first launch pays the match cost.
    """
    names = df["filename"].astype(str)

    if os.path.exists(cache_path):
        try:
            cache = pd.read_csv(cache_path).set_index("filename")["subsurvey"]
            if names.isin(cache.index).all():
                return names.map(cache).fillna("NONE").to_numpy()
        except Exception:
            pass

    if not os.path.exists(cat_path):
        print(f"  Note: {os.path.basename(cat_path)} not found — no subsurvey views")
        return None
    try:
        from astropy.io import fits
        from astropy.coordinates import SkyCoord
        import astropy.units as u
    except Exception as exc:
        print(f"  Note: astropy unavailable ({exc}) — no subsurvey views")
        return None

    ra  = pd.to_numeric(df["RA"],  errors="coerce").to_numpy()
    dec = pd.to_numeric(df["DEC"], errors="coerce").to_numpy()
    ok  = np.isfinite(ra) & np.isfinite(dec)
    out = np.full(len(df), "NONE", dtype=object)

    t0 = time.time()
    with fits.open(cat_path) as h:
        c    = h[1].data
        css  = np.array(c["SUBSURVEY"]).astype(str)
        catc = SkyCoord(np.array(c["RA"], float) * u.deg,
                        np.array(c["DEC"], float) * u.deg)
    obs = SkyCoord(ra[ok] * u.deg, dec[ok] * u.deg)
    idx, sep, _ = obs.match_to_catalog_sky(catc)
    hit = sep.arcsec < SUBSURVEY_MATCH_ARCSEC
    out[np.where(ok)[0]] = np.where(hit, css[idx], "NONE")
    print(f"  subsurveys assigned via {os.path.basename(cat_path)} "
          f"({int(hit.sum())}/{int(ok.sum())} matched, {time.time()-t0:.1f}s)")

    try:
        pd.DataFrame({"filename": names, "subsurvey": out}).to_csv(cache_path,
                                                                   index=False)
    except Exception:
        pass
    return out


def load_xpca_redshifts(files=XPCA_FILES, base=BASE):
    """Return a DataFrame indexed by filename with 4XP redshift estimates
    (XPCA_Z, XPCA_ZPROB, XPCA_ZTYPE), joined by PROV == L1 filename."""
    try:
        from astropy.io import fits
    except Exception:
        return None
    prov, z, zp, zt = [], [], [], []
    found = False
    for f in files:
        path = f if os.path.isabs(f) else os.path.join(base, f)
        if not os.path.exists(path):
            continue
        found = True
        with fits.open(path) as h:
            d = h[1].data
            prov.extend(np.array(d["PROV"]).astype(str))
            z.extend(np.array(d["zBest"], dtype=float))
            zp.extend(np.array(d["zBestProb"], dtype=float))
            zt.extend(np.array(d["zBestType"]).astype(str))
    if not found:
        print("  Note: no *_xpca.fits found — redshifts unavailable")
        return None
    tbl = pd.DataFrame({"XPCA_Z": z, "XPCA_ZPROB": zp, "XPCA_ZTYPE": zt},
                       index=pd.Index(prov, name="filename"))
    tbl = tbl[~tbl.index.duplicated(keep="first")]
    print(f"  {len(tbl):,} xpca redshifts loaded from {len(files)} catalogues")
    return tbl


def _empty_tag():
    return {"quality": None, "issues": set(), "z": None}


def load_tags(path=TAGS_CSV):
    """filename → {'quality': str|None, 'issues': set, 'z': float|None}."""
    tags = {}
    if not os.path.exists(path):
        return tags
    try:
        t = pd.read_csv(path)
    except Exception:
        return tags
    for _, r in t.iterrows():
        fn = str(r.get("filename", "")).strip()
        if not fn:
            continue
        q = r.get("quality")
        q = str(q).strip() if (q is not None and str(q).strip()
                               and str(q).lower() != "nan") else None
        q = QUAL_MIGRATE.get(q, q)   # remap legacy good/bad values
        iss = str(r.get("issues") or "")
        issues = {x for x in iss.split("|") if x and x != "nan"}
        try:
            z = float(r.get("redshift"))
            if not np.isfinite(z):
                z = None
        except (TypeError, ValueError):
            z = None
        tags[fn] = {"quality": q, "issues": issues, "z": z}
    print(f"  {len(tags):,} spectra tags loaded from {os.path.basename(path)}")
    return tags


def save_tags(tags, path=TAGS_CSV):
    """Write non-empty tag entries to the tag CSV. Returns the row count."""
    rows = []
    for fn, d in tags.items():
        q   = d.get("quality")
        iss = d.get("issues") or set()
        z   = d.get("z")
        if not q and not iss and z is None:
            continue
        rows.append({"filename": fn, "quality": q or "",
                     "issues": "|".join(sorted(iss)),
                     "redshift": "" if z is None else f"{z:.5f}"})
    pd.DataFrame(rows, columns=["filename", "quality", "issues",
                                "redshift"]).to_csv(path, index=False)
    return len(rows)


# ── drawing helpers ────────────────────────────────────────────────────────

def _fv(v, d=2):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "–"
    try:
        return f"{float(v):.{d}f}"
    except (TypeError, ValueError):
        return str(v)

def _style_ax(ax):
    ax.set_facecolor("#16213e")
    for sp in ax.spines.values():
        sp.set_color("#333")
    ax.tick_params(colors="white", labelsize=8, direction="in", length=2)
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")
    ax.title.set_color("white")


def draw_one_spectrum(ax, wave, flux, err, sky, row, sigma, cat_color,
                      sky_scale=1.0,
                      show_raw=True, show_smooth=True,
                      show_sky=True, show_err=True,
                      desi_wave=None, desi_flux=None, show_desi=True,
                      desi_sigma=None, show_lines=True, z_override=None,
                      alt_wave=None, alt_flux=None, show_alt=True,
                      title_fontsize=13, empty_msg="Read error"):
    """
    Plot a single spectrum (or concatenated HR arms) with the colour scheme:
      light gray = raw spectrum
      blue       = smoothed spectrum (the L1 default)
      green      = sky / sky_scale  (displayed at natural vertical position)
      red        = error spectrum (positive, 1σ)
      orange     = "alt" overlay = my LR combination

    LRS arm splice lines and photometric magnitude markers are overlaid.
    """
    ax.cla()
    _style_ax(ax)

    if wave is None:
        ax.text(0.5, 0.5, empty_msg, ha="center", va="center",
                transform=ax.transAxes, color="#f44", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        _set_title(ax, row, cat_color, title_fontsize)
        return

    mask = (np.isfinite(flux) & np.isfinite(err)
            & (err > 0) & (np.abs(flux) < 1e12))
    if mask.sum() < 5:
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="#fa0")
        _set_title(ax, row, cat_color, title_fontsize)
        return

    # --- valid data ---
    w_all = wave        # may contain NaN gaps (HR arms)
    f_all = flux
    e_all = err

    w = wave[mask]; f = flux[mask]; e = err[mask]

    lo, hi = np.nanpercentile(f, 2), np.nanpercentile(f, 98)
    if lo == hi:
        lo -= 1; hi += 1
    frange = hi - lo
    pad    = frange * 0.12

    fs = gaussian_smooth(f, sigma)    # smooth only valid points

    # Need to smooth the full (gapped) array for display
    # Build a smoothed version of the whole array (NaN-aware)
    fs_all = np.full_like(f_all, np.nan)
    fs_all[mask] = fs

    # ── sky display (green, scaled down by sky_scale) ────────────────────
    sky_display = None
    if sky is not None:
        sky_valid = sky[mask]
        sky_finite = sky_valid[np.isfinite(sky_valid)]
        if len(sky_finite) > 10 and np.nanstd(sky_finite) > 0:
            scale = max(float(sky_scale), 1.0)
            sky_display = np.full_like(f_all, np.nan)
            sky_display[mask] = sky_valid / scale

    # ── error display (red line) ─────────────────────────────────────────
    err_display = np.full_like(f_all, np.nan)
    err_display[mask] = e

    # ── y-range (accommodate sky and err) ───────────────────────────────
    ymin = lo - pad
    ymax = hi + pad * 2

    if sky_display is not None:
        sky_plot_vals = sky_display[np.isfinite(sky_display)]
        if len(sky_plot_vals):
            ymin = min(ymin, np.nanpercentile(sky_plot_vals, 2) - frange * 0.05)

    e_p98 = np.nanpercentile(e, 98) if len(e) else 0
    ymax  = max(ymax, e_p98 * 1.08)

    # ── x-range (handle NaN gaps in HR arms) ────────────────────────────
    w_finite = w_all[np.isfinite(w_all)]
    xlo, xhi = w_finite[0], w_finite[-1]

    # ── plots ─────────────────────────────────────────────────────────────
    ax.axhline(0, color="#555555", lw=0.4, alpha=0.6, zorder=0)

    # Sky (behind everything)
    if sky_display is not None and show_sky:
        ax.plot(w_all, sky_display,
                color=C_SKY, alpha=0.65, linewidth=0.5, zorder=1)

    # Error
    if show_err:
        ax.plot(w_all, err_display,
                color=C_ERR, alpha=0.55, linewidth=0.45, zorder=2)

    # Raw spectrum
    if show_raw:
        ax.plot(w_all, f_all,
                color=C_RAW, alpha=0.65, linewidth=0.4, zorder=3)

    # Smoothed spectrum (on top)
    if show_smooth:
        ax.plot(w_all, fs_all,
                color=C_SMOOTH, alpha=0.92, linewidth=0.9, zorder=4)

    # "alt" overlay = my combination (LR), orange, smoothed with the same kernel
    if show_alt and alt_wave is not None and alt_flux is not None:
        amask = np.isfinite(alt_flux) & (np.abs(alt_flux) < 1e12)
        if amask.sum() >= 5:
            afs = np.full_like(alt_flux, np.nan)
            afs[amask] = gaussian_smooth(alt_flux[amask], sigma)
            ax.plot(alt_wave, afs,
                    color=C_ALT, alpha=0.9, linewidth=0.9, zorder=5,
                    label="alt: LR (mine)")
            ax.text(0.985, 0.90, "alt ✓", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6, color=C_ALT, zorder=8)

    # DESI DR1 overlay (pink, smoothed with the same kernel)
    if show_desi and desi_wave is not None and desi_flux is not None:
        dmask = np.isfinite(desi_flux) & (np.abs(desi_flux) < 1e12)
        if dmask.sum() >= 5:
            _ds = sigma if desi_sigma is None else desi_sigma
            dfs = np.full_like(desi_flux, np.nan)
            dfs[dmask] = gaussian_smooth(desi_flux[dmask], _ds)
            ax.plot(desi_wave, dfs,
                    color=C_DESI, alpha=0.85, linewidth=0.8, zorder=5,
                    label="DESI DR1")
            ax.text(0.985, 0.96, "DESI ✓", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6, color=C_DESI,
                    zorder=8)

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("λ (Å)", fontsize=7)
    ax.set_ylabel("Flux", fontsize=7)

    # ── LRS arm splice markers ─────────────────────────────────────────
    cat = str(row.get("category", ""))
    if cat in LRS_CATS:
        for wsp in LRS_SPLICE_WAVES:
            if xlo < wsp < xhi:
                ax.axvline(wsp, color="#ffffff", lw=1.8, ls="--", alpha=0.55)

    # ── Photometric magnitude markers ─────────────────────────────────
    _C_AA_S = 2.998e18

    def _ab_to_flam(mag_ab, wave_aa):
        f_nu = 10.0 ** (-0.4 * (mag_ab + 48.6))
        return f_nu * _C_AA_S / wave_aa ** 2

    _phot_points = []
    for _col_mag, _col_id, _label in [
        ("CAL_MAG_BLUE",  "CAL_MAG_ID_BLUE",  "B"),
        ("CAL_MAG_GREEN", "CAL_MAG_ID_GREEN",  "G"),
        ("CAL_MAG_RED",   "CAL_MAG_ID_RED",    "R"),
    ]:
        _mag = row.get(_col_mag)
        _id  = str(row.get(_col_id) or "").strip()
        try:
            _mag = float(_mag)
        except (TypeError, ValueError):
            _mag = np.nan
        if np.isfinite(_mag) and _mag > 0 and _id:
            _wv = _filter_wave(_id)
            if _wv is not None:
                _phot_points.append((_wv, _ab_to_flam(_mag, _wv),
                                     f"{_label}={_mag:.1f}"))

    if not _phot_points:
        _mag_main = row.get("CAT_MAG")
        _mtype    = str(row.get("CAT_MAG_TYPE") or "").strip()
        try:
            _mag_main = float(_mag_main)
        except (TypeError, ValueError):
            _mag_main = np.nan
        if np.isfinite(_mag_main) and _mag_main > 0 and _mtype:
            _wv = _filter_wave(_mtype)
            if _wv is not None:
                _short = (_mtype.replace("_AB", "").replace("_Vega", "")
                          .replace("_PSF", ""))
                _phot_points.append((_wv, _ab_to_flam(_mag_main, _wv),
                                     f"{_short}={_mag_main:.1f}"))

    _ymin, _ymax = ax.get_ylim()
    _wspan = xhi - xlo
    _phot_in_range = [(wv, fl, lb) for (wv, fl, lb) in _phot_points
                      if xlo <= wv <= xhi]
    if _phot_in_range:
        _flam_min = min(fl for _, fl, _ in _phot_in_range)
        _flam_max = max(fl for _, fl, _ in _phot_in_range)
        _new_ymin = min(_ymin, _flam_min * 0.7)
        _new_ymax = max(_ymax, _flam_max * 1.4)
        if _new_ymin != _ymin or _new_ymax != _ymax:
            ax.set_ylim(_new_ymin, _new_ymax)
            _ymin, _ymax = _new_ymin, _new_ymax

    for _wv, _fl, _lbl in _phot_in_range:
        _mc = _arm_color_for_wave(_wv)
        ax.plot(_wv, _fl, marker="D", ms=5, color=_mc,
                mec="white", mew=0.5, zorder=6, alpha=0.9)
        ax.text(_wv + _wspan * 0.015, _fl, _lbl,
                color=_mc, fontsize=6, va="center", ha="left",
                zorder=7, alpha=0.9)

    # ── redshifted spectral-line markers ─────────────────────────────────
    if show_lines:
        zval = z_override if z_override is not None else row.get("XPCA_Z")
        draw_line_markers(ax, zval, xlo, xhi)

    _set_title(ax, row, cat_color, title_fontsize)


def draw_line_markers(ax, z, xlo, xhi):
    """Tick + label every major emission/absorption line at its observed
    wavelength rest·(1+z), along the bottom of the axes. No-op if z is invalid."""
    try:
        z = float(z)
    except (TypeError, ValueError):
        return
    if not np.isfinite(z):
        return
    ymin, ymax = ax.get_ylim()
    span     = ymax - ymin
    tick_top = ymin + span * 0.50   # extend ~halfway up the panel
    lab_y    = ymin + span * 0.03
    for lines, col, ls in ((EM_LINES, C_EM, "-"), (ABS_LINES, C_ABS, ":")):
        for rest, lab in lines:
            obs = rest * (1.0 + z)
            if xlo <= obs <= xhi:
                ax.vlines(obs, ymin, tick_top, color=col, lw=0.7,
                          ls=ls, alpha=0.75, zorder=6)
                ax.text(obs, lab_y, lab, color=col, fontsize=7, rotation=90,
                        ha="center", va="bottom", alpha=0.95, zorder=7)


def _set_title(ax, row, color, fontsize=13):
    sg  = _fv(row.get("SNR_G"), 1)
    sr  = _fv(row.get("SNR_R"), 1)
    si  = _fv(row.get("SNR_I"), 1)
    sz  = _fv(row.get("SNR_Z"), 1)
    # Redshift: prefer xpca (4XP) estimate, fall back to QZC_Z
    z_raw   = row.get("XPCA_Z")
    zprob   = row.get("XPCA_ZPROB")
    cls     = str(row.get("XPCA_ZTYPE") or "").strip()
    try:
        z_raw = float(z_raw)
    except (TypeError, ValueError):
        z_raw = np.nan
    if np.isfinite(z_raw):
        z = f"{z_raw:.3f}"
        try:
            if np.isfinite(float(zprob)):
                z += f" (p={float(zprob):.2f})"
        except (TypeError, ValueError):
            pass
    else:
        z   = _fv(row.get("QZC_Z"), 3)
        cls = cls or str(row.get("QSC_CLS") or "").strip()
    mag = _fv(row.get("CAT_MAG"), 1)
    mtype = str(row.get("CAT_MAG_TYPE") or "").strip().replace("_AB","").replace("_Vega","")
    sub   = str(row.get("CAT_SUBSURVEY") or "").strip()
    name  = str(row.get("filename", ""))
    # Exposure time
    try:
        texp_s = float(row.get("TEXPTIME") or 0)
        texp_str = f"  texp={texp_s:.0f}s" if texp_s > 0 else ""
    except (TypeError, ValueError):
        texp_str = ""
    ax.set_title(
        f"{name}\n"
        f"g/r/i/z SNR={sg}/{sr}/{si}/{sz}  z={z}  {cls}{texp_str}\n"
        f"{sub}  {mtype}={mag}",
        fontsize=fontsize, pad=3, color=color)


# ── statistics dashboard ───────────────────────────────────────────────────

def draw_stats_page(fig, df):
    fig.clf()
    fig.patch.set_facecolor("#1a1a2e")
    fig.suptitle("4MOST Spectra – Statistics Dashboard",
                 color="white", fontsize=13, y=0.97)

    gs = gridspec.GridSpec(3, 3, figure=fig,
                           left=0.07, right=0.97,
                           top=0.92, bottom=0.06,
                           hspace=0.58, wspace=0.38)

    cats = [c for c in ["LR", "HR", "LRD", "HIZ"]
            if c in df["category"].unique()]

    def sax(pos, title):
        ax = fig.add_subplot(pos)
        ax.set_facecolor("#16213e")
        ax.tick_params(colors="white", labelsize=8)
        for sp in ax.spines.values():
            sp.set_color("#444")
        ax.set_title(title, color="white", fontsize=9, pad=4)
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")
        return ax

    BAND_DEFS = [
        ("SNR_G", "SDSS-g  (4000–5540 Å)", "#4FC3F7"),
        ("SNR_R", "SDSS-r  (5540–6900 Å)", "#EF9A9A"),
        ("SNR_I", "SDSS-i  (6900–8350 Å)", "#A5D6A7"),
        ("SNR_Z", "SDSS-z  (8350–9500 Å)", "#FFE082"),
    ]
    for (bcol, blabel, bcolor), gpos in zip(BAND_DEFS,
                                             [gs[0,0], gs[0,1], gs[0,2], gs[1,0]]):
        ax = sax(gpos, f"SNR vs mag  [{blabel}]")
        for cat in cats:
            sub = df[df["category"] == cat]
            snr = pd.to_numeric(sub[bcol], errors="coerce")
            mag = pd.to_numeric(sub["CAT_MAG"], errors="coerce")
            ok  = np.isfinite(snr) & np.isfinite(mag) & (snr > 0)
            if ok.sum() > 0:
                ax.scatter(mag[ok], snr[ok], s=5, alpha=0.55,
                           color=CAT_COLORS[cat], label=cat, rasterized=True)
        ax.set_xlabel("Catalogue magnitude", fontsize=7)
        ax.set_ylabel("Median per-pixel SNR", fontsize=7)
        ax.set_yscale("log"); ax.set_ylim(0.5, 30)
        ax.axhline(1, color=bcolor, lw=0.8, ls="--", alpha=0.7)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.legend(fontsize=6, framealpha=0.3, facecolor="#333",
                  labelcolor="white", markerscale=2)

    ax = sax(gs[1, 1], "Measured Redshift Distribution")
    for cat in cats:
        sub = df[df["category"] == cat]
        z   = pd.to_numeric(sub["QZC_Z"], errors="coerce")
        ok  = np.isfinite(z) & (z > -0.5) & (z < 8)
        if ok.sum() > 0:
            ax.hist(z[ok], bins=50, histtype="step",
                    color=CAT_COLORS[cat], label=cat, linewidth=1.3)
    ax.set_xlabel("Redshift z"); ax.set_ylabel("Count")
    ax.legend(fontsize=7, framealpha=0.3, facecolor="#333", labelcolor="white")

    ax = sax(gs[1, 2], "Object Classification (pipeline)")
    cls_counts = df["QSC_CLS"].value_counts().dropna()
    if len(cls_counts):
        colors = plt.cm.Set2(np.linspace(0, 1, len(cls_counts)))
        _, _, autotexts = ax.pie(
            cls_counts.values, labels=cls_counts.index,
            autopct=lambda p: f"{p:.0f}%" if p > 3 else "",
            colors=colors,
            textprops={"color": "white", "fontsize": 7},
            pctdistance=0.78)
        for at in autotexts:
            at.set_fontsize(6)
        ax.legend([f"{k}: {v:,}" for k, v in cls_counts.items()],
                  loc="lower center", bbox_to_anchor=(0.5, -0.15),
                  fontsize=6, ncol=3, framealpha=0.3, facecolor="#333",
                  labelcolor="white")

    ax = sax(gs[2, 0], "Total Exposure Time")
    for cat in cats:
        sub = df[df["category"] == cat]
        t   = pd.to_numeric(sub["TEXPTIME"], errors="coerce") / 3600.0
        ok  = np.isfinite(t) & (t > 0)
        if ok.sum() > 0:
            ax.hist(t[ok], bins=30, histtype="step",
                    color=CAT_COLORS[cat], label=cat, linewidth=1.3)
    ax.set_xlabel("Total exposure time (hours)"); ax.set_ylabel("Count")
    ax.legend(fontsize=7, framealpha=0.3, facecolor="#333", labelcolor="white")

    ax = sax(gs[2, 1], "Sky Coverage (RA / Dec)")
    for cat in cats:
        sub = df[df["category"] == cat]
        ra  = pd.to_numeric(sub["RA"],  errors="coerce")
        dec = pd.to_numeric(sub["DEC"], errors="coerce")
        ok  = np.isfinite(ra) & np.isfinite(dec)
        if ok.sum() > 0:
            ax.scatter(ra[ok], dec[ok], s=2, alpha=0.4,
                       color=CAT_COLORS[cat], label=cat, rasterized=True)
    ax.set_xlabel("RA (deg)"); ax.set_ylabel("Dec (deg)")
    ax.legend(fontsize=7, framealpha=0.3, facecolor="#333",
              labelcolor="white", markerscale=3)
    ax.invert_xaxis()

    ax = sax(gs[2, 2], "Positional Offsets: Spectra − Catalogue")
    dra  = pd.to_numeric(df["DELTA_RA_AS"],  errors="coerce")
    ddec = pd.to_numeric(df["DELTA_DEC_AS"], errors="coerce")
    sep  = pd.to_numeric(df["SEP_ARCSEC"],   errors="coerce")
    for cat in cats:
        mask = (df["category"] == cat) & np.isfinite(dra) & np.isfinite(ddec)
        if mask.sum() > 0:
            ax.scatter(dra[mask], ddec[mask], s=4, alpha=0.45,
                       color=CAT_COLORS[cat], label=cat, rasterized=True)
    theta = np.linspace(0, 2 * np.pi, 200)
    for r, ls in [(0.1, "--"), (0.5, ":")]:
        ax.plot(r * np.cos(theta), r * np.sin(theta),
                color="white", lw=0.8, ls=ls, alpha=0.5, label=f'{r}"')
    ax.set_xlabel(r'$\Delta$RA · cos(Dec)  (arcsec)', fontsize=7)
    ax.set_ylabel(r'$\Delta$Dec  (arcsec)', fontsize=7)
    ax.set_aspect("equal", adjustable="datalim")
    ax.axhline(0, color="white", lw=0.4, alpha=0.3)
    ax.axvline(0, color="white", lw=0.4, alpha=0.3)
    n_beyond_half = int((sep > 0.5).sum())
    n_beyond_01   = int((sep > 0.1).sum())
    ax.text(0.02, 0.98,
            f'>{0.5}": {n_beyond_half}/{len(df)}\n'
            f'>{0.1}": {n_beyond_01}/{len(df)}\n'
            f'max: {sep.max():.3f}"',
            transform=ax.transAxes, fontsize=6.5, color="white",
            va="top", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="#0d1117", alpha=0.7))
    ax.legend(fontsize=6, framealpha=0.3, facecolor="#333",
              labelcolor="white", markerscale=2, loc="lower right")

    fig.canvas.draw_idle()


# ── viewer state ───────────────────────────────────────────────────────────

class ViewerState:
    def __init__(self, df, desi_index=None, l1_index=None):
        self.df         = df
        self.views      = self._build_views(df)   # list of (key, label, color, mask)
        self.view       = 0                        # index into self.views
        self.page       = 0
        self.n_per_page = 9
        self.sigma      = 7.0
        self.sigma_desi = 5.0
        self.sky_scale  = 1.0
        self.show_raw    = False
        self.show_smooth = True
        self.show_sky    = False
        self.show_err    = False
        self.show_desi   = True
        self.show_lines  = True
        self.show_alt    = True
        self.tab        = "spectra"
        self._cache     = {}     # filepath → (wave, flux, err, sky)
        self.desi_index = desi_index or {}   # coord_key → DESI json path
        self._desi_cache = {}    # json path → (wave, flux)
        self.l1_index   = l1_index or {}     # {category: {coord_key: LR4 path}}
        self.tags       = load_tags()        # filename → {quality, issues, z}
        self.tags_dirty = False
        self.use_custom_z = False            # single-mode redshift override
        self.qual_filter = None              # None or a quality string

    def get_tag(self, fn):
        return self.tags.get(fn, _empty_tag())

    def set_quality(self, fn, q):
        """Set quality, or clear it if the same value is re-selected."""
        cur = self.tags.setdefault(fn, _empty_tag())
        cur["quality"] = None if cur["quality"] == q else q
        self.tags_dirty = True

    def toggle_issue(self, fn, key):
        cur = self.tags.setdefault(fn, _empty_tag())
        issues = cur["issues"]
        issues.discard(key) if key in issues else issues.add(key)
        self.tags_dirty = True

    def set_redshift(self, fn, z):
        """Store a custom redshift (or None to clear) for this spectrum."""
        cur = self.tags.setdefault(fn, _empty_tag())
        cur["z"] = None if z is None else float(z)
        self.tags_dirty = True

    @staticmethod
    def _build_views(df):
        """Ordered list of (key, label, color, boolean-mask) tab filters:
        resolution arms, then populated subsurveys, then a DESI-match view."""
        views = []
        for cat in ["LR", "HR", "LRD", "HIZ"]:
            m = (df["category"] == cat)
            if m.any():
                views.append((cat, cat, CAT_COLORS.get(cat, "#555"), m))
        if "subsurvey" in df.columns:
            for ss, label in SUBSURVEY_VIEWS:
                m = (df["subsurvey"] == ss)
                if m.any():
                    views.append((ss, label, SUBSURVEY_COLORS.get(ss, "#607D8B"), m))
        if "has_desi" in df.columns:
            m = df["has_desi"].astype(bool)
            if m.any():
                views.append(("DESI", "DESI", C_DESI, m))
        return views

    @property
    def view_color(self):
        return self.views[self.view][2]

    def _effective_df(self):
        sub = self.df[self.views[self.view][3]]
        if self.qual_filter:
            keep = {fn for fn, d in self.tags.items()
                    if d.get("quality") == self.qual_filter}
            sub = sub[sub["filename"].isin(keep)]
        return sub.reset_index(drop=True)

    def n_pages(self):
        return max(1, int(np.ceil(len(self._effective_df()) / self.n_per_page)))

    def page_rows(self):
        edf   = self._effective_df()
        start = self.page * self.n_per_page
        return edf.iloc[start : start + self.n_per_page]

    def _get_raw(self, filepath):
        if filepath not in self._cache:
            if len(self._cache) > 80:
                self._cache.pop(next(iter(self._cache)))
            self._cache[filepath] = read_spectrum(filepath)
        return self._cache[filepath]

    def get_spectrum(self, filepath):
        """Return (wave, flux, err, sky) for a spectrum."""
        return self._get_raw(filepath)

    def get_desi(self, filename):
        """Return (wave, flux) of the matched DESI spectrum, or (None, None)."""
        path = self.desi_index.get(_coord_key(filename))
        if not path:
            return None, None
        if path not in self._desi_cache:
            if len(self._desi_cache) > 80:
                self._desi_cache.pop(next(iter(self._desi_cache)))
            self._desi_cache[path] = read_desi_spectrum(path)
        return self._desi_cache[path]

    def l1_path(self, category, filename):
        """Path to the L1-combination (default, LR4/…) file, or None."""
        return self.l1_index.get(category, {}).get(_coord_key(filename))

    def get_default(self, category, filename, fp):
        """(wave, flux, err, sky, is_l1) for the default (blue) trace = the L1
        combination (LR4/…). Returns all-None when this source has no L1 file,
        so the panel shows nothing rather than a blue non-L1 spectrum."""
        p = self.l1_path(category, filename)
        if p:
            return (*self._get_raw(p), True)
        return (None, None, None, None, False)

    def get_overlay(self, fp):
        """(wave, flux) of the 'alt' overlay = my combination (LR CSV filepath)."""
        w, f, _e, _s = self._get_raw(fp)
        return w, f


# ── grid layout helper ─────────────────────────────────────────────────────

def _grid_dims(n):
    if n <= 1:  return 1, 1
    if n <= 2:  return 1, 2
    if n <= 3:  return 1, 3
    if n <= 4:  return 2, 2
    if n <= 6:  return 2, 3
    if n <= 9:  return 3, 3
    if n <= 12: return 3, 4
    return 4, 4


# ── main viewer ────────────────────────────────────────────────────────────

class FourMostViewer:
    def __init__(self, df, desi_index=None, l1_index=None):
        self.state     = ViewerState(df, desi_index=desi_index,
                                     l1_index=l1_index)
        self.fig       = plt.figure(figsize=(16, 10))
        self.fig.patch.set_facecolor("#1a1a2e")
        self._install_resize_bug_guard()
        self.spec_axes = []
        self._build_spectra_ui()
        self._draw_spectra()
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    def _install_resize_bug_guard(self):
        """Silence a matplotlib bug where a ResizeEvent (no `.inaxes`) reaches a
        widget mouse-handler and spams `AttributeError: … 'inaxes'`. Only that
        specific error is swallowed; everything else keeps its normal handler."""
        cb = self.fig.canvas.callbacks
        prev = getattr(cb, "exception_handler", None)

        def _handler(exc):
            if isinstance(exc, AttributeError) and "inaxes" in str(exc):
                return
            if prev is not None:
                return prev(exc)
            raise exc

        cb.exception_handler = _handler

    # ── UI construction ────────────────────────────────────────────────

    def _build_spectra_ui(self):
        self.fig.clf()
        self.spec_axes = []
        st = self.state

        nr, nc = _grid_dims(st.n_per_page)

        # ── spectra grid (top bar for filters, bottom 22% for controls) ─
        #    single-spectrum mode shrinks the grid to make room for the
        #    tagging sidebar on the right.
        self._single = (st.n_per_page == 1)
        grid_right   = 0.79 if self._single else 0.98
        # more vertical gap for denser grids so the 3-line titles clear the
        # axis labels of the row above
        hspace = 0.42 + 0.14 * max(0, nr - 2)
        self._gs = gridspec.GridSpec(
            nr, nc,
            left=0.04, right=grid_right,
            top=0.84,  bottom=0.22,
            hspace=hspace, wspace=0.22)
        for i in range(nr * nc):
            ax = self.fig.add_subplot(self._gs[i // nc, i % nc])
            _style_ax(ax)
            self.spec_axes.append(ax)

        # ── layout constants ──────────────────────────────────────────
        # Row A: category buttons + nav + per-page + stats  (y=0.14)
        # Row B: smooth slider (narrow)                     (y=0.075)
        # Row C: colour legend                              (y=0.035)
        rA_y = 0.14;  rA_h = 0.040
        rB_y = 0.075; rB_h = 0.028

        # ── filter tabs (top bar): arms + subsurveys + DESI ───────────
        self._view_btns = []
        n_v   = max(1, len(st.views))
        fb_x0, fb_x1, fb_h = 0.035, 0.965, 0.032
        fb_y  = 0.925
        gap   = 0.006
        w     = min(0.090, (fb_x1 - fb_x0) / n_v - gap)
        for i, (key, label, color, _mask) in enumerate(st.views):
            bx  = fb_x0 + i * (w + gap)
            axb = self.fig.add_axes([bx, fb_y, w, fb_h])
            active = (i == st.view)
            btn = Button(axb, label,
                         color=color if active else "#2a2a2a",
                         hovercolor=color)
            btn.label.set_fontsize(8)
            btn.label.set_color("white")
            for sp in axb.spines.values():
                sp.set_color("white" if active else "#444")
                sp.set_linewidth(1.6 if active else 0.6)
            btn.on_clicked(lambda e, idx=i: self._select_view(idx))
            self._view_btns.append((btn, color))

        # ── prev / next / page-jump (left side) ───────────────────────
        nav_x = 0.04
        btn_w = 0.062
        tb_w  = 0.105   # page jump text box — wide enough for 4-digit numbers
        gap   = 0.006

        ax_prev = self.fig.add_axes([nav_x, rA_y, btn_w, rA_h])
        self._btn_prev = Button(ax_prev, "◀  Prev",
                                color="#263238", hovercolor="#37474F")
        self._btn_prev.label.set_fontsize(9)
        self._btn_prev.label.set_color("white")
        self._btn_prev.on_clicked(self._prev)

        # Page-jump TextBox
        tb_x = nav_x + btn_w + gap
        ax_tb = self.fig.add_axes([tb_x, rA_y + rA_h * 0.05, tb_w, rA_h * 0.90])
        ax_tb.set_facecolor("#0d1117")
        for sp in ax_tb.spines.values():
            sp.set_color("#555")
        self._page_tb = TextBox(ax_tb, "", initial="1",
                                color="#0d1117", hovercolor="#1a2a3a",
                                label_pad=0.02)
        self._page_tb.label.set_color("white")
        self._page_tb.text_disp.set_color("white")
        self._page_tb.on_submit(self._on_page_jump)

        # Page total label
        ax_pg = self.fig.add_axes([tb_x + tb_w + gap, rA_y, 0.055, rA_h])
        ax_pg.set_xticks([]); ax_pg.set_yticks([])
        ax_pg.set_facecolor("#0d1117")
        for sp in ax_pg.spines.values():
            sp.set_color("#333")
        self._page_total_txt = ax_pg.text(
            0.5, 0.5, "", ha="center", va="center",
            color="#aaa", fontsize=8, transform=ax_pg.transAxes)

        next_x = tb_x + tb_w + gap + 0.055 + gap
        ax_next = self.fig.add_axes([next_x, rA_y, btn_w, rA_h])
        self._btn_next = Button(ax_next, "Next  ▶",
                                color="#263238", hovercolor="#37474F")
        self._btn_next.label.set_fontsize(9)
        self._btn_next.label.set_color("white")
        self._btn_next.on_clicked(self._next)

        # ── per-page slider (shifted right to clear the Next button) ──
        pp_x = next_x + btn_w + 0.12
        ax_pp = self.fig.add_axes([pp_x, rA_y + rA_h * 0.1, 0.11, rA_h * 0.7])
        ax_pp.set_facecolor("#16213e")
        self._slider_pp = Slider(ax_pp, "Per page", 1, 16,
                                 valinit=st.n_per_page, valstep=1,
                                 color="#4CAF50")
        self._slider_pp.label.set_color("white")
        self._slider_pp.valtext.set_color("white")
        self._slider_pp.on_changed(self._on_per_page)

        # ── 📊 Stats button — far right, clearly separated ────────────
        ax_st = self.fig.add_axes([0.87, rA_y, 0.09, rA_h])
        self._btn_stats = Button(ax_st, "[Stats]",
                                 color="#37474F", hovercolor="#546E7A")
        self._btn_stats.label.set_fontsize(9)
        self._btn_stats.label.set_color("white")
        self._btn_stats.on_clicked(self._show_stats)

        # ── 4MOST smoothing slider ────────────────────────────────────
        ax_sl = self.fig.add_axes([0.055, rB_y, 0.22, rB_h])
        ax_sl.set_facecolor("#16213e")
        self._slider_sigma = Slider(
            ax_sl, "4MOST σ", 0, 20,
            valinit=st.sigma, valstep=0.5,
            color="#FF9800")
        self._slider_sigma.label.set_color("white")
        self._slider_sigma.valtext.set_color("white")
        self._slider_sigma.on_changed(self._on_sigma)

        # ── DESI smoothing slider ─────────────────────────────────────
        ax_sd = self.fig.add_axes([0.35, rB_y, 0.22, rB_h])
        ax_sd.set_facecolor("#16213e")
        self._slider_sigma_desi = Slider(
            ax_sd, "DESI σ", 0, 20,
            valinit=st.sigma_desi, valstep=0.5,
            color=C_DESI)
        self._slider_sigma_desi.label.set_color("white")
        self._slider_sigma_desi.valtext.set_color("white")
        self._slider_sigma_desi.on_changed(self._on_sigma_desi)

        # ── sky scale slider (÷1 … ÷10) ──────────────────────────────
        ax_sky = self.fig.add_axes([0.64, rB_y, 0.18, rB_h])
        ax_sky.set_facecolor("#16213e")
        self._slider_sky = Slider(
            ax_sky, "Sky ÷", 1, 30,
            valinit=st.sky_scale, valstep=0.5,
            color="#44cc55")
        self._slider_sky.label.set_color("white")
        self._slider_sky.valtext.set_color("white")
        self._slider_sky.on_changed(self._on_sky_scale)

        # ── component toggle buttons (replace static legend) ──────────
        rC_y = 0.030; rC_h = 0.030
        _toggle_defs = [
            ("raw",    C_RAW,    "show_raw"),
            ("smooth", C_SMOOTH, "show_smooth"),
            ("sky",    C_SKY,    "show_sky"),
            ("err",    C_ERR,    "show_err"),
            ("desi",   C_DESI,   "show_desi"),
            ("alt",    C_ALT,    "show_alt"),
            ("lines",  C_EM,     "show_lines"),
        ]
        self._toggle_btns = {}
        tbtn_w = 0.072
        tbtn_gap = 0.010
        tbtn_x0 = 0.04
        for i, (lbl, on_col, attr) in enumerate(_toggle_defs):
            bx = tbtn_x0 + i * (tbtn_w + tbtn_gap)
            is_on = getattr(st, attr, True)
            fc = on_col if is_on else "#2a2a2a"
            axb = self.fig.add_axes([bx, rC_y, tbtn_w, rC_h])
            btn = Button(axb, lbl, color=fc,
                         hovercolor=on_col if is_on else "#3a3a3a")
            btn.label.set_fontsize(8)
            btn.label.set_color("white")
            btn.on_clicked(lambda e, a=attr, c=on_col: self._on_toggle(a, c))
            self._toggle_btns[attr] = (btn, on_col)

        # ── keyboard hint ─────────────────────────────────────────────
        self.fig.text(0.86, 0.044,
                      "← →  pages   1-9  tab   s  smooth   q  quit",
                      color="#555", fontsize=6.5, ha="left", va="bottom")

        # ── active quality-filter indicator (visible in any mode) ─────
        self._qfilter_txt = self.fig.text(0.04, 0.185, "", fontsize=8,
                                          ha="left", va="bottom")

        # ── tagging sidebar (single-spectrum mode only) ───────────────
        if self._single:
            self._build_tag_sidebar()


    def _build_tag_sidebar(self):
        """Right-hand quality/issue/redshift panel; single-spectrum mode only."""
        sx, sw = 0.815, 0.165
        self.fig.text(sx, 0.852, "Tag spectrum", color="white",
                      fontsize=10, fontweight="bold", va="top")

        # ── quality (each row: tag button + "only" view filter) ────────
        self._qual_btns = {}
        self._qfilter_btns = {}
        self.fig.text(sx, 0.828, "Quality        (only = filter view)",
                      color="#aaa", fontsize=7.5, va="top")
        qy, qh, qg = 0.802, 0.028, 0.005
        for i, q in enumerate(QUALITIES):
            y = qy - i * (qh + qg)
            axb = self.fig.add_axes([sx, y, sw * 0.70, qh])
            btn = Button(axb, q, color="#2a2a2a", hovercolor=QUAL_COLORS[q])
            btn.label.set_fontsize(8); btn.label.set_color("white")
            btn.on_clicked(lambda e, qq=q: self._on_quality(qq))
            self._qual_btns[q] = btn

            axf = self.fig.add_axes([sx + sw * 0.72, y, sw * 0.28, qh])
            fbtn = Button(axf, "only", color="#2a2a2a", hovercolor="#455A64")
            fbtn.label.set_fontsize(7); fbtn.label.set_color("white")
            fbtn.on_clicked(lambda e, qq=q: self._on_qual_filter(qq))
            self._qfilter_btns[q] = fbtn

        # ── issues ────────────────────────────────────────────────────
        self._issue_btns = {}
        iy0 = qy - len(QUALITIES) * (qh + qg) - 0.012
        self.fig.text(sx, iy0, "Issues", color="#aaa", fontsize=8, va="top")
        iy, ih, ig = iy0 - 0.026, 0.022, 0.005
        for i, (key, label) in enumerate(ISSUE_DEFS):
            axb = self.fig.add_axes([sx, iy - i * (ih + ig), sw, ih])
            btn = Button(axb, label, color="#2a2a2a", hovercolor="#455A64")
            btn.label.set_fontsize(7.5); btn.label.set_color("white")
            btn.on_clicked(lambda e, k=key: self._on_issue(k))
            self._issue_btns[key] = btn

        # ── redshift override ─────────────────────────────────────────
        ry = iy - len(ISSUE_DEFS) * (ih + ig) - 0.012
        self._z_hdr = self.fig.text(sx, ry, "Redshift", color="#aaa",
                                    fontsize=8, va="top")
        ax_zt = self.fig.add_axes([sx, ry - 0.036, sw, 0.028])
        self._z_toggle_btn = Button(ax_zt, "custom z: off",
                                    color="#2a2a2a", hovercolor="#455A64")
        self._z_toggle_btn.label.set_fontsize(8)
        self._z_toggle_btn.label.set_color("white")
        self._z_toggle_btn.on_clicked(self._on_toggle_customz)

        ax_zs = self.fig.add_axes([sx + 0.028, ry - 0.072, sw - 0.028, 0.022])
        ax_zs.set_facecolor("#16213e")
        self._z_slider = Slider(ax_zs, "z", 0.0, 6.5, valinit=0.0,
                                valstep=0.001, color=C_EM)
        self._z_slider.label.set_color("white"); self._z_slider.label.set_fontsize(8)
        self._z_slider.valtext.set_color("white"); self._z_slider.valtext.set_fontsize(7)
        self._z_slider.on_changed(self._on_z_slider)

        self.fig.text(sx, ry - 0.096, "set z:", color="#aaa", fontsize=7.5,
                      va="top")
        ax_zb = self.fig.add_axes([sx + 0.05, ry - 0.104, sw - 0.05, 0.026])
        ax_zb.set_facecolor("#0d1117")
        for sp in ax_zb.spines.values():
            sp.set_color("#555")
        self._z_box = TextBox(ax_zb, "", initial="",
                              color="#0d1117", hovercolor="#1a2a3a")
        self._z_box.text_disp.set_color("white")
        self._z_box.on_submit(self._on_z_box)

        # ── save ──────────────────────────────────────────────────────
        save_y = ry - 0.150
        axs = self.fig.add_axes([sx, save_y, sw, 0.034])
        self._btn_save_tags = Button(axs, "Save tags",
                                     color="#1565C0", hovercolor="#1E88E5")
        self._btn_save_tags.label.set_fontsize(9)
        self._btn_save_tags.label.set_color("white")
        self._btn_save_tags.on_clicked(self._on_save_tags)
        self._tag_status = self.fig.text(sx, save_y - 0.016, "", color="#888",
                                         fontsize=7, va="top")


    # ── drawing ────────────────────────────────────────────────────────

    def _draw_spectra(self):
        st   = self.state
        rows = st.page_rows()
        col  = st.view_color
        n    = len(rows)

        if len(self.spec_axes) != st.n_per_page:
            self._build_spectra_ui()

        # scale the per-panel title font with the grid column count so the
        # filename/stats lines don't overflow into neighbouring panels
        _nc = _grid_dims(st.n_per_page)[1]
        title_fs = max(5.5, min(13.0, 26.0 / _nc))

        self._page_total_txt.set_text(f"/ {st.n_pages()}")
        if hasattr(self, "_qfilter_txt"):
            if st.qual_filter:
                self._qfilter_txt.set_text(
                    f"● showing only quality = {st.qual_filter}"
                    f"  ({len(st._effective_df())})")
                self._qfilter_txt.set_color(QUAL_COLORS.get(st.qual_filter, "#aaa"))
            else:
                self._qfilter_txt.set_text("")

        for i, ax in enumerate(self.spec_axes):
            ax.set_visible(i < n)
            if i >= n:
                continue
            row = rows.iloc[i]
            fp  = str(row.get("filepath", ""))
            cat = str(row.get("category", ""))
            fn  = row.get("filename", "")
            # default (blue) trace = L1 combination (LR4); nothing if absent
            w, f, e, sky, is_l1 = st.get_default(cat, fn, fp)
            dw, df_ = st.get_desi(fn)
            # "alt" overlay (orange) = my combination (LR)
            aw, af_ = st.get_overlay(fp) if is_l1 else (None, None)
            empty_msg = "No LR4 spectrum" if not is_l1 else "Read error"
            # redshift override applies only in single-spectrum mode
            z_over = None
            if self._single and st.use_custom_z:
                z_over = st.get_tag(str(row.get("filename", ""))).get("z")
            draw_one_spectrum(ax, w, f, e, sky, row, st.sigma, col,
                              sky_scale=st.sky_scale,
                              show_raw=st.show_raw,
                              show_smooth=st.show_smooth,
                              show_sky=st.show_sky,
                              show_err=st.show_err,
                              desi_wave=dw, desi_flux=df_,
                              show_desi=st.show_desi,
                              desi_sigma=st.sigma_desi,
                              show_lines=st.show_lines,
                              z_override=z_over,
                              alt_wave=aw, alt_flux=af_,
                              show_alt=st.show_alt,
                              title_fontsize=title_fs,
                              empty_msg=empty_msg)

        # single-spectrum mode: track current file + sync tag buttons
        if getattr(self, "_single", False):
            if n:
                r0 = rows.iloc[0]
                self._current_fn = str(r0.get("filename", ""))
                try:
                    self._current_zbest = float(r0.get("XPCA_Z"))
                except (TypeError, ValueError):
                    self._current_zbest = float("nan")
            else:
                self._current_fn = ""
                self._current_zbest = float("nan")
            self._refresh_tag_buttons()

        self.fig.canvas.draw_idle()

    def _draw_stats(self):
        draw_stats_page(self.fig, self.state.df)
        # Back button
        ax_back = self.fig.add_axes([0.88, 0.005, 0.10, 0.038])
        btn = Button(ax_back, "← Spectra",
                     color="#37474F", hovercolor="#546E7A")
        btn.label.set_fontsize(8); btn.label.set_color("white")
        btn.on_clicked(self._back_to_spectra)
        self._stats_back_btn = btn

    # ── event handlers ──────────────────────────────────────────────────

    def _select_view(self, idx):
        st = self.state
        if not (0 <= idx < len(st.views)):
            return
        st.view = idx
        st.page = 0
        st.tab  = "spectra"
        # restyle tab buttons to reflect the active selection
        for i, (btn, color) in enumerate(self._view_btns):
            active = (i == idx)
            fc = color if active else "#2a2a2a"
            btn.ax.set_facecolor(fc)
            btn.color = fc
            for sp in btn.ax.spines.values():
                sp.set_color("white" if active else "#444")
                sp.set_linewidth(1.6 if active else 0.6)
        self._sync_page_tb()
        self._draw_spectra()

    def _prev(self, *_):
        if self.state.tab != "spectra":
            return
        self.state.page = max(0, self.state.page - 1)
        self._sync_page_tb()
        self._draw_spectra()

    def _next(self, *_):
        if self.state.tab != "spectra":
            return
        self.state.page = min(self.state.n_pages() - 1, self.state.page + 1)
        self._sync_page_tb()
        self._draw_spectra()

    def _on_page_jump(self, text):
        """Called when user presses Enter in the page TextBox."""
        try:
            pg = int(text.strip()) - 1   # convert from 1-indexed
            pg = max(0, min(pg, self.state.n_pages() - 1))
        except (ValueError, AttributeError):
            self._sync_page_tb()
            return
        self.state.page = pg
        self._draw_spectra()

    def _sync_page_tb(self):
        """Update the TextBox to show the current page number."""
        try:
            self._page_tb.set_val(str(self.state.page + 1))
        except Exception:
            pass

    def _on_sigma(self, val):
        self.state.sigma = float(val)
        if self.state.tab == "spectra":
            self._draw_spectra()

    def _on_sigma_desi(self, val):
        self.state.sigma_desi = float(val)
        if self.state.tab == "spectra":
            self._draw_spectra()

    def _on_sky_scale(self, val):
        self.state.sky_scale = float(val)
        if self.state.tab == "spectra":
            self._draw_spectra()

    def _on_toggle(self, attr, on_color):
        """Toggle a spectral component on/off and recolor its button."""
        new_val = not getattr(self.state, attr, True)
        setattr(self.state, attr, new_val)
        if attr in self._toggle_btns:
            btn, col = self._toggle_btns[attr]
            fc = col if new_val else "#2a2a2a"
            hc = col if new_val else "#3a3a3a"
            btn.ax.set_facecolor(fc)
            btn.color = fc
            btn.hovercolor = hc
        if self.state.tab == "spectra":
            self._draw_spectra()

    # ── tagging ─────────────────────────────────────────────────────────

    def _on_quality(self, q):
        fn = getattr(self, "_current_fn", "")
        if not fn:
            return
        self.state.set_quality(fn, q)
        self._refresh_tag_buttons()
        self.fig.canvas.draw_idle()

    def _on_issue(self, key):
        fn = getattr(self, "_current_fn", "")
        if not fn:
            return
        self.state.toggle_issue(fn, key)
        self._refresh_tag_buttons()
        self.fig.canvas.draw_idle()

    def _on_qual_filter(self, q):
        """Restrict the current view to spectra tagged with quality q (toggle)."""
        st = self.state
        st.qual_filter = None if st.qual_filter == q else q
        st.page = 0
        self._sync_page_tb()
        self._draw_spectra()

    def _on_save_tags(self, *_):
        n = save_tags(self.state.tags)
        self.state.tags_dirty = False
        if hasattr(self, "_tag_status"):
            self._tag_status.set_text(f"saved {n} → {os.path.basename(TAGS_CSV)}")
        self.fig.canvas.draw_idle()

    # ── redshift override ────────────────────────────────────────────────

    def _current_eff_z(self):
        """Effective redshift for the shown spectrum (custom if set, else best)."""
        tag = self.state.get_tag(getattr(self, "_current_fn", ""))
        if self.state.use_custom_z and tag.get("z") is not None:
            return tag["z"]
        return getattr(self, "_current_zbest", float("nan"))

    def _on_toggle_customz(self, *_):
        fn = getattr(self, "_current_fn", "")
        if not fn:
            return
        self.state.use_custom_z = not self.state.use_custom_z
        # seed the custom z from the best fit the first time it is enabled
        if self.state.use_custom_z and self.state.get_tag(fn).get("z") is None:
            zb = getattr(self, "_current_zbest", float("nan"))
            self.state.set_redshift(fn, 0.0 if not np.isfinite(zb) else zb)
        self._refresh_tag_buttons()
        self._draw_spectra()

    def _on_z_slider(self, val):
        if getattr(self, "_z_updating", False) or not self.state.use_custom_z:
            return
        fn = getattr(self, "_current_fn", "")
        if not fn:
            return
        self.state.set_redshift(fn, float(val))
        self._sync_z_box()
        self._draw_spectra()

    def _on_z_box(self, text):
        if getattr(self, "_z_updating", False):
            return
        text = text.strip()
        if not text:
            return
        try:
            z = max(0.0, min(6.5, float(text)))
        except ValueError:
            return
        fn = getattr(self, "_current_fn", "")
        if not fn:
            return
        self.state.use_custom_z = True
        self.state.set_redshift(fn, z)
        self._refresh_tag_buttons()
        self._draw_spectra()

    def _sync_z_box(self):
        if hasattr(self, "_z_box"):
            self._z_updating = True
            self._z_box.set_val(f"{self._current_eff_z():.4f}"
                                if np.isfinite(self._current_eff_z()) else "")
            self._z_updating = False

    def _refresh_tag_buttons(self):
        """Sync quality/issue/redshift widgets to the current spectrum's tags."""
        if not getattr(self, "_single", False) or not hasattr(self, "_qual_btns"):
            return
        tag = self.state.get_tag(getattr(self, "_current_fn", ""))
        for q, btn in self._qual_btns.items():
            on = (tag["quality"] == q)
            fc = QUAL_COLORS[q] if on else "#2a2a2a"
            btn.ax.set_facecolor(fc); btn.color = fc
        for q, btn in self._qfilter_btns.items():
            on = (self.state.qual_filter == q)
            fc = QUAL_COLORS[q] if on else "#2a2a2a"
            btn.ax.set_facecolor(fc); btn.color = fc
        for key, btn in self._issue_btns.items():
            on = (key in tag["issues"])
            fc = "#455A64" if on else "#2a2a2a"
            btn.ax.set_facecolor(fc); btn.color = fc

        # redshift widgets
        if hasattr(self, "_z_toggle_btn"):
            on = self.state.use_custom_z
            self._z_toggle_btn.label.set_text(f"custom z: {'on' if on else 'off'}")
            fc = C_EM if on else "#2a2a2a"
            self._z_toggle_btn.ax.set_facecolor(fc); self._z_toggle_btn.color = fc
            zeff = self._current_eff_z()
            self._z_updating = True
            if np.isfinite(zeff):
                self._z_slider.set_val(max(0.0, min(6.5, zeff)))
            self._z_updating = False
            self._sync_z_box()
            zb = getattr(self, "_current_zbest", float("nan"))
            self._z_hdr.set_text(f"Redshift  (best={zb:.3f})"
                                 if np.isfinite(zb) else "Redshift")

        if hasattr(self, "_tag_status"):
            dirty = "  •unsaved" if self.state.tags_dirty else ""
            zc = tag.get("z")
            zstr = f"  z*={zc:.3f}" if zc is not None else ""
            self._tag_status.set_text(
                (f"{tag['quality'] or 'untagged'}"
                 f"  [{len(tag['issues'])} issue(s)]{zstr}{dirty}"))

    def _on_per_page(self, val):
        st = self.state
        first_idx = st.page * st.n_per_page   # global index of first shown spectrum
        st.n_per_page = int(val)
        st.page = first_idx // st.n_per_page  # keep that spectrum on-screen
        self._sync_page_tb()
        if st.tab == "spectra":
            self._draw_spectra()

    def _show_stats(self, *_):
        self.state.tab = "stats"
        self._draw_stats()

    def _back_to_spectra(self, *_):
        self.state.tab = "spectra"
        self._build_spectra_ui()
        self._sync_page_tb()
        self._draw_spectra()

    def _on_key(self, event):
        # Suppress all shortcuts while any TextBox is being typed in
        for tb in (getattr(self, "_page_tb", None), getattr(self, "_z_box", None)):
            if tb is not None and getattr(tb, "capturekeystrokes", False):
                return
        k = event.key
        if k in ("right", "n"):
            self._next()
        elif k in ("left", "p"):
            self._prev()
        elif k == "s":
            new_s = (self.state.sigma + 2) % 22
            self._slider_sigma.set_val(new_s)
        elif k in ("q", "escape"):
            plt.close("all")
        elif k == "0":
            if self.state.tab == "stats":
                self._back_to_spectra()
            else:
                self._show_stats()
        elif k in ("1", "2", "3", "4", "5", "6", "7", "8", "9"):
            idx = int(k) - 1
            if idx < len(self.state.views):
                self._select_view(idx)

    def show(self):
        if _interactive:
            plt.show()
        else:
            print("\nNon-interactive mode: saving PNG previews …")
            for idx, (key, label, _c, _m) in enumerate(self.state.views):
                self._select_view(idx)
                out = os.path.join(BASE, f"preview_{key}.png")
                self.fig.savefig(out, dpi=120, bbox_inches="tight",
                                 facecolor="#1a1a2e")
                print(f"  Saved {out}")
            self.state.tab = "stats"
            self._draw_stats()
            out = os.path.join(BASE, "preview_stats.png")
            self.fig.savefig(out, dpi=120, bbox_inches="tight",
                             facecolor="#1a1a2e")
            print(f"  Saved {out}")
            print("Done.")


# ── entry point ──────────────────────────────────────────────────────────

def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else CSV_DEFAULT

    if not os.path.exists(csv_path):
        print(f"ERROR: {csv_path!r} not found.")
        print("Run   python collect_headers.py   first.")
        sys.exit(1)

    print(f"Loading {csv_path} …")
    df = pd.read_csv(csv_path, low_memory=False)

    csv_dir = os.path.dirname(os.path.abspath(csv_path))
    _walk_index = {}          # filename → full path, built once on first miss
    _index_built = [False]

    def _build_walk_index():
        for cat in ["LR", "HR", "LRD", "HIZ"]:
            d = os.path.join(csv_dir, cat)
            if not os.path.isdir(d):
                continue
            for root, _dirs, files in os.walk(d):
                for fn in files:
                    _walk_index.setdefault(fn, os.path.join(root, fn))
        _index_built[0] = True

    def _resolve(p):
        p = str(p)
        if os.path.isabs(p) and os.path.exists(p):
            return p
        rel = os.path.join(csv_dir, p)
        if os.path.exists(rel):
            return rel
        if not _index_built[0]:
            _build_walk_index()
        return _walk_index.get(os.path.basename(p), rel)

    df["filepath"] = df["filepath"].apply(_resolve)

    for col in ["RA", "DEC", "MEDSNR", "QZC_Z", "QSC_CLS",
                "CAT_MAG", "SEP_ARCSEC"]:
        if col not in df.columns:
            df[col] = np.nan

    cats    = dict(df["category"].value_counts())
    n_z     = int(pd.to_numeric(df["QZC_Z"], errors="coerce").notna().sum())
    n_match = int((pd.to_numeric(df["SEP_ARCSEC"], errors="coerce") < 3).sum())
    print(f"  {len(df):,} spectra  |  categories: {cats}")
    print(f"  {n_match:,} catalogue matches  |  {n_z:,} redshifts measured")

    # Auxiliary data (DESI, L1 combination, S16 catalogue, xpca) is resolved
    # relative to the CSV's directory, so a copy of viewer.py elsewhere still
    # finds it when pointed at the data CSV. Defaults to BASE when identical.
    import glob as _glob
    data_dir = csv_dir

    desi_index = build_desi_index(
        match_fits=os.path.join(data_dir, "SPV_DESI_match.fits"),
        desi_dir=os.path.join(data_dir, "DESI_spectra"))
    df["has_desi"] = df["filename"].apply(_coord_key).map(
        lambda k: k in desi_index)

    l1_index = build_l1_index(base=data_dir)

    _s16 = sorted(_glob.glob(os.path.join(data_dir,
                                          "S16_*target_catalogue*.fits*")))
    ss = assign_subsurvey(df,
                          cat_path=_s16[0] if _s16 else S16_CAT,
                          cache_path=os.path.join(data_dir,
                                                  ".subsurvey_cache.csv"))
    if ss is not None:
        df["subsurvey"] = ss

    xpca = load_xpca_redshifts(base=data_dir)
    if xpca is not None:
        df = df.join(xpca, on="filename")

    print()
    print("Controls:")
    print("  ← →     navigate pages     1-9      select tab (arm/subsurvey/DESI)")
    print("  TextBox  type page number  s        cycle smoothing")
    print("  0        toggle stats      q/Esc    quit")
    print("  Per page = 1  →  quality/issue tagging sidebar (Save tags to disk)")
    print()

    viewer = FourMostViewer(df, desi_index=desi_index, l1_index=l1_index)
    viewer.show()


if __name__ == "__main__":
    main()
