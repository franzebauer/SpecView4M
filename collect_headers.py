#!/usr/bin/env python3
"""
collect_headers.py  –  Extract metadata from all 4MOST L1 FITS spectra
and cross-match against a target catalogue.

Current catalogue: M_SED_cat.fits  (M_SED subsurvey, 230,346 targets)

Outputs:
  spectra_headers.csv   – one row per spectrum, headers + catalogue match
"""

import os, glob, struct, gzip, time, sys, argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--cats", nargs="+", default=["LR", "HR", "LRD", "HIZ"],
                help="Subdirectories to scan (default: LR HR LRD HIZ)")
ap.add_argument("--csv", default=None,
                help="Output CSV path (default: spectra_headers.csv)")
args = ap.parse_args()

# ── paths ──────────────────────────────────────────────────────────────────
BASE     = os.path.dirname(os.path.abspath(__file__))
CATS     = args.cats
CAT_FILE = os.path.join(BASE, "M_SED_cat.fits")
OUT_CSV  = args.csv or os.path.join(BASE, "spectra_headers.csv")

# Subsurvey priority (index = priority rank, lower = higher priority)
# M_SED_cat.fits contains only M_SED; list is kept for forward compatibility.
SUBSURVEY_PRIORITY = ["M_SED", "W_SED", "W_VARL", "W_VARZ", "W_VARG", "W_HIZ"]

# ── FITS header helpers ────────────────────────────────────────────────────

def _parse_cards(block):
    """Parse 80-char FITS cards from a 2880-byte block."""
    cards, end_found = {}, False
    for i in range(0, min(len(block), 2880), 80):
        card = block[i:i+80].decode("ascii", errors="replace")
        key  = card[:8].strip()
        rest = card[8:]
        if key == "END":
            end_found = True
            break
        if "=" in rest[:2]:
            val_comment = rest[2:]
            val = (val_comment.split("'")[1].strip()
                   if "'" in val_comment
                   else val_comment.split("/")[0].strip())
            if key:
                cards[key] = val
    return cards, end_found


def read_fits_headers(filepath):
    """
    Read all FITS header cards from a spectrum file.
    Returns a flat dict merging primary and first-extension headers.
    Uses astropy for fpack-compressed .fits.fz files.
    """
    if filepath.endswith(".fz"):
        import astropy.io.fits as afits
        merged = {}
        with afits.open(filepath) as hdul:
            for hdu in hdul[:2]:
                merged.update({k: v for k, v in hdu.header.items()
                                if k not in ("", "COMMENT", "HISTORY")})
        return merged

    opener = gzip.open if filepath.endswith(".gz") else open
    merged = {}
    in_primary = True
    with opener(filepath, "rb") as fh:
        while True:
            block = fh.read(2880)
            if not block or len(block) < 80:
                break
            cards, end_found = _parse_cards(block)
            merged.update(cards)
            if end_found:
                if in_primary:
                    in_primary = False     # move to extension header
                    continue
                else:
                    break                  # done after extension
    return merged


# ── M_SED_cat.fits catalogue reader ───────────────────────────────────────
#
# M_SED_cat.fits: NAXIS1=209 bytes/row, TFIELDS=35, 230,346 rows
# Key differences from S16:
#   RA/DEC are float32 (E) not float64 (D)
#   DATE_EARLIEST/DATE_LATEST are float32 (E) not float64 (D)
#   NAME=20A, SUBSURVEY=5A, RULESET=5A, MAG_TYPE=10A, CLASSIFICATION=3A
#   CAL_MAG_ID_BLUE/GREEN/RED=14A (full filter names, e.g. "DECAM_G_AB_PSF")
# Data starts at byte 11520 (4 × 2880-byte FITS header blocks)
CAT_DTYPE = np.dtype([
    ("NAME",              "S20"),
    ("RA",                ">f4"),   # float32 in this catalogue
    ("DEC",               ">f4"),   # float32 in this catalogue
    ("PMRA",              ">f4"),
    ("PMDEC",             ">f4"),
    ("PARALLAX",          ">f4"),
    ("EPOCH",             ">f4"),
    ("RESOLUTION",        ">i2"),
    ("SUBSURVEY",         "S5"),
    ("CADENCE",           ">i4"),
    ("TEMPLATE_REDSHIFT", ">f4"),
    ("RULESET",           "S5"),
    ("EXTENT_FLAG",       ">i2"),
    ("EXTENT_PARAMETER",  ">f4"),
    ("EXTENT_INDEX",      ">f4"),
    ("MAG_TYPE",          "S10"),
    ("MAG",               ">f4"),
    ("MAG_ERR",           ">f4"),
    ("DATE_EARLIEST",     ">f4"),   # float32 in this catalogue
    ("DATE_LATEST",       ">f4"),   # float32 in this catalogue
    ("REDDENING",         ">f4"),
    ("REDSHIFT_ESTIMATE", ">f4"),
    ("REDSHIFT_ERROR",    ">f4"),
    ("CLASSIFICATION",    "S3"),
    ("COMPLETENESS",      ">f4"),
    ("CAL_MAG_BLUE",      ">f4"),
    ("CAL_MAG_ERR_BLUE",  ">f4"),
    ("CAL_MAG_ID_BLUE",   "S14"),   # full filter name, e.g. "DECAM_G_AB_PSF"
    ("CAL_MAG_GREEN",     ">f4"),
    ("CAL_MAG_ERR_GREEN", ">f4"),
    ("CAL_MAG_ID_GREEN",  "S14"),
    ("CAL_MAG_RED",       ">f4"),
    ("CAL_MAG_ERR_RED",   ">f4"),
    ("CAL_MAG_ID_RED",    "S14"),
    ("TEMPLATE",          "S24"),
])
assert CAT_DTYPE.itemsize == 209, f"dtype size mismatch: {CAT_DTYPE.itemsize}"
CAT_NROWS  = 230_346
CAT_HDRBLK = 4   # 4 × 2880 = 11520 bytes before binary data


def read_catalogue(filepath):
    """Read S16 target catalogue → structured dict of numpy arrays."""
    print(f"Reading M_SED catalogue ({CAT_NROWS:,} rows, "
          f"{CAT_NROWS * CAT_DTYPE.itemsize / 1e6:.0f} MB) …", flush=True)
    t0 = time.time()
    with open(filepath, "rb") as fh:
        fh.read(CAT_HDRBLK * 2880)           # skip FITS header blocks
        raw = fh.read(CAT_NROWS * CAT_DTYPE.itemsize)
    data = np.frombuffer(raw, dtype=CAT_DTYPE, count=CAT_NROWS)
    print(f"  Done in {time.time()-t0:.1f}s. "
          f"RA [{data['RA'].min():.2f}°, {data['RA'].max():.2f}°]  "
          f"Dec [{data['DEC'].min():.2f}°, {data['DEC'].max():.2f}°]", flush=True)

    def _str(arr):
        return np.array([v.decode("ascii", "replace").strip() for v in arr])

    return {
        "NAME":           _str(data["NAME"]),
        "RA":             data["RA"].astype(np.float64),
        "DEC":            data["DEC"].astype(np.float64),
        "MAG":            data["MAG"].astype(np.float32),
        "MAG_ERR":        data["MAG_ERR"].astype(np.float32),
        "MAG_TYPE":       _str(data["MAG_TYPE"]),
        "CLASSIFICATION": _str(data["CLASSIFICATION"]),
        "RULESET":        _str(data["RULESET"]),
        "REDSHIFT_ESTIMATE": data["REDSHIFT_ESTIMATE"].astype(np.float32),
        "TEMPLATE":       _str(data["TEMPLATE"]),
        "SUBSURVEY":      _str(data["SUBSURVEY"]),
        "CAL_MAG_BLUE":   data["CAL_MAG_BLUE"].astype(np.float32),
        "CAL_MAG_ERR_BLUE": data["CAL_MAG_ERR_BLUE"].astype(np.float32),
        "CAL_MAG_ID_BLUE":  _str(data["CAL_MAG_ID_BLUE"]),
        "CAL_MAG_GREEN":  data["CAL_MAG_GREEN"].astype(np.float32),
        "CAL_MAG_ERR_GREEN": data["CAL_MAG_ERR_GREEN"].astype(np.float32),
        "CAL_MAG_ID_GREEN":  _str(data["CAL_MAG_ID_GREEN"]),
        "CAL_MAG_RED":    data["CAL_MAG_RED"].astype(np.float32),
        "CAL_MAG_ERR_RED": data["CAL_MAG_ERR_RED"].astype(np.float32),
        "CAL_MAG_ID_RED":  _str(data["CAL_MAG_ID_RED"]),
    }


# ── fast cross-matching ────────────────────────────────────────────────────

def crossmatch_radec(spec_ra, spec_dec, cat_ra, cat_dec, max_sep_arcsec=3.0):
    """
    Match (spec_ra, spec_dec) against (cat_ra, cat_dec) using sorted-Dec.
    Returns (match_idx, match_sep_arcsec).  match_idx == -1 means no match.
    match_idx is an index into cat_ra/cat_dec.
    """
    max_sep_deg = max_sep_arcsec / 3600.0
    n_spec = len(spec_ra)

    sort_idx   = np.argsort(cat_dec)
    sorted_dec = cat_dec[sort_idx]
    sorted_ra  = cat_ra[sort_idx]

    match_idx = np.full(n_spec, -1, dtype=np.int64)
    match_sep = np.full(n_spec, np.inf)

    for i in range(n_spec):
        sdec = spec_dec[i];  sra = spec_ra[i]
        if not (np.isfinite(sdec) and np.isfinite(sra)):
            continue

        lo = int(np.searchsorted(sorted_dec, sdec - max_sep_deg))
        hi = int(np.searchsorted(sorted_dec, sdec + max_sep_deg))
        if lo >= hi:
            continue

        cos_dec  = max(abs(np.cos(np.radians(sdec))), 1e-4)
        ra_range = max_sep_deg / cos_dec
        ra_diff  = np.abs(sorted_ra[lo:hi] - sra)
        ra_diff  = np.minimum(ra_diff, 360.0 - ra_diff)
        ra_mask  = ra_diff <= ra_range
        if not ra_mask.any():
            continue

        c_dec = np.radians(sorted_dec[lo:hi][ra_mask])
        c_ra  = np.radians(sorted_ra[lo:hi][ra_mask])
        s_dec = np.radians(sdec);  s_ra = np.radians(sra)
        cos_sep = np.clip(
            np.sin(s_dec)*np.sin(c_dec) +
            np.cos(s_dec)*np.cos(c_dec)*np.cos(c_ra - s_ra),
            -1.0, 1.0)
        sep_deg = np.degrees(np.arccos(cos_sep))

        best = int(np.argmin(sep_deg))
        if sep_deg[best] <= max_sep_deg:
            local_idx     = np.where(ra_mask)[0][best]
            match_idx[i]  = sort_idx[lo + local_idx]
            match_sep[i]  = sep_deg[best] * 3600.0

    return match_idx, match_sep


def crossmatch_with_priority(spec_ra, spec_dec, cat, max_sep_arcsec=3.0):
    """
    Cross-match spectra against the catalogue, respecting SUBSURVEY_PRIORITY.
    For each spectrum, the match from the highest-priority subsurvey wins.
    Returns (match_idx_in_cat, match_sep_arcsec).
    """
    n_spec = len(spec_ra)
    match_idx_global = np.full(n_spec, -1, dtype=np.int64)
    match_sep_global = np.full(n_spec, np.inf)

    for ss in SUBSURVEY_PRIORITY:
        # How many are still unmatched?
        unmatched = match_idx_global < 0
        n_unmatched = unmatched.sum()
        if n_unmatched == 0:
            break

        # Filter catalogue to this subsurvey
        ss_mask  = cat["SUBSURVEY"] == ss
        n_ss     = ss_mask.sum()
        if n_ss == 0:
            print(f"  [{ss}] not found in catalogue — skipping", flush=True)
            continue

        orig_idx = np.where(ss_mask)[0]   # positions in full cat arrays
        sub_ra   = cat["RA"][ss_mask]
        sub_dec  = cat["DEC"][ss_mask]

        print(f"  [{ss}] {n_ss:,} cat rows, "
              f"matching {n_unmatched:,} unmatched spectra …",
              end="", flush=True)
        t0 = time.time()

        um_ra  = spec_ra[unmatched]
        um_dec = spec_dec[unmatched]
        sub_match, sub_sep = crossmatch_radec(um_ra, um_dec, sub_ra, sub_dec,
                                              max_sep_arcsec)

        um_positions = np.where(unmatched)[0]
        newly = 0
        for j, gpos in enumerate(um_positions):
            if sub_match[j] >= 0:
                match_idx_global[gpos] = orig_idx[sub_match[j]]
                match_sep_global[gpos] = sub_sep[j]
                newly += 1

        print(f"  matched {newly:,} new  ({time.time()-t0:.1f}s)", flush=True)

    return match_idx_global, match_sep_global


# ── header keys to extract from spectra ───────────────────────────────────

PRIMARY_KEYS = [
    "RA", "DEC", "OBJECT", "EXPTIME", "TEXPTIME",
    "MJD-OBS", "DATE-OBS", "PROG_ID", "FLUXCAL",
]
EXT_KEYS = [
    "MEDSNR", "MINSNR", "MAXSNR", "SNR_4L1",
    "WAVELMIN", "WAVELMAX", "NELEM", "SPEC_RES",
    "QZC_Z", "QZC_ZERR", "QZC_ZPRB", "QZC_ZTPL",
    "QSC_CLS", "QSC_CPRB",
    "OBJ_UID", "FIB_ID", "PATH", "NCOMBINE",
    "EXT_OBJ", "HELIOCR",
]
NUMERIC_KEYS = {
    "RA", "DEC", "EXPTIME", "TEXPTIME", "MJD-OBS",
    "MEDSNR", "MINSNR", "MAXSNR", "SNR_4L1",
    "WAVELMIN", "WAVELMAX", "NELEM", "SPEC_RES",
    "QZC_Z", "QZC_ZERR", "QZC_ZPRB", "QSC_CPRB",
    "NCOMBINE", "HELIOCR",
}


def extract_row(filepath, category):
    try:
        cards = read_fits_headers(filepath)
    except Exception as e:
        return {"filepath": os.path.relpath(filepath, BASE),
                "category": category, "_error": str(e)}

    row = {"filepath": os.path.relpath(filepath, BASE),
           "category": category,
           "filename": os.path.basename(filepath)}
    for k in PRIMARY_KEYS + EXT_KEYS:
        v = cards.get(k, None)
        if k in NUMERIC_KEYS:
            try:
                row[k] = float(v) if v not in (None, "", "N/A") else np.nan
            except (ValueError, TypeError):
                row[k] = np.nan
        else:
            row[k] = v
    return row


# ── collection ────────────────────────────────────────────────────────────

def collect_headers():
    rows = []
    t0   = time.time()
    for cat in CATS:
        files = sorted(
            glob.glob(os.path.join(BASE, cat, "**", "*.fits"),    recursive=True) +
            glob.glob(os.path.join(BASE, cat, "**", "*.fits.fz"), recursive=True)
        )
        print(f"  {cat}: {len(files)} files", flush=True)
        n = len(files)
        for i, fp in enumerate(files):
            rows.append(extract_row(fp, cat))
            if (i + 1) % 500 == 0 or i == n - 1:
                elapsed = time.time() - t0
                rate    = (i + 1) / elapsed
                print(f"    [{cat}] {i+1}/{n}  ({rate:.0f} files/s)",
                      end="\r", flush=True)
        print(f"    [{cat}] done ({n} files)              ", flush=True)

    df = pd.DataFrame(rows)
    print(f"Collected {len(df):,} headers in {time.time()-t0:.1f}s")
    return df


# ── main ──────────────────────────────────────────────────────────────────

def main():
    # ── Step 1: collect headers (always full re-ingest) ────────────────
    print("Collecting FITS headers (full re-ingest) …")
    df = collect_headers()
    df.to_csv(OUT_CSV, index=False)
    print(f"Saved header-only CSV → {OUT_CSV}")

    if len(df) == 0:
        print("No files collected — exiting.")
        return

    # ── Step 2: load M_SED catalogue ───────────────────────────────────
    try:
        cat = read_catalogue(CAT_FILE)
    except Exception as e:
        print(f"WARNING: Could not read catalogue: {e}")
        print("Saving without crossmatch.")
        df.to_csv(OUT_CSV, index=False)
        return

    # ── Step 3: priority crossmatch ────────────────────────────────────
    print(f"\nCross-matching {len(df):,} spectra (priority order: "
          f"{' > '.join(SUBSURVEY_PRIORITY)}) …", flush=True)
    t0       = time.time()
    spec_ra  = pd.to_numeric(df["RA"],  errors="coerce").values
    spec_dec = pd.to_numeric(df["DEC"], errors="coerce").values

    match_idx, match_sep = crossmatch_with_priority(spec_ra, spec_dec, cat,
                                                    max_sep_arcsec=3.0)
    n_matched = (match_idx >= 0).sum()
    print(f"\n  Total matched {n_matched:,}/{len(df):,} spectra "
          f"in {time.time()-t0:.1f}s", flush=True)

    # ── Attach catalogue columns ───────────────────────────────────────
    df["SEP_ARCSEC"]        = match_sep
    df["CAT_NAME"]          = ""
    df["CAT_RA"]            = np.nan
    df["CAT_DEC"]           = np.nan
    df["CAT_MAG"]           = np.nan
    df["CAT_MAG_ERR"]       = np.nan
    df["CAT_MAG_TYPE"]      = ""
    df["CAT_CLASS"]         = ""
    df["CAT_RULESET"]       = ""
    df["CAT_Z_EST"]         = np.nan
    df["CAT_TEMPLATE"]      = ""
    df["CAT_SUBSURVEY"]     = ""
    df["CAL_MAG_BLUE"]      = np.nan
    df["CAL_MAG_ERR_BLUE"]  = np.nan
    df["CAL_MAG_ID_BLUE"]   = ""
    df["CAL_MAG_GREEN"]     = np.nan
    df["CAL_MAG_ERR_GREEN"] = np.nan
    df["CAL_MAG_ID_GREEN"]  = ""
    df["CAL_MAG_RED"]       = np.nan
    df["CAL_MAG_ERR_RED"]   = np.nan
    df["CAL_MAG_ID_RED"]    = ""

    good = match_idx >= 0
    for col, arr in [
        ("CAT_NAME",          cat["NAME"]),
        ("CAT_RA",            cat["RA"]),
        ("CAT_DEC",           cat["DEC"]),
        ("CAT_MAG",           cat["MAG"]),
        ("CAT_MAG_ERR",       cat["MAG_ERR"]),
        ("CAT_MAG_TYPE",      cat["MAG_TYPE"]),
        ("CAT_CLASS",         cat["CLASSIFICATION"]),
        ("CAT_RULESET",       cat["RULESET"]),
        ("CAT_Z_EST",         cat["REDSHIFT_ESTIMATE"]),
        ("CAT_TEMPLATE",      cat["TEMPLATE"]),
        ("CAT_SUBSURVEY",     cat["SUBSURVEY"]),
        ("CAL_MAG_BLUE",      cat["CAL_MAG_BLUE"]),
        ("CAL_MAG_ERR_BLUE",  cat["CAL_MAG_ERR_BLUE"]),
        ("CAL_MAG_ID_BLUE",   cat["CAL_MAG_ID_BLUE"]),
        ("CAL_MAG_GREEN",     cat["CAL_MAG_GREEN"]),
        ("CAL_MAG_ERR_GREEN", cat["CAL_MAG_ERR_GREEN"]),
        ("CAL_MAG_ID_GREEN",  cat["CAL_MAG_ID_GREEN"]),
        ("CAL_MAG_RED",       cat["CAL_MAG_RED"]),
        ("CAL_MAG_ERR_RED",   cat["CAL_MAG_ERR_RED"]),
        ("CAL_MAG_ID_RED",    cat["CAL_MAG_ID_RED"]),
    ]:
        df.loc[good, col] = arr[match_idx[good]]

    # ── Signed positional offsets (arcsec) ────────────────────────────
    spec_ra_v  = pd.to_numeric(df["RA"],      errors="coerce").values
    spec_dec_v = pd.to_numeric(df["DEC"],     errors="coerce").values
    cat_ra_v   = pd.to_numeric(df["CAT_RA"],  errors="coerce").values
    cat_dec_v  = pd.to_numeric(df["CAT_DEC"], errors="coerce").values
    cos_dec    = np.cos(np.radians(np.where(np.isfinite(spec_dec_v), spec_dec_v, 0.0)))
    df["DELTA_RA_AS"]  = (spec_ra_v  - cat_ra_v)  * cos_dec * 3600.0
    df["DELTA_DEC_AS"] = (spec_dec_v - cat_dec_v) * 3600.0

    n_half = int((pd.to_numeric(df["SEP_ARCSEC"], errors="coerce") > 0.5).sum())
    print(f"  Spectra with no match within 0.5\": {n_half:,}/{len(df):,}")

    # Show subsurvey breakdown
    print("\nMatch breakdown by subsurvey:")
    for ss in SUBSURVEY_PRIORITY:
        n = (df["CAT_SUBSURVEY"] == ss).sum()
        print(f"  {ss:8s}: {n:5,}")
    print(f"  {'unmatched':8s}: {(~good).sum():5,}")

    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved matched CSV → {OUT_CSV}")
    print()
    print(df[["filename","category","RA","DEC",
              "QZC_Z","QSC_CLS","CAT_MAG","CAT_SUBSURVEY",
              "SEP_ARCSEC","DELTA_RA_AS","DELTA_DEC_AS"]].head(10).to_string())


if __name__ == "__main__":
    main()
