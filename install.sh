#!/usr/bin/env bash
#
# SpecView4M installer — fetches the spectra, DESI data and catalogues from
# Google Drive and unpacks them next to viewer.py so the tabs/overlays work.
#
# Usage:
#   ./install.sh          # default data: L1 combination (LR4, HIZ4) + DESI
#   ./install.sh --full   # also the LR overlay, subsurvey catalogue, xpca, conditions
#
# Requires: python, tar, network access. Uses `gdown` to pull from Drive.
#
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FULL=0
[ "${1:-}" = "--full" ] && FULL=1

echo "==> Installing Python dependencies"
python -m pip install -q -r requirements.txt "gdown>=5"

# fetch <file_id> <output>   — download a single Drive file (skip if present)
fetch () {
  if [ -e "$2" ]; then echo "  have $2"; return; fi
  echo "  downloading $2"
  gdown "https://drive.google.com/uc?id=$1" -O "$2"
}
# fetch_tar <file_id> <archive> — download + extract whatever dir it contains
# (marker-based skip, so it does not assume the internal directory name)
fetch_tar () {
  local marker=".installed_$2"
  if [ -f "$marker" ]; then echo "  have $2 (unpacked)"; return; fi
  fetch "$1" "$2"
  echo "  unpacking $2"
  tar xzf "$2" && rm -f "$2" && touch "$marker"
}

echo "==> Default data (spectra index + L1 combination + DESI)"

# REQUIRED index — SPV_objects.fits (same columns as spectra_headers.csv, but a
# compact FITS; viewer.py reads either). Auto-detected when you run the viewer.
fetch     1HmYz19x-9KuF7nNaPxx1ktPs0JnQMpIh SPV_objects.fits

# L1 combination = the default (blue) trace
fetch_tar 1MxGJEEOpSXIs868UKPbX-7XxOvI-DdxJ SPV_LR_L1ver.tar.gz               # -> LR4/
fetch_tar 1e2LwXTITLl8TIi3VDITBS5vNrVNv-n8T SPV_LRD_L1ver.tar.gz              # -> LRD4/
fetch_tar 1RgLW3W5rTel7YPNSSlc3R_rNuMDLrjQw SPV_HIZ_L1ver.tar.gz              # -> HIZ4/
fetch_tar 1d3TUeT0ava5WCNha-8zioIASv4Dc0RR9 SPV_HR_L1ver.tar.gz               # -> HR4/

# DESI overlay + tab
fetch_tar 13F7mOH-ZGY4LLRv_LcQ0UtFFtggt5_B5 DESI_spectra.tar.gz               # -> DESI_spectra/
fetch     1HS7pEP78L7uNF8DHPPCcCOkVaIFDX97d SPV_DESI_match.fits

if [ "$FULL" -eq 1 ]; then
  echo "==> Optional data (my LR/HIZ overlays, subsurvey catalogue, redshifts, conditions)"
  fetch_tar 1gZVBUJM9PYxiuGJFeVZ7iLWmvX66xhlw SPV_LR.tar.gz                   # -> LR/   (my combination, orange overlay)
  fetch_tar 18Gu0iPxqpbVhooeq-Rq7K168NyN8GwGn SPV_HIZ.tar.gz                  # -> HIZ/  (my combination, orange overlay)
  fetch     1O8kDRbDxml3MQgagfp5PCYg5PqJrdccb LR_xpca.fits
  fetch     1zD9i3hA1N4EozdKv0voGLY106Qyt7QcM HIZ_xpca.fits
  fetch     1ffwPkrFQEfop36Sd4rr9NE6YT8BfSuPt LRD_xpca.fits                   # Drive file is "LRD.xpca.fits"
  # (seeing/airmass/moon are now columns in SPV_objects.fits — no conditions.csv)

  # S16/SEP catalogue → subsurvey tabs (M_SED/W_SED/VAR_Z/VAR_G/VAR_L/W_HIZ)
  fetch     1BcOE8Uc1tCSnPNYm6e7lcpxf338KHc6c S16_20250429T0456Z_target_catalogue.fits.gz
fi

echo
echo "==> Done."
echo "Run:  python viewer.py       # auto-detects SPV_objects.fits"
