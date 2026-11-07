#!/bin/bash
# Mirror what process_track.py needs to seasat (2026-10-27: GMTSAR processing runs there). Same relative paths under
# /zdata/datasets/IGML/code as under the repo, so REPO-relative paths in the script hold. The macOS MATLAB link in
# gmtsar_tools is left out (seasat's matlab is /usr/local/bin/matlab). Run from anywhere on the Mac after any change.
set -eu
R=$(cd "$(dirname "$0")/../.." && pwd)
DST=erapaportbruck@seasat.ucsd.edu:/zdata/datasets/IGML/code
cd "$R"
rsync -a --relative --exclude 'gmtsar_tools/matlab' \
    scratchpad/s1_gapfill/process_track.py \
    scratchpad/s1_gapfill/merge_seams.py \
    scratchpad/s1_gapfill/parity_align.py \
    scratchpad/s1_gapfill/rescue_test.py \
    scratchpad/s1_gapfill/prep_tropo_igml.sh \
    scratchpad/s1_gapfill/coseismic_events.json \
    scratchpad/s1_gapfill/gmtsar_tools \
    scratchpad/s1_gapfill/size_audit/needed_frames.json \
    scratchpad/s1_gapfill/size_audit/scene_timing.json \
    scratchpad/s1_inventory_20261015/s1_acquisitions.json \
    scratchpad/s1_inventory_20261015/extra_tracks/subswaths.json \
    "$DST/"
echo "synced to $DST"
