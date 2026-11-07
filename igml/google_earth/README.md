# GMTSAR geolocation in Google Earth (T173, 2016-10-06 x 10-30)

Open with **Google Earth Pro** (desktop): File > Open. Radio lists behave best there.

| File | What it shows |
|---|---|
| `gmtsar_geolocation_t173_conv_stage.kmz` | Subswath grids straight from `conv`, before merging |
| `gmtsar_geolocation_t173_merged_stage.kmz` | Grids after `merge_swath`, plus two sites on the F1/F2 seam and the whole track |
| `gmtsar_geolocation_t173_merged_stage_sites_only.kmz` | The same without the whole-track layers, small enough to send |

**How to use:** each site has a "flicker these" list where one layer shows at a time. Double-click the folder to fly there, then click between Stock and Fixed. "go to" items zoom to a good spot.

**Suggested order:** Pinon Flat in the conv file, then in the merged file (the bright reflector moves off the surveyed point); Salton Sea west shore against the yellow line; Mexicali seam in the merged file (the canal steps at the white line).

## The numbers (stock GMTSAR's labels minus their content)

| Stage | F1 | F2 | F3 |
|---|---|---|---|
| conv | +16 samples, +3 lines | +16, +3 | +16, +3 |
| merge_swath | +16, +3 | +49, +1 | +29, +5 |
| Measured, merged, amplitude | +15.4, +3.0 | +48.6, +1.0 | +28.2, +4.9 |
| Ground, conv | 65-85 m | 65-85 m | 65-85 m |
| Ground, merged | 70-85 m | 180-190 m | 115-135 m |

Ground distance varies with terrain. conv's offset is by construction: a cell holds a filter centred on sample 32 m, line 6 r, and its node is declared at 32 m + 16, 6 r + 3. merge_swath adds whole-cell placement, and F2 lands one column too far.

**Pinon Flat reflectors from the survey (display_amp centroid):** stock conv 66 / 60 m, stock merged 116 / 110 m, fixed 2 / 21 m (D1 / D2).

## Caveats

- The Salton Sea has receded since 2016. Compare with the yellow Sentinel-2 line (2016-10-13) or Google Earth's historical imagery, not today's shore.
- Cross-correlation in a small window is noisy where bright point targets dominate (+5 to +17 samples for a true +16). The balloons quote the exact offsets; each window's measurements are in `measured_offsets.json`.
- Layers are interpolated between declared nodes through GMTSAR's `proj_ra2ll.csh` lookup tables, without its 100-m `blockmedian`. Both builds use the same `trans.dat`.
- "Stock" is the lab's archived product. "Fixed" also has the burst-alignment fixes, which do not move positions.

Rebuild: `/zdata/datasets/IGML/work/geoloc_kml_seasat.sh` on seasat, then `src/make_geolocation_kmz.py` and `src/preview_geolocation_kmz.py` on the Mac.
