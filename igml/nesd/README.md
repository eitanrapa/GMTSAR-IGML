# Network spectral diversity (NESD)

Geometric alignment leaves a small azimuth misregistration in each date, and TOPS turns that into a phase step at every burst boundary. GMTSAR's `preproc_batch_tops_esd.csh` estimates the shift of each date against the super-master alone, which is weak for dates years away, because the burst overlaps have decorrelated.

`nesd.py` measures each date against its next 4 dates instead, then solves one shift per date frame from all those pairs, with the super-master fixed at zero. The shift goes back into a second alignment.

## How to run it

Per subswath, on a stack with one super-master:

```csh
# 1. Pass 1: geometry only, keeping the burst overlaps
env IGML_ESD_COLLECT=1 preproc_batch_tops_esd.csh data.in dem.grd 2 1
nesd.py collect ALIGN_DIR ESD_DIR --master YYYYMMDD

# 2. Measure the pairs, solve, write ESD_DIR/shifts.txt
nesd.py run ESD_DIR --master YYYYMMDD --dates good_dates.txt

# 3. Pass 2: every date aligned again with its shift
env IGML_ESD_SHIFTS=ESD_DIR/shifts.txt preproc_batch_tops_esd.csh data.in dem.grd 2 1
```

It needs numpy, and this fork's `spectral_diversity` and `gmtsar_sharedir.csh` on PATH. `--dates` should list the dates that passed your alignment checks.

## The method

- **Measurement:** GMTSAR's `spectral_diversity` on a pair of dates' overlap images, taking the median of the double-difference phase over 2 pi times the Doppler separation.
- **Pair uncertainty:** the spread of those phases, with a floor of 0.02 rad.
- **Solve:** weighted least squares for x_b - x_a = y, with Huber reweighting (k = 1.345), separately for each super-master frame. The super-master is fixed at 0.
- **Output:** one constant per date frame, per subswath. A frame the network cannot reach, such as one a whole burst out of step with the super-master, keeps its estimate against the super-master and is listed as a fallback.

The network idea comes from Fattahi, Agram & Simons (2017), "A Network-Based Enhanced Spectral Diversity Approach for TOPS Time-Series Analysis," *IEEE TGRS* 55(2), 777-786, doi:10.1109/TGRS.2016.2614925. This code was not built from that paper's text, so the details above are its own choices and may differ from theirs.

## Results

| Test | Result |
|---|---|
| T166 F1 pair, burst step | -77.6 deg geometry only; -28.4 deg with GMTSAR's ESD against the super-master; -1.1 deg with NESD |
| T173, 431 dates, each subswath | 3,428 of 3,428 frame pairs used, 0 fallbacks; residual median 0.05 / 0.06 / 0.08 deg, p95 0.39 / 0.36 / 0.52 deg (F1 / F2 / F3) |
| T071, 428 dates, each subswath | 1,702 / 1,702 / 1,694 frame pairs all used, 1 fallback frame in F3; residual median 0.08 / 0.08 / 0.17 deg, p95 0.47 / 0.56 / 1.11 deg |

Residuals show that the network is self-consistent, not that the shifts are correct. The evidence for that is the burst steps disappearing.

## Phase closure

One shift per date cancels around any triangle, so NESD leaves closure alone. On 16 T173 triplets (`audit/nesd_closure.py`), wrapped closure rms was 5.1 -> 5.0 deg on mixed-satellite triplets and 5.5 -> 5.5 deg on S1A-only ones, without and with NESD.

This has two consequences:

- Closure cannot detect these steps or confirm the correction, because the uncorrected pairs close just as well.
- Pairs made from different solutions do not close with each other. Do not mix pre- and post-NESD pairs in one network, and do not re-solve existing dates when adding new ones. Keeping the existing shifts fixed while adding dates is not implemented yet.

## Relation to the pipeline

`nesd.py` is `process_track.py`'s `nesd()` and its helpers without the pipeline around them. The pipeline still runs its own copy.

Checked against the production runs on T071 and T173: `solve` on the cached pairs reproduces `shifts.txt` byte for byte on all three subswaths, and `measure` reproduces the `pairs.jsonl` records of sample pairs.
