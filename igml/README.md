# IGML's additions to GMTSAR

This folder holds code written around GMTSAR for the IGML project. It is separate from GMTSAR's sources and build, so nothing here changes what `make` produces.

The fixes to GMTSAR itself are ordinary commits on this branch:

| Commit | Fix |
|---|---|
| 95a75c4 | A0: `stitch_tops` subtracts the last frame's overlap line offset once |
| 2286ce8 | I1, I2, I3, M1, M6: `conv` cells centred, subswaths placed exactly by `merge_swath`, no zero edges |
| 4476ff4 | A2: elevation antenna pattern removed for IPF < 002.43 scenes |
| f759576 | F6, A4: `assemble_tops` locates each scene's bursts with its own linesPerBurst and epoch |
| 3105ac0 | A1: every date cut between bursts at the super-master's lines (`IGML_A1`) |
| e94dbe8 | A3: a resampled date keeps its far-range columns (`IGML_A3`; not used in production) |
| 94beffb | A5: `preproc_batch_tops_esd.csh` collects the burst overlaps and applies a network's shifts (`IGML_ESD_COLLECT`, `IGML_ESD_SHIFTS`) |
| 90b6b68 | A6: spectral diversity's Doppler separation over the burst cycle, not a whole burst |
| 680971e | `filter.csh` reuses each reference date's decimated amplitude (`IGML_AMP_CACHE`; byte-identical, faster) |

## What is here

| Folder | Contents |
|---|---|
| `nesd/` | **Network spectral diversity**, standalone: measure the pairs, solve, write the shifts. Start here for NESD; see its README |
| `pipeline/` | `process_track.py`, the driver that runs a whole track from raw scenes to replaced stacks, with its helpers. NESD is also inside it (`nesd`, `_nesd_solve`, `_sd_pair`, `start_pass2`) |
| `patches/` | The scripts that patched the sources before the fork existed. The commits above supersede them |
| `audit/` | The checks and figures behind the reported numbers |
| `runs/` | The launch and chain scripts as they ran on seasat, including the gates between stages |
| `google_earth/` | The KMZ builders that show the geolocation error at the `conv` and `merge_swath` stages |

`pipeline/`, `audit/` and `runs/` are the record of what ran, not portable tools: they hold this site's paths (`/zdata`, `/palsar`, seasat) and need the data there. `nesd/nesd.py` has no such paths.

`pipeline/gmtsar_tools/` holds copies of three existing csh scripts that the driver puts on its PATH; they were not written for this project.
