"""Our patched copy of GMTSAR's make_s1a_tops.c compute_eap (2026-10-29; audit A2: the elevation-antenna-pattern phase
of IPF 002.36 scenes is never removed, -76 deg at the F1/F2 seam). GMTSAR removes it when the scene's manifest.safe
and the AUX_CAL xml are appended to the annotation (process_track.py does that for IPF < 002.43 dates); the routine
itself had, at master 3dd43d3:
  - the aux pattern picked by imageNumber: a dual-pol (1SDV) product numbers its VV annotations 4-6, which landed
    on an EW pattern; now by the annotation's swath (IW1-3), and the aux entry's swath / polarisation are checked
  - both linear interpolations (elevation angle in slant-range time; the pattern in elevation angle) with swapped
    weights
  - the second table's bound checked against the first table's length, and unbounded table searches
  - the range-time axis one sample early (tau0 + (i - 1) / fs)
src/make_s1a_tops.orig.c = the lab's source; src/make_s1a_tops.c is written here; build_eap.sh builds bin_eap/.
"""
from pathlib import Path

F = Path("/zdata/datasets/IGML/gmtsar_fix/src")
s = (F / "make_s1a_tops.orig.c").read_text()
reps = [
    ('	search_tree(xml_tree, "/product/adsHeader/imageNumber/", tmp_str, 1, 0, 1);\n'
     '	ii = (int)str2double(tmp_str);\n'
     '	mode = mode + (ii - 1) * 4;\n',
     '	/* IGML 2026-10-29: by swath, not imageNumber (a 1SDV product numbers its VV annotations 4-6) */\n'
     '	search_tree(xml_tree, "/product/adsHeader/swath/", tmp_str, 1, 0, 1);\n'
     '	if (strncmp(tmp_str, "IW", 2) != 0 || tmp_str[2] < \'1\' || tmp_str[2] > \'3\')\n'
     '		die("compute_eap: not an IW swath: ", tmp_str);\n'
     '	ii = tmp_str[2] - \'0\';\n'
     '	strncpy(eap_swath, tmp_str, 3);\n'
     '	eap_swath[3] = \'\\0\';\n'
     '	mode = mode + (ii - 1) * 4;\n'),
    ('	mode = mode + 24;\n',
     '	mode = mode + 24;\n'
     '	strncpy(eap_pol, tmp_str, 2);\n'
     '	eap_pol[2] = \'\\0\';\n'
     '	/* IGML: the aux entry must be this swath and polarisation */\n'
     '	search_tree(xml_tree, "/auxiliaryCalibration/calibrationParamsList/calibrationParams/swath/", tmp_str, 1, 3, mode);\n'
     '	if (strncmp(tmp_str, eap_swath, 3) != 0)\n'
     '		die("compute_eap: aux pattern is not for swath ", eap_swath);\n'
     '	search_tree(xml_tree, "/auxiliaryCalibration/calibrationParamsList/calibrationParams/polarisation/", tmp_str, 1, 3, mode);\n'
     '	if (strncmp(tmp_str, eap_pol, 2) != 0)\n'
     '		die("compute_eap: aux pattern is not for polarisation ", eap_pol);\n'
     '	fprintf(stderr, "(EAP %s %s, aux pattern %d)", eap_swath, eap_pol, mode);\n'),
    ('	double H[5], Torb, Worb, phi[5], real, imag;\n',
     '	double H[5], Torb, Worb, phi[5], real, imag;\n'
     '	char eap_swath[8], eap_pol[4];\n'),
    ('		tau[ii] = tau0 + (double)(ii - 1) / fs;\n',
     '		tau[ii] = tau0 + (double)ii / fs; /* IGML: sample ii is at tau0 + ii / fs */\n'),
    ('	while (tau_sub[jj] < tau0)\n		jj++;\n',
     '	while (jj < srtcount && tau_sub[jj] < tau0)\n		jj++;\n'),
    ('		while (tau_sub[jj] < tau[ii])\n			jj++;\n'
     '		if (jj < 1 || jj > srtcount) {\n',
     '		while (jj < srtcount && tau_sub[jj] < tau[ii])\n			jj++;\n'
     '		if (jj < 1 || jj >= srtcount) {\n'),
    ('			theta[ii] = ((tau[ii] - tau_sub[jj - 1]) * theta_sub[jj - 1] + (tau_sub[jj] - tau[ii]) * theta_sub[jj]) /\n',
     '			theta[ii] = ((tau_sub[jj] - tau[ii]) * theta_sub[jj - 1] + (tau[ii] - tau_sub[jj - 1]) * theta_sub[jj]) /\n'),
    ('	while (theta_eap[jj] < theta[0])\n		jj++;\n',
     '	while (jj < n && theta_eap[jj] < theta[0])\n		jj++;\n'),
    ('		while (theta_eap[jj] < theta[ii])\n			jj++;\n'
     '		if (jj < 1 || jj > srtcount) {\n',
     '		while (jj < n && theta_eap[jj] < theta[ii])\n			jj++;\n'
     '		if (jj < 1 || jj >= n) { /* IGML: the EAP table is n long (was checked against srtcount) */\n'),
    ('			real = (Geap[(jj - 1) * 2] * (theta[ii] - theta_eap[jj - 1]) + Geap[jj * 2] * (theta_eap[jj] - theta[ii])) / dtheta;\n'
     '			imag = (Geap[(jj - 1) * 2 + 1] * (theta[ii] - theta_eap[jj - 1]) + Geap[jj * 2 + 1] * (theta_eap[jj] - theta[ii])) /\n'
     '			       dtheta;\n',
     '			/* IGML: weights were swapped */\n'
     '			real = (Geap[(jj - 1) * 2] * (theta_eap[jj] - theta[ii]) + Geap[jj * 2] * (theta[ii] - theta_eap[jj - 1])) / dtheta;\n'
     '			imag = (Geap[(jj - 1) * 2 + 1] * (theta_eap[jj] - theta[ii]) + Geap[jj * 2 + 1] * (theta[ii] - theta_eap[jj - 1])) /\n'
     '			       dtheta;\n'),
]
for old, new in reps:
    assert s.count(old) == 1, (s.count(old), old[:90])
    s = s.replace(old, new)
(F / "make_s1a_tops.c").write_text(s)
print("patched make_s1a_tops.c")
