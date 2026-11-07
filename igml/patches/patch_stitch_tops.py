"""Make our patched copy of GMTSAR's stitch_tops.c (2026-10-28): drop the second subtraction of line_offset for the last
image. The lab's GMTSAR (/usr/local/GMTSAR) is not touched; src/stitch_tops.orig.c is its copy."""
from pathlib import Path

F = Path("/zdata/datasets/IGML/gmtsar_fix/src")
s = (F / "stitch_tops.orig.c").read_text()
old = ('    if (line_offset!= 0) {\n'
       '        //printf("Need to shift lines by %d\\n",line_offset);\n'
       '        n_start = n_start - line_offset;\n'
       '    }\n')
assert s.count(old) == 1, s.count(old)
new = ('    /* IGML 2026-10-28: line_offset is already in line_offset_total, applied just above; subtracting it again\n'
       '       placed the last image |line_offset| lines off whenever the overlap cross-correlation found a shift. */\n')
(F / "stitch_tops.c").write_text(s.replace(old, new))
print("patched")
