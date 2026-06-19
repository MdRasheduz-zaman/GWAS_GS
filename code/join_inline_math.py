#!/usr/bin/env python3
"""Join inline math that wraps across a soft line break onto ONE line.
GitHub only renders inline `$...$` that lives on a single source line; a span
broken by a newline is shown as raw text (looks like the subscript `_` was eaten).
We merge consecutive prose lines until every inline `$...$` is balanced. Markdown
already treats a soft break as a space, so rejoining does not change the prose."""
import glob, os, re

def dollars(s):
    return re.sub(r"`[^`]*`", "", s).replace("$$", "").count("$")

def fix(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    out, i, fence = [], 0, False
    while i < len(lines):
        l = lines[i]
        st = l.strip()
        if st.startswith("```"):
            fence = not fence; out.append(l); i += 1; continue
        if fence or st == "":
            out.append(l); i += 1; continue
        cur, cnt = l, dollars(l)
        while (cnt % 2 == 1 and i + 1 < len(lines)
               and lines[i+1].strip() != "" and not lines[i+1].strip().startswith("```")):
            i += 1
            cur = cur + " " + lines[i].strip()
            cnt += dollars(lines[i])
        out.append(cur); i += 1
    new = "\n".join(out)
    if new != "\n".join(lines):
        open(path, "w", encoding="utf-8").write(new); return True
    return False

course = os.path.join(os.path.dirname(__file__), "..", "course")
changed = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(course, "*.md"))) if fix(f)]
print("joined inline math in:", ", ".join(changed) if changed else "(none)")
