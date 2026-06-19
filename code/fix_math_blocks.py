#!/usr/bin/env python3
"""Make display-math GitHub-safe: ensure a blank line before an opening `$$`
and after a closing `$$` (and around single-line `$$...$$`). VS Code is lenient;
GitHub's MathJax requires the blank lines or it renders the LaTeX as raw text."""
import sys, glob, os

def fix(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    out, inside = [], False
    for i, line in enumerate(lines):
        t = line.strip()
        nxt = lines[i+1] if i+1 < len(lines) else ""
        if t == "$$" and not inside:                      # opening delimiter
            if out and out[-1].strip() != "": out.append("")
            out.append(line); inside = True
        elif t == "$$" and inside:                        # closing delimiter
            out.append(line); inside = False
            if nxt.strip() != "": out.append("")
        elif (not inside) and t.startswith("$$") and t.endswith("$$") and len(t) > 3:
            if out and out[-1].strip() != "": out.append("")  # single-line block
            out.append(line)
            if nxt.strip() != "": out.append("")
        else:
            out.append(line)
    new = "\n".join(out)
    if new != "\n".join(lines):
        open(path, "w", encoding="utf-8").write(new)
        return True
    return False

course = os.path.join(os.path.dirname(__file__), "..", "course")
changed = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(course, "*.md"))) if fix(f)]
print("normalized:", ", ".join(changed) if changed else "(none)")
