#!/usr/bin/env python3
"""Convert every display block `$$ ... $$` to a GitHub ```math fenced block.
Why: GitHub markdown-processes the *content* of $$...$$ (eating \\, \#, \;, and
requiring blank lines), which keeps breaking equations. A ```math code fence is
rendered verbatim by GitHub's math engine — no markdown mangling, no blank-line
fuss, and \\ row-breaks work (so real matrices render). Inline $...$ is untouched."""
import glob, os, re

def convert(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    out, i = [], 0
    while i < len(lines):
        st = lines[i].strip()
        # pass existing code fences through untouched (incl. already-converted ```math)
        if st.startswith("```"):
            out.append(lines[i]); i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                out.append(lines[i]); i += 1
            if i < len(lines): out.append(lines[i]); i += 1
            continue
        # standalone multi-line $$ ... $$
        if st == "$$":
            body, j = [], i + 1
            while j < len(lines) and lines[j].strip() != "$$":
                body.append(lines[j].strip()); j += 1
            content = " ".join(x for x in body if x)
            if out and out[-1].strip() != "": out.append("")
            out += ["```math", content, "```", ""]
            i = j + 1
            if i < len(lines) and lines[i].strip() == "": i += 1
            continue
        # single-line $$ ... $$ on its own line
        m = re.match(r'^\$\$(.+?)\$\$$', st)
        if m:
            if out and out[-1].strip() != "": out.append("")
            out += ["```math", m.group(1).strip(), "```", ""]
            i += 1
            if i < len(lines) and lines[i].strip() == "": i += 1
            continue
        out.append(lines[i]); i += 1
    new = "\n".join(out)
    if new != "\n".join(lines):
        open(path, "w", encoding="utf-8").write(new); return True
    return False

course = os.path.join(os.path.dirname(__file__), "..", "course")
changed = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(course, "*.md"))) if convert(f)]
print("converted to ```math:", ", ".join(changed) if changed else "(none)")
