#!/usr/bin/env python3
"""Make display math GitHub-safe.
GitHub's MathJax reliably renders a block ONLY as a single-line `$$ ... $$` (its
documented form). Multi-line `$$\n...\n$$` blocks are left as raw text on GitHub
(underscores escaped, `\\;` shown as `;`). So we:
  1. collapse every multi-line `$$` block to ONE line: `$$ <joined> $$`,
     preserving `\\` (matrix/array row breaks);
  2. drop cosmetic spacing macros `\\;` and `\\,` (these are what render as the
     stray ';' and ',' the reader sees) — replace with a normal space;
  3. keep a blank line before and after each block.
VS Code / our PDF build are unaffected (both handle either form)."""
import glob, os, re

def fix(path):
    src = open(path, encoding="utf-8").read()
    lines = src.split("\n")
    out, i = [], 0
    while i < len(lines):
        if lines[i].strip() == "$$":                      # opening delimiter
            body, j = [], i + 1
            while j < len(lines) and lines[j].strip() != "$$":
                body.append(lines[j].strip()); j += 1
            content = " ".join(x for x in body if x)
            if out and out[-1].strip() != "": out.append("")
            out.append("$$ " + content + " $$")
            out.append("")
            i = j + 1
            if i < len(lines) and lines[i].strip() == "": i += 1   # absorb one trailing blank
        else:
            out.append(lines[i]); i += 1
    text = "\n".join(out)
    # drop cosmetic spacing macros (they live only in math); keep \\ row breaks.
    # NOTE: do NOT collapse other whitespace — would corrupt code-block indentation.
    text = text.replace("\\;", " ").replace("\\,", " ")
    if text != src:
        open(path, "w", encoding="utf-8").write(text)
        return True
    return False

course = os.path.join(os.path.dirname(__file__), "..", "course")
changed = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(course, "*.md"))) if fix(f)]
print("normalized:", ", ".join(changed) if changed else "(none)")
