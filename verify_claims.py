#!/usr/bin/env python3
"""Verify that written claims are supported by quotes that exist in their sources.

Portable. Drop into any repo, point it at a directory of claim files, and it will tell you
which claims rest on something checkable and which rest on nothing. No repo-specific
assumptions beyond the file format below.

    uv run --with pydantic --with pyyaml verify_claims.py claims/
    uv run --with pydantic --with pyyaml verify_claims.py claims/ --strict   # exit 1 on failure

FILE FORMAT --- one YAML per source document:

    source:
      citation: wang2023interpretability
      local: papers/wang_2023.pdf        # relative to --root, default cwd
      sha256: 3f9a...                    # optional; fails loudly if the file changes
      extract_cmd: pdftotext -layout     # pinned, so line numbers mean something
    claims:
      any-key-you-like:
        statement: "The circuit recovers 87% of the logit difference."
        verified: true                   # rejected unless quotes are present
        quotes:
          - exact: "C achieves 87% of the performance of M."
            section: "§4"                # at least one of section / page / line
        notes: "anything, no length limit"

WHY IT IS BUILT THIS WAY

Quotes are anchored as a W3C Web Annotation TextQuoteSelector: `exact`, plus `prefix` and
`suffix` when `exact` alone is ambiguous. Offsets are recorded but advisory --- they move
when the extraction tool changes, the selector does not.

Matching is two-tier. Exact first, on text normalised for the differences extraction
invents: combining diacritics, dotless i, smart quotes, the Unicode minus, hyphenation
across line breaks, spaces LaTeX pads into inline math. If that fails, the alphanumeric
skeleton, which removes everything else --- and a match found only there is reported as
`loose`, never as exact, because the skeleton cannot tell `a - b` from `a + b`.

A quote that resolves nowhere is a failure, not a warning, and is never deleted.
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import re
import subprocess
import sys
import unicodedata

import yaml
from pydantic import BaseModel, Field, model_validator

CACHE = pathlib.Path("/tmp/verify_claims_extracts")


def fold(s: str) -> str:
    """Normalise away differences no reader would call differences."""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    for a, b in (("‘", "'"), ("’", "'"), ("“", '"'), ("”", '"'),
                 ("−", "-"), ("–", "-"), ("—", "-"), (" ", " "),
                 ("ı", "i"), ("ﬁ", "fi"), ("ﬂ", "fl"), ("_", ""), ("-", "")):
        s = s.replace(a, b)
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\s*([(){}\[\]|])\s*", r"\1", s)
    return s.strip().lower()


def skeleton(s: str) -> str:
    """Letters and digits only. Last resort; loses operators, so matches are flagged."""
    return re.sub(r"[^0-9a-z]", "", fold(s))


class Source(BaseModel):
    citation: str
    local: str | None = None
    sha256: str | None = None
    extract_cmd: str = "pdftotext -layout"

    def text(self, root: pathlib.Path) -> str:
        if not self.local:
            raise ValueError(f"{self.citation}: no local source")
        f = root / self.local
        if not f.exists():
            raise ValueError(f"{self.citation}: missing {self.local}")
        if self.sha256:
            got = hashlib.sha256(f.read_bytes()).hexdigest()
            if got != self.sha256:
                raise ValueError(f"{self.citation}: source changed "
                                 f"({self.sha256[:12]} -> {got[:12]})")
        if f.suffix in (".txt", ".md"):
            return f.read_text(errors="replace")
        CACHE.mkdir(exist_ok=True)
        out = CACHE / (hashlib.sha256(str(f).encode()).hexdigest()[:16] + ".txt")
        if not out.exists():
            subprocess.run(self.extract_cmd.split() + [str(f), str(out)],
                           check=True, capture_output=True)
        return out.read_text(errors="replace")


class Quote(BaseModel):
    exact: str = Field(min_length=12)
    prefix: str = ""
    suffix: str = ""
    section: str | None = None
    page: int | None = None
    line: int | None = None

    @model_validator(mode="after")
    def _needs_a_locator(self):
        if self.section is None and self.page is None and self.line is None:
            raise ValueError("a quote needs at least one of section, page or line")
        return self

    def resolve(self, text: str) -> str:
        flat = fold(text)
        n = flat.count(fold(self.exact))
        if n == 1:
            return "ok"
        if n > 1:
            anchor = fold(self.prefix + self.exact + self.suffix)
            return "ok" if flat.count(anchor) == 1 else "ambiguous"
        bones = skeleton(text)
        k = bones.count(skeleton(self.exact))
        return "loose" if k == 1 else ("ambiguous" if k > 1 else "not found")


class Claim(BaseModel):
    statement: str
    verified: bool = False
    quotes: list[Quote] = []
    notes: str = ""

    @model_validator(mode="after")
    def _verified_needs_quotes(self):
        if self.verified and not self.quotes:
            raise ValueError("verified with no quotes; that is an assertion, not a check")
        return self


class Document(BaseModel):
    source: Source
    claims: dict[str, Claim]

    @classmethod
    def load(cls, path: pathlib.Path) -> "Document":
        return cls.model_validate(yaml.safe_load(path.read_text()))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("directory", help="directory of claim YAML files")
    ap.add_argument("--root", default=".", help="base for `local` paths (default: cwd)")
    ap.add_argument("--strict", action="store_true", help="exit 1 if any quote fails")
    args = ap.parse_args()
    root = pathlib.Path(args.root).resolve()

    print(f"{'file':<28}{'claims':>7}{'quoted':>8}{'exact':>7}{'loose':>7}{'FAIL':>6}  notes")
    failures = 0
    for path in sorted(pathlib.Path(args.directory).glob("*.yaml")):
        try:
            doc = Document.load(path)
            text = doc.source.text(root)
        except Exception as e:                                  # noqa: BLE001
            print(f"{path.stem:<28}{'':>7}{'':>8}{'':>7}{'':>7}{'':>6}  {e}")
            failures += 1
            continue
        tally = {"ok": 0, "loose": 0, "not found": 0, "ambiguous": 0}
        for c in doc.claims.values():
            for q in c.quotes:
                tally[q.resolve(text)] += 1
        bad = tally["not found"] + tally["ambiguous"]
        failures += bad
        quoted = sum(1 for c in doc.claims.values() if c.quotes)
        note = "" if quoted else "no claim is supported by a quote"
        print(f"{path.stem:<28}{len(doc.claims):>7}{quoted:>8}{tally['ok']:>7}"
              f"{tally['loose']:>7}{bad:>6}  {note}")
    print(f"\n{failures} failing")
    return 1 if (args.strict and failures) else 0


if __name__ == "__main__":
    sys.exit(main())
