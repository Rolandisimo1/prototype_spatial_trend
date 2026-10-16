"""Build the related-work review for the fleet trend paper.

Every in-text citation is generated from lit_refs.json, which holds DOIs
resolved against OpenAlex -- no citation is typed by hand, because a mistyped
DOI in a manuscript is both invisible and permanent. One reference in the
first pass resolved to the wrong paper on a short title fragment ("Observing
the Observers"), so the resolver checks the expected first author before
accepting a match.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REFS = os.path.join(os.path.dirname(HERE), "lit_refs.json")


def load_refs(path=REFS):
    with open(path) as f:
        return json.load(f)


def cite(refs, key, text=None):
    """Markdown link whose target is the resolved DOI. Raises on an unknown key."""
    if key not in refs:
        raise KeyError(f"no resolved reference for {key!r}; "
                       f"have {sorted(refs)}")
    r = refs[key]
    # Surnames with particles ("van Strien", "de Vries") lose the particle if we
    # just take the last token, so a reference may carry an explicit override.
    surname = r.get("surname") or (r["author"].split()[-1] if r["author"] else "Anon")
    label = text or f"{surname} et al. {r['year']}"
    # Older Ecology DOIs embed parentheses and square brackets
    # (e.g. 10.1890/0012-9658(2002)083[2248:esorwd]2.0.co;2). Bare in a markdown
    # link target those close the link early and pandoc drops the hyperlink
    # silently. Angle-bracket form is literal and safe for every DOI.
    return f"[{label}](<https://doi.org/{r['doi']}>)"


def bibliography(refs):
    lines = []
    for k in sorted(refs, key=lambda k: (refs[k]["author"].split()[-1] if refs[k]["author"] else "")):
        r = refs[k]
        lines.append(f"- {r['author']} et al. ({r['year']}). *{r['title']}*. "
                     f"{r['venue'] or 'n/a'}. "
                     f"[doi:{r['doi']}](<https://doi.org/{r['doi']}>)")
    return "\n".join(lines)
