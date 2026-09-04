#!/usr/bin/env python
"""Render every figure for the diel-activity paper and report.

    python figures/make_figures.py --outdir figs

Inputs are resolved from the artifact store by filename, so this runs from a clean checkout
without any hand-staging. Add a new figure by writing a function in dielfigs.py and registering
it in dielfigs.FIGURES; do not draw figures ad hoc in a session.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dielfigs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="figs")
    ap.add_argument("--only", nargs="*", default=None,
                    help="figure names to render; default all")
    ap.add_argument("--no-style", action="store_true",
                    help="skip the figure-style hook (use plain matplotlib defaults)")
    args = ap.parse_args()

    try:
        host  # provided by the Claude Science kernel
    except NameError:
        raise SystemExit("This runner needs the `host` object to resolve input artifacts. "
                         "Run it inside the analysis kernel, or pass explicit paths to "
                         "dielfigs.load_all().")

    style = None
    if not args.no_style:
        try:
            style = lambda: apply_figure_style(sizes=(9, 8, 7))  # noqa: E731
        except NameError:
            style = None

    made = dielfigs.render_all(host, outdir=args.outdir, style=style, only=args.only)
    for name, path in sorted(made.items()):
        print(f"{name:<32} {path}")
    print(f"\n{len(made)} figures written to {args.outdir}")


if __name__ == "__main__":
    main()
