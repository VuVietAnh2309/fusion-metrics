# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Command-line interface.

Score one fused image:
    fusion-metrics A.png B.png fused.png [--preset feature | --metrics Q_ABF Q_MI] [--json]

Score a whole method (files matched by name across three folders) and write a CSV:
    fusion-metrics --dirs srcA/ srcB/ fused/ --preset all --csv results.csv
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

from . import METRICS, PRESETS, __version__, evaluate, read_gray, read_gray_cv

EXTS = (".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff")


def _parser():
    p = argparse.ArgumentParser(prog="fusion-metrics",
                                description="Objective quality metrics for image fusion (MATLAB-faithful).")
    p.add_argument("images", nargs="*", metavar="IMG", help="source A, source B and fused image")
    p.add_argument("--dirs", nargs=3, metavar=("DIR_A", "DIR_B", "DIR_F"),
                   help="score every fused image in DIR_F against same-named files in DIR_A / DIR_B")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--preset", choices=sorted(PRESETS), help="metric group (default: all)")
    g.add_argument("--metrics", nargs="+", choices=list(METRICS), metavar="NAME", help="explicit metric list")
    p.add_argument("--reader", choices=["matlab", "cv"], default="matlab",
                   help="gray conversion: MATLAB rgb2gray (default) or OpenCV")
    p.add_argument("--json", action="store_true", help="print JSON instead of a table")
    p.add_argument("--csv", metavar="FILE", help="with --dirs: write per-image scores to FILE")
    p.add_argument("--list", action="store_true", help="list metrics and presets, then exit")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def _find(folder, stem):
    for ext in EXTS:
        path = os.path.join(folder, stem + ext)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"no image named '{stem}' in {folder}")


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.list:
        print("metrics:", " ".join(METRICS))
        for k, v in PRESETS.items():
            print(f"preset {k:12s}", " ".join(v))
        return 0
    names = args.metrics or PRESETS[args.preset or "all"]
    read = read_gray if args.reader == "matlab" else read_gray_cv

    if args.dirs:
        da, db, df = args.dirs
        stems = sorted(os.path.splitext(f)[0] for f in os.listdir(df) if f.lower().endswith(EXTS))
        rows = []
        for i, s in enumerate(stems, 1):
            r = evaluate(read(_find(da, s)), read(_find(db, s)), read(_find(df, s)), names)
            rows.append({"image": s, **r})
            print(f"\r[{i}/{len(stems)}] {s}", end="", file=sys.stderr)
        print(file=sys.stderr)
        if args.csv:
            with open(args.csv, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["image", *names])
                w.writeheader()
                w.writerows(rows)
        mean = {n: float(np.mean([r[n] for r in rows])) for n in names}
        result = {"n_images": len(rows), "mean": mean}
    elif len(args.images) == 3:
        result = evaluate(*(read(p) for p in args.images), names)
    else:
        _parser().error("give three images (A B fused) or --dirs DIR_A DIR_B DIR_F")

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        scores = result.get("mean", result)
        if "n_images" in result:
            print(f"mean over {result['n_images']} images")
        for k, v in scores.items():
            print(f"{k:12s} {v:.6f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
