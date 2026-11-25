# fusion-metrics — Copyright (c) 2026 Vu Viet Anh (Vũ Việt Anh) <vietanhresearcher@gmail.com>
# MIT License (third-party-derived metrics: see NOTICE). https://github.com/VuVietAnh2309/fusion-metrics
"""Minimal example: score one fused image with the Python API.

    python examples/score_images.py MRI.png PET.png fused.png [preset]

preset is one of: image, information, feature, similarity, all (default).
The same thing from the shell:  fusion-metrics MRI.png PET.png fused.png --preset feature
"""
import sys

import fusion_metrics as fm

a, b, f = sys.argv[1:4]
preset = sys.argv[4] if len(sys.argv) > 4 else "all"

scores = fm.evaluate_files(a, b, f, preset)   # reads with MATLAB imread + rgb2gray
for name, value in scores.items():
    print(f"{name:12s} {value:.6f}")
