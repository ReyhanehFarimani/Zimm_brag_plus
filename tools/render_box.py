#!/usr/bin/env python
"""Headless snapshots of a PERIODIC-BOX trajectory (user 2026-09-24: "add a few trajectory of the dense system"):
positions wrapped into the cell, the residues drawn with the house preset (tools/ovito_preset.py: helices as red R /
amber L spherocylinders, coils as blue spheres, registry arrows), orthographic view down z of the whole box, and
optionally a slab of thickness SLAB a through the middle so individual chains can be followed.
  ~/venv/bin/python tools/render_box.py <traj.xyz> --frames 0 3 6 --out-prefix runs/t45/box500/snapshots/J9 [--slab 8] [--size 1400]
"""
import argparse
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("preset", os.path.join(HERE, "ovito_preset.py")); preset = importlib.util.module_from_spec(spec); spec.loader.exec_module(preset)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traj"); ap.add_argument("--frames", type=int, nargs="+", default=[-1]); ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--slab", type=float, default=0.0, help="also render a slab |z - L/2| < slab/2 (a)"); ap.add_argument("--size", type=int, default=1400)
    a = ap.parse_args()
    from ovito.io import import_file
    from ovito.modifiers import WrapPeriodicImagesModifier, ExpressionSelectionModifier, DeleteSelectedModifier, PythonScriptModifier
    from ovito.vis import TachyonRenderer, Viewport
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    pipe = import_file(a.traj)                     # extended xyz: the Properties= header names the columns, Lattice= the cell
    pipe.modifiers.append(WrapPeriodicImagesModifier())
    pipe.modifiers.append(PythonScriptModifier(function=preset.modify))
    nfr = pipe.source.num_frames
    data0 = pipe.compute(0); L = float(data0.cell[0, 0])
    renderer = TachyonRenderer(shadows=False, direct_light_intensity=1.1, ambient_occlusion=True, antialiasing_samples=12)
    for fr in a.frames:
        fr = fr if fr >= 0 else nfr + fr
        data = pipe.compute(fr)
        sweep = data.attributes.get("sweep", fr)
        for tag, slab in (("box", 0.0), ("slab", a.slab)):
            if tag == "slab" and slab <= 0:
                continue
            if slab > 0:
                pipe.modifiers.append(ExpressionSelectionModifier(expression=f"abs(Position.Z - {L / 2}) > {slab / 2}"))
                pipe.modifiers.append(DeleteSelectedModifier())
            pipe.add_to_scene()
            vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(0.0, 0.0, -1.0), camera_pos=(L / 2, L / 2, L + 10.0), fov=L * 0.56)
            out = f"{a.out_prefix}_{tag}_sweep{int(sweep)}.png"
            vp.render_image(filename=out, size=(a.size, a.size), renderer=renderer, background=(1.0, 1.0, 1.0), frame=fr)
            pipe.remove_from_scene()
            if slab > 0:
                pipe.modifiers.pop(); pipe.modifiers.pop()
            print(out, flush=True)


if __name__ == "__main__":
    main()
