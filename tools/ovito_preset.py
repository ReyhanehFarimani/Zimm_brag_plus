#!/usr/bin/env python
"""OVITO preset for the Zimm-Bragg+ trajectories (<out_prefix>_conf.xyz, extended-xyz written by
logging.cpp: species C/R/L, pos, orientation quaternion, tangent, registry, rod_L, rod_D, spin), in
the style of the CG figure of the grant (render_chain_cg.py): helices as spherocylinders coloured
red (R) / amber (L) with the dumped rod length and diameter, coils as blue spheres, and the registry
direction as a dark arrow on every helical rod (user 2026-09-23).

Three ways to use it:

  1. OVITO GUI (Pro): Add modification -> Python script -> load this file.  The modify() function
     below becomes the modifier; nothing else to set.  On this machine start the GUI with
     QT_QPA_PLATFORM=xcb (the Wayland viewport crashes).
  2. Session state: `ovitos tools/ovito_preset.py --save-state zimm_preset.ovito <traj>` writes a
     .ovito file with the trajectory loaded and the modifier applied; open it in the GUI and, to look at
     another trajectory, just change the file of the pipeline source.
  3. Headless render: `ovitos tools/ovito_preset.py <traj> --out frame.png [--frame -1] [--size 1600 900]`
     (never opens a viewport; the same camera logic as render_chain_cg.py: principal frame of the
     chain, view down the flattest axis).

Colours: R #e34948, L #eda100, C #2a78d6, registry arrow #262626 -- the house palette.
"""
import sys

import numpy as np

COLOR = {"R": (0.890, 0.286, 0.282), "L": (0.929, 0.631, 0.000), "C": (0.165, 0.471, 0.839),
         "X": (0.55, 0.55, 0.55)}       # X = the star core (a sphere of radius rod_D / 2 = core_radius)
COLOR_REG = (0.15, 0.15, 0.15)
REG_LEN = 1.6          # arrow protrusion beyond the rod surface [code length units = a]
REG_WIDTH = 0.18
BEAD = 1.0             # the bead's own diameter, added to the dumped screw length and diameter as in the
                       # CG figure (rodL = axial span + a, rodD = 2 R_screw + a; render_chain_cg.py)
COIL_D = 2.4           # coil sphere diameter (about one C-C bond: consecutive coil spheres touch)


def modify(frame, data):
    """OVITO Python modifier: shapes, colours and registry arrows from the dumped columns."""
    from ovito.data import ParticleType
    from ovito.vis import ParticlesVis, VectorVis
    parts = data.particles_
    n = parts.count
    types = parts.particle_types_
    name_of = {t.id: t.name for t in types.types}
    species = np.array([name_of.get(int(i), "C") for i in types])
    helix = np.isin(species, ("R", "L"))
    # one particle type per class with its shape and colour
    for t in types.types:
        nm = t.name if t.name in COLOR else "C"
        tm = types.make_mutable(t)                                  # types are shared with the source
        tm.color = COLOR[nm]
        tm.shape = ParticlesVis.Shape.Spherocylinder if nm in ("R", "L") else ParticlesVis.Shape.Sphere
    rod_L = np.asarray(parts["rod_L"]); rod_D = np.asarray(parts["rod_D"])
    # spherocylinder: OVITO takes the cylinder radius from Aspherical Shape.X and the cylinder length
    # from .Z, capped by hemispheres, so the axial span is rod_L when Z = rod_L - rod_D
    L = rod_L + BEAD; D = rod_D + BEAD
    shape = np.zeros((n, 3))
    shape[helix, 0] = shape[helix, 1] = 0.5 * D[helix]
    shape[helix, 2] = np.maximum(L[helix] - D[helix], 0.0)
    parts.create_property("Aspherical Shape", data=shape)
    core = species == "X"
    parts.create_property("Radius", data=np.where(helix, 0.5 * D, np.where(core, 0.5 * rod_D, 0.5 * COIL_D)))
    parts.create_property("Color", data=np.array([COLOR[s if s in COLOR else "C"] for s in species]))
    # registry arrow from the rod axis, helices only (a coil's registry is a dummy)
    reg = np.asarray(parts["registry"], dtype=float).copy()
    reg[~helix] = 0.0
    reg[helix] *= (0.5 * D[helix] + REG_LEN)[:, None]
    p = parts.create_property("Registry arrow", data=reg)
    p.vis = VectorVis(alignment=VectorVis.Alignment.Base, color=COLOR_REG, width=REG_WIDTH,
                      scaling=1.0, enabled=True)
    data.attributes["n_helix"] = int(helix.sum())


def _pipeline(traj):
    from ovito.io import import_file
    from ovito.modifiers import PythonScriptModifier
    pipeline = import_file(traj)
    pipeline.modifiers.append(PythonScriptModifier(function=modify))
    return pipeline


def render(traj, out, frame=-1, size=(1600, 900)):
    """Headless render of one frame: the chain is rotated into its principal frame (x = long axis,
    z = flattest, a proper rotation so no handedness flips) and viewed down z, as in render_chain_cg.py."""
    from ovito.modifiers import AffineTransformationModifier
    from ovito.vis import TachyonRenderer, Viewport
    pipeline = _pipeline(traj)
    nframes = pipeline.source.num_frames
    fr = nframes + frame if frame < 0 else frame
    data = pipeline.compute(fr)
    pos = np.asarray(data.particles.positions)
    cen = pos.mean(0)
    _, _, vt = np.linalg.svd(pos - cen, full_matrices=False)
    rot = vt.copy()
    if np.linalg.det(rot) < 0:
        rot[2] *= -1.0                          # rotation, not reflection (would flip every handedness)
    M = np.zeros((3, 4)); M[:, :3] = rot; M[:, 3] = -rot @ cen
    pipeline.modifiers.append(AffineTransformationModifier(transformation=M.tolist(), operate_on={"particles", "vector_properties"}))
    half = np.abs((pos - cen) @ rot.T).max(0)
    aspect = size[0] / size[1]
    pipeline.add_to_scene()
    vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(0.0, 0.0, -1.0), camera_pos=(0.0, 0.0, float(half[2] + 50.0)))
    vp.fov = float(max(half[1], half[0] / aspect) * 1.12 + 2.0)
    renderer = TachyonRenderer(shadows=False, direct_light_intensity=1.1, ambient_occlusion=True,
                               ambient_occlusion_brightness=0.85)
    vp.render_image(filename=out, size=tuple(size), renderer=renderer, background=(1.0, 1.0, 1.0), frame=fr)
    pipeline.remove_from_scene()
    print(f"wrote {out}  (frame {fr} of {nframes}, {int(data.attributes['n_helix'])} helical residues)")


def save_state(traj, path):
    import ovito
    pipeline = _pipeline(traj)
    pipeline.add_to_scene()
    ovito.scene.save(path)
    print(f"wrote session state {path}  (source: {traj}; change the file in the GUI to view another run)")


def _run_as_script():
    """True when executed as a file by ovitos/python.  OVITO's GUI Python-script modifier also executes
    this file with __name__ == "__main__" (sys.argv = [its path], __file__ = None), and the argparse
    block below must not run there: its sys.exit aborts the script before the GUI picks up modify()."""
    return __name__ == "__main__" and bool(globals().get("__file__"))


if _run_as_script():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("traj")
    ap.add_argument("--out", default="")
    ap.add_argument("--save-state", default="")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--size", type=int, nargs=2, default=(1600, 900))
    a = ap.parse_args()
    if a.save_state:
        save_state(a.traj, a.save_state)
    if a.out:
        render(a.traj, a.out, a.frame, tuple(a.size))
    if not a.out and not a.save_state:
        sys.exit("give --out <png> and/or --save-state <file.ovito>")
