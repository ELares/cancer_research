"""Render an animated axial-slice GIF (and MP4) of sim-tme-3d's
per-step trajectory captured by `sim-tme-3d --snapshot` (#193).

Reads:
    simulations/output/tme-3d/trajectory_dead.npy   (u8, n_steps x layers x rows x cols)
    simulations/output/tme-3d/trajectory_damp.npy   (f32, same shape)
    simulations/output/tme-3d/trajectory_lp.npy     (f32, same shape)
    simulations/output/tme-3d/trajectory_meta.json  (condition + grid descriptor)

Writes:
    simulations/output/tme-3d/trajectory_axial.gif
    simulations/output/tme-3d/trajectory_axial.mp4  (if ffmpeg available)

Visualization:
    Three panels showing a central mid-plane slice (the row axis fixed
    at rows/2) of the spheroid at each step:
      1. Dead-cell mask    (grayscale background + red dead cells)
      2. DAMP field        (inferno colormap, dynamic range)
      3. LP (lipid perox.) (viridis colormap, dynamic range)
    180 frames at 15 fps = 12s animation.

Run:
    python3 scripts/render_tme_3d_trajectory.py
    python3 scripts/render_tme_3d_trajectory.py --fps 30 --no-mp4
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap


REPO_ROOT = Path(__file__).resolve().parent.parent
TRAJ_DIR = REPO_ROOT / "simulations" / "output" / "tme-3d"
EXPECTED_SCHEMA_VERSION = 1


def _load_trajectory(
    traj_dir: Path,
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, "np.ndarray | None", "np.ndarray | None", dict
]:
    """Load the trajectory arrays + metadata; assert schema version.

    Returns (dead, damp, lp, persister, subclone, meta). `persister` (#241) and
    `subclone` (#242) are optional overlays, present only when the run used the
    corresponding model (`--snapshot=persister` / `--snapshot=clonal`); `None`
    otherwise. Raises SystemExit on missing required files or schema mismatch
    with a clear, actionable message.
    """
    required = [
        "trajectory_dead.npy",
        "trajectory_damp.npy",
        "trajectory_lp.npy",
        "trajectory_meta.json",
    ]
    missing = [f for f in required if not (traj_dir / f).exists()]
    if missing:
        raise SystemExit(
            f"ERROR: missing trajectory file(s) in {traj_dir}: {missing}\n"
            f"Run `cargo run --release -p sim-tme-3d -- --snapshot` first."
        )

    meta = json.loads((traj_dir / "trajectory_meta.json").read_text())
    v = meta.get("schema_version")
    if v != EXPECTED_SCHEMA_VERSION:
        raise SystemExit(
            f"ERROR: trajectory_meta.json schema_version={v!r}, expected "
            f"{EXPECTED_SCHEMA_VERSION}. Bump the renderer or regenerate "
            f"the trajectory with a matching binary."
        )

    dead = np.load(traj_dir / "trajectory_dead.npy")
    damp = np.load(traj_dir / "trajectory_damp.npy")
    lp = np.load(traj_dir / "trajectory_lp.npy")

    # Shape sanity — the renderer assumes 4-D (steps, layers, rows, cols).
    for name, a in [("dead", dead), ("damp", damp), ("lp", lp)]:
        if a.ndim != 4:
            raise SystemExit(
                f"ERROR: trajectory_{name}.npy has ndim={a.ndim}, expected 4 "
                f"(steps, layers, rows, cols)."
            )
    if not (dead.shape == damp.shape == lp.shape):
        raise SystemExit(
            f"ERROR: trajectory arrays disagree on shape: "
            f"dead={dead.shape}, damp={damp.shape}, lp={lp.shape}."
        )

    # Optional persister-fraction overlay (#241). Only written when the run
    # enabled the persister model (e.g. `--snapshot=persister`); a plain
    # snapshot has no such file and renders the original three panels.
    persister = None
    persister_path = traj_dir / "trajectory_persister.npy"
    if persister_path.exists():
        persister = np.load(persister_path)
        if persister.ndim != 4 or persister.shape != dead.shape:
            raise SystemExit(
                f"ERROR: trajectory_persister.npy shape {persister.shape} does "
                f"not match the other fields {dead.shape}."
            )

    # Optional static subclone-id map (#242). One 3D frame (no time axis),
    # matching a trajectory frame's spatial axes. Present only for the
    # `--snapshot=clonal` preset.
    subclone = None
    subclone_path = traj_dir / "subclone.npy"
    if subclone_path.exists():
        subclone = np.load(subclone_path)
        if subclone.ndim != 3 or subclone.shape != dead.shape[1:]:
            raise SystemExit(
                f"ERROR: subclone.npy shape {subclone.shape} does not match a "
                f"trajectory frame {dead.shape[1:]}."
            )

    return dead, damp, lp, persister, subclone, meta


def _make_dead_cmap() -> ListedColormap:
    """Binary cmap: 0 → very light grey (alive/background), 1 → red (dead).

    Distinct from `Reds` so the spheroid is visible even when there
    are zero dead cells in a slice.
    """
    return ListedColormap([(0.92, 0.92, 0.92, 1.0), (0.86, 0.14, 0.14, 1.0)])


def _render(
    dead: np.ndarray,
    damp: np.ndarray,
    lp: np.ndarray,
    persister: "np.ndarray | None",
    subclone: "np.ndarray | None",
    meta: dict,
    out_dir: Path,
    fps: int,
    skip_mp4: bool,
) -> list[Path]:
    """Build the FuncAnimation and save GIF (+ MP4 unless skip_mp4)."""
    # The Rust grid is stored C-order as (row, col, layer), so np.load
    # yields axes (step, row, col, layer). The spheroid is a centered,
    # isotropic 60³ cube, so fixing any spatial axis at its midpoint
    # gives an equivalent central cross-section through the core.
    n_steps, n_rows, _n_cols, _n_layers = dead.shape
    mid = n_rows // 2  # central mid-plane; the slice spans (col, layer)

    # Compute global color ranges so the animation has stable color scales
    # (otherwise each frame's vmin/vmax shifts and the animation flickers).
    damp_max = max(float(damp.max()), 1e-6)
    lp_max = max(float(lp.max()), 1e-6)
    show_persister = persister is not None
    show_subclone = subclone is not None
    # Panel indices: 0..2 are dead/damp/lp; the two optional panels follow.
    pers_idx = 3 if show_persister else None
    sub_idx = (3 + int(show_persister)) if show_subclone else None

    # Dose-administration steps (#239). Empty for steady-state Constant
    # presets; non-empty for multi-dose / bolus / infusion snapshots, where
    # we annotate the frames so the viewer can see death waves sync to doses.
    dose_steps = set(int(s) for s in meta.get("dose_steps", []))
    # Highlight the dose frame plus a few following frames (the death wave
    # lags the dose), so the marker is visible for a beat rather than 1 frame.
    dose_window = set()
    for d in dose_steps:
        for k in range(d, d + 5):
            dose_window.add(k)

    n_panels = 3 + int(show_persister) + int(show_subclone)
    fig, axes = plt.subplots(
        1, n_panels, figsize=(4.3 * n_panels, 4.5), constrained_layout=True
    )
    cond = meta.get("condition", {})
    dose_caption = (
        f"   doses@{sorted(dose_steps)}" if dose_steps else "   (constant dosing)"
    )
    fig.suptitle(
        f"sim-tme-3d  axial mid-slice  ({cond.get('treatment', '?')}, "
        f"immune={cond.get('immune_mode', '?')}, "
        f"stromal={cond.get('stromal_mode') or 'off'}, "
        f"ph={cond.get('ph_mode') or 'off'}, "
        f"λ_O₂={cond.get('o2_lambda_um', '?')}µm){dose_caption}",
        fontsize=10,
    )

    dead_cmap = _make_dead_cmap()
    im_dead = axes[0].imshow(dead[0, mid], cmap=dead_cmap, vmin=0, vmax=1, origin="lower")
    im_damp = axes[1].imshow(
        damp[0, mid], cmap="inferno", vmin=0, vmax=damp_max, origin="lower"
    )
    im_lp = axes[2].imshow(
        lp[0, mid], cmap="viridis", vmin=0, vmax=lp_max, origin="lower"
    )

    axes[0].set_title("Dead-cell mask")
    axes[1].set_title(f"DAMP field (max={damp_max:.2f})")
    axes[2].set_title(f"LP field (max={lp_max:.2f})")

    # Optional persister panel: drug-tolerant persister fraction (#241), fixed
    # 0..1 scale (it is a fraction) so the build-up of tolerance over the
    # run is read directly off the colorbar.
    im_pers = None
    if show_persister:
        im_pers = axes[pers_idx].imshow(
            persister[0, mid], cmap="magma", vmin=0.0, vmax=1.0, origin="lower"
        )
        axes[pers_idx].set_title("Persister fraction")

    # Optional subclone panel (#242): STATIC (no time axis) — the Voronoi
    # subclone map is fixed for the run. 0 = stroma (grey); 1..K = distinct
    # discrete colors so the patch structure is legible.
    im_sub = None
    if show_subclone:
        k = int(subclone.max())
        sub_colors = [(0.92, 0.92, 0.92, 1.0)] + [
            plt.get_cmap("tab10")(i % 10) for i in range(k)
        ]
        sub_cmap = ListedColormap(sub_colors)
        # vmax = k + 1 so each integer id 0..k maps to one whole color band;
        # tick at each band center for a discrete id legend.
        im_sub = axes[sub_idx].imshow(
            subclone[mid], cmap=sub_cmap, vmin=0, vmax=k + 1, origin="lower"
        )
        axes[sub_idx].set_title(f"Subclone id (K={k})")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(im_damp, ax=axes[1], fraction=0.045, pad=0.04)
    fig.colorbar(im_lp, ax=axes[2], fraction=0.045, pad=0.04)
    if show_persister:
        fig.colorbar(im_pers, ax=axes[pers_idx], fraction=0.045, pad=0.04)
    if show_subclone:
        # Discrete legend: a tick per id, centered in its color band.
        # 0 = stroma; 1..K = subclones.
        k = int(subclone.max())
        cbar = fig.colorbar(im_sub, ax=axes[sub_idx], fraction=0.045, pad=0.04)
        cbar.set_ticks([i + 0.5 for i in range(k + 1)])
        cbar.set_ticklabels(["stroma"] + [str(i) for i in range(1, k + 1)])

    step_text = fig.text(0.5, 0.02, "", ha="center", fontsize=9, family="monospace")

    def update(step: int):
        im_dead.set_data(dead[step, mid])
        im_damp.set_data(damp[step, mid])
        im_lp.set_data(lp[step, mid])
        if im_pers is not None:
            im_pers.set_data(persister[step, mid])
        # Count cumulative dead cells in this slice for a quantitative cue.
        n_dead_slice = int(dead[step, mid].sum())
        # Mark dose frames so multi-dose death waves are visually attributable.
        if step in dose_steps:
            dose_marker = "  💉 DOSE"
        elif step in dose_window:
            dose_marker = "  💉 ···"
        else:
            dose_marker = ""
        step_text.set_text(
            f"step {step + 1:3d}/{n_steps}    dead-in-slice={n_dead_slice}{dose_marker}"
        )
        # Red frame border on the dose step itself — a hard-to-miss cue.
        border_on = step in dose_steps
        for ax in axes:
            for spine in ax.spines.values():
                spine.set_edgecolor("red" if border_on else "none")
                spine.set_linewidth(3.0 if border_on else 0.0)
        artists = [im_dead, im_damp, im_lp, step_text]
        if im_pers is not None:
            artists.append(im_pers)
        return artists

    interval_ms = max(1, int(1000.0 / fps))
    anim = animation.FuncAnimation(
        fig, update, frames=n_steps, interval=interval_ms, blit=False
    )

    written: list[Path] = []
    gif_path = out_dir / "trajectory_axial.gif"
    print(f"Writing {gif_path.relative_to(REPO_ROOT)} ({n_steps} frames @ {fps} fps)…")
    anim.save(gif_path, writer=animation.PillowWriter(fps=fps))
    written.append(gif_path)

    if not skip_mp4:
        mp4_path = out_dir / "trajectory_axial.mp4"
        try:
            print(f"Writing {mp4_path.relative_to(REPO_ROOT)} (ffmpeg)…")
            anim.save(mp4_path, writer=animation.FFMpegWriter(fps=fps, bitrate=2400))
            written.append(mp4_path)
        except (FileNotFoundError, RuntimeError) as e:
            # ffmpeg not on PATH or matplotlib's writer rejected the
            # invocation. GIF still landed, so this is degraded-output,
            # not a hard fail.
            print(f"  skipped MP4 ({type(e).__name__}: {e}); GIF still produced.")

    plt.close(fig)
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--traj-dir",
        type=Path,
        default=TRAJ_DIR,
        help=f"Trajectory directory (default: {TRAJ_DIR.relative_to(REPO_ROOT)})",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=15,
        help="Frames per second (default 15 → 12s for 180 steps)",
    )
    parser.add_argument(
        "--no-mp4",
        action="store_true",
        help="Skip MP4 output (GIF only). Useful when ffmpeg is unavailable.",
    )
    args = parser.parse_args()

    dead, damp, lp, persister, subclone, meta = _load_trajectory(args.traj_dir)
    written = _render(
        dead,
        damp,
        lp,
        persister,
        subclone,
        meta,
        args.traj_dir,
        fps=args.fps,
        skip_mp4=args.no_mp4,
    )
    print(f"\nDone. Wrote {len(written)} file(s):")
    for p in written:
        print(f"  {p.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
