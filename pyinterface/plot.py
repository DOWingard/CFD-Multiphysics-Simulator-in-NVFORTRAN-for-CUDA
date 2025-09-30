#!/usr/bin/env python3

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio.v2 as imageio

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    data_dir = os.path.join(repo_root, "kernel", "data")
    plot_dir = os.path.join(repo_root, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(data_dir, "fort.*")))
    if not files:
        print(f"No fort.* files found in {data_dir}")
        return

    print(f"Found {len(files)} files in {data_dir}")
    frame_paths = []

    all_data = []
    for f in files:
        try:
            data = np.loadtxt(f)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] >= 4:
                all_data.append(data[:, :4])
        except:
            continue
    all_data = np.vstack(all_data)
    i_all, j_all, k_all, _ = all_data.T
    i_min, i_max = i_all.min(), i_all.max()
    k_min, k_max = k_all.min(), k_all.max()

    for idx, f in enumerate(files):
        try:
            data = np.loadtxt(f)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] < 4:
                print(f"Skipping {f}: expected 4 columns (i j k pressure)")
                continue
        except Exception as e:
            print(f"Failed to read {f}: {e}")
            continue

        i, j, k, p = data.T

        mask = (j >= 20) & (j <= 40)
        i, j, k, p = i[mask], j[mask], k[mask], p[mask]
        if len(p) == 0:
            continue

        p_mean = np.mean(p)
        p_std = np.std(p)
        if p_std == 0:
            p_std = 1e-12  
        vmin = p_mean - 2 * p_std  
        vmax = p_mean + 2 * p_std

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        scatter = ax.scatter(i, j, k, c=p, cmap="seismic", s=20, alpha=0.7,
                             vmin=vmin, vmax=vmax)  
        fig.colorbar(scatter, ax=ax, label="Pressure deviation from mean")

        ax.set_xlabel("i")
        ax.set_ylabel("j")
        ax.set_zlabel("k")
        ax.set_title(f"3D Pressure Distribution (20 ≤ j ≤ 40)\n{os.path.basename(f)}")

        ax.set_xlim(i_min, i_max)
        ax.set_ylim(20, 40)
        ax.set_zlim(k_min, k_max)

        plt.tight_layout()

        frame_path = os.path.join(plot_dir, f"frame_{idx:04d}.png")
        plt.savefig(frame_path, dpi=150)
        plt.close(fig)
        frame_paths.append(frame_path)

    gif_path = os.path.join(plot_dir, "pressure_evolution.gif")
    with imageio.get_writer(gif_path, mode="I", duration=0.5) as writer:
        for frame in frame_paths:
            writer.append_data(imageio.imread(frame))

    print(f"Saved in plots/")

if __name__ == "__main__":
    main()
