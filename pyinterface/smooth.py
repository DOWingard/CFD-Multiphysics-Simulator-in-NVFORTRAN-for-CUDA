#!/usr/bin/env python3

import os
import glob
import numpy as np
from collections import deque
from scipy.ndimage import label

def smooth_nan_3d(data, mesh_size):
    """
    Smooth NaNs in a 3D grid using the average of 26 neighbors.
    If neighbors are NaN, include their neighbors recursively.
    Returns the smoothed data and number of NaNs/NaN clusters.
    """
    pressure_grid = np.full((mesh_size, mesh_size, mesh_size), np.nan, dtype=float)
    for row in data:
        i, j, k, p = int(row[0])-1, int(row[1])-1, int(row[2])-1, row[3]
        pressure_grid[i, j, k] = p

    nan_indices = np.argwhere(np.isnan(pressure_grid))
    total_nans = len(nan_indices)
    
    # Count NaN clusters using 26-cell containing domain
    if total_nans > 0:
        binary_grid = np.isnan(pressure_grid).astype(int)
        structure = np.ones((3,3,3), dtype=int)
        _, num_clusters = label(binary_grid, structure=structure)
    else:
        num_clusters = 0

    neighbor_offsets = [(di,dj,dk) 
                        for di in [-1,0,1] 
                        for dj in [-1,0,1] 
                        for dk in [-1,0,1] 
                        if not (di==0 and dj==0 and dk==0)]

    for idx in nan_indices:
        i, j, k = idx
        visited = set()
        values = []

        queue = deque([(i,j,k)])
        while queue:
            ci, cj, ck = queue.popleft()
            if (ci,cj,ck) in visited:
                continue
            visited.add((ci,cj,ck))

            for di,dj,dk in neighbor_offsets:
                ni, nj, nk = ci+di, cj+dj, ck+dk
                if 0 <= ni < mesh_size and 0 <= nj < mesh_size and 0 <= nk < mesh_size:
                    val = pressure_grid[ni,nj,nk]
                    if np.isnan(val):
                        if (ni,nj,nk) not in visited:
                            queue.append((ni,nj,nk))
                    else:
                        values.append(val)
        if values:
            pressure_grid[i,j,k] = np.mean(values)
        else:
            pressure_grid[i,j,k] = 10000.0  # fallback

    smoothed_data = []
    for k in range(mesh_size):
        for j in range(mesh_size):
            for i in range(mesh_size):
                smoothed_data.append([i+1, j+1, k+1, pressure_grid[i,j,k]])

    return np.array(smoothed_data), total_nans, num_clusters

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    data_dir = os.path.join(repo_root, "kernel", "data")

    files = sorted(glob.glob(os.path.join(data_dir, "fort.*")))
    if not files:
        print(f"No fort.* files found in {data_dir}")
        return

    first_data = np.loadtxt(files[0])
    if first_data.ndim == 1:
        mesh_size = int(first_data[:3].max())
    else:
        mesh_size = int(np.max(first_data[:, :3]))

    for f in files:
        data = np.loadtxt(f)
        if data.ndim == 1:
            data = data.reshape(1,-1)
        smoothed, nans, clusters = smooth_nan_3d(data, mesh_size)
        np.savetxt(f, smoothed, fmt="%d %d %d %.10e")
        print(f"Processed {f}: NaNs={nans}, NaN clusters={clusters}")

if __name__ == "__main__":
    main()
