import argparse
import json
import numpy as np
import re
from pathlib import Path
import matplotlib.pyplot as plt


def parse_xyz(path):
    """Parse an extended XYZ trajectory.

    Returns
    -------
    positions : list of (N, 3) arrays
        Cartesian coordinates for each time step.
    elements : list of lists
        Atomic symbols for each time step.
    box : ndarray shape (3,)
        Box lengths assuming an orthogonal cell.
    """
    positions = []
    elements = []
    box = None
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            natoms = int(line.strip())
            comment = fh.readline()
            if box is None:
                match = re.search(r'Lattice="([^"]+)"', comment)
                if not match:
                    raise ValueError("Lattice information missing in XYZ header")
                vals = [float(x) for x in match.group(1).split()]
                box = np.array([abs(vals[0]), abs(vals[4]), abs(vals[8])])
            step_pos = []
            step_ele = []
            for _ in range(natoms):
                atom_line = fh.readline()
                if not atom_line:
                    break
                parts = atom_line.split()
                step_ele.append(parts[1])
                step_pos.append([float(x) for x in parts[2:5]])
            if len(step_pos) != natoms:
                break
            positions.append(np.asarray(step_pos, dtype=float))
            elements.append(step_ele)
    if not positions:
        raise ValueError("No frames parsed from XYZ file")
    return positions, elements, box


def compute_rdf(
    positions,
    elements,
    box,
    groups,
    types,
    bins,
    include=None,
    exclude=None,
    pair_include=None,
    pair_exclude=None,
):
    """Calculate cumulative RDF for atom groups with optional pair filters."""

    include = set(include or groups.keys())
    exclude = set(exclude or [])
    pair_include = pair_include or {}
    pair_exclude = pair_exclude or {}

    # Pre-convert group indices to zero-based arrays for efficiency
    groups_idx = {g: np.array(idxs, dtype=int) - 1 for g, idxs in groups.items()}
    hist = {g: np.zeros(len(bins) - 1, dtype=int) for g in groups}

    for step_pos, step_ele in zip(positions, elements):
        step_ele = np.array(step_ele)
        for gname, idxs in groups_idx.items():
            if gname not in include or gname in exclude:
                continue

            sel_pos = step_pos[idxs]

            allowed = pair_include.get(gname)
            if allowed is None:
                allowed = list(groups_idx.keys())
            allowed = [a for a in allowed if a not in pair_exclude.get(gname, [])]

            # gather indices of allowed neighbour groups
            target_indices = []
            for aname in allowed:
                if aname not in groups_idx:
                    raise ValueError(f"Unknown group '{aname}' in pair filter")
                target_indices.append(groups_idx[aname])
            if not target_indices:
                continue
            target_indices = np.concatenate(target_indices)

            if "All" in types:
                mask = np.ones(len(target_indices), dtype=bool)
            else:
                mask = np.isin(step_ele[target_indices], types)
            target_pos = step_pos[target_indices][mask]

            for p in sel_pos:
                delta = np.abs(target_pos - p)
                delta = np.where(delta < 0.5 * box, delta, np.abs(delta - box))
                dist = np.sqrt((delta ** 2).sum(axis=1))
                hist[gname] += np.histogram(dist, bins=bins)[0]

    return hist


def save_results(hist, bins, prefix="group"):
    for name, counts in hist.items():
        with open(f"{prefix}_{name}.txt", "w") as fh:
            fh.write("Bin\tCount\n")
            for b, c in zip(bins[:-1], counts):
                fh.write(f"{b:.4f}\t{c}\n")


def plot_results(hist, bins, plot_individual=False):
    if plot_individual:
        for name, counts in hist.items():
            plt.figure()
            plt.plot(bins[:-1], np.cumsum(counts))
            plt.xlabel("Distance")
            plt.ylabel("Cumulative count")
            plt.title(f"RDF for {name}")
            plt.tight_layout()
            plt.savefig(f"rdf_{name}.png")
            plt.close()
    else:
        plt.figure()
        for name, counts in hist.items():
            plt.plot(bins[:-1], np.cumsum(counts), label=name)
        plt.xlabel("Distance")
        plt.ylabel("Cumulative count")
        plt.legend()
        plt.tight_layout()
        plt.savefig("rdf_combined.png")
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="Compute cumulative RDF for atom groups")
    parser.add_argument("config", nargs="?", default="rdf_config.json", help="Path to config JSON")
    args = parser.parse_args()
    with open(args.config) as fh:
        cfg = json.load(fh)
    filename = cfg["filename"]
    groups = cfg["groups"]
    bins = np.arange(cfg["Limit"], cfg["Range"] + cfg["Limit"], cfg["Limit"])
    types = cfg.get("Type", ["All"])
    include = cfg.get("include_groups")
    exclude = cfg.get("exclude_groups")
    pair_in = cfg.get("pair_include")
    pair_ex = cfg.get("pair_exclude")
    plot_individual = cfg.get("plot_individual", False)

    positions, elements, box = parse_xyz(filename)
    hist = compute_rdf(
        positions,
        elements,
        box,
        groups,
        types,
        bins,
        include,
        exclude,
        pair_in,
        pair_ex,
    )
    save_results(hist, bins, prefix="rdf")
    plot_results(hist, bins, plot_individual)


if __name__ == "__main__":
    main()
