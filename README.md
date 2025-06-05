# Trajectory Analysis

A small collection of Python utilities for analysing atomic trajectories from ab-initio or classical molecular dynamics. The tools operate on `.xyz` files exported from [Ovito](https://www.ovito.org/) with particle indices enabled.

The repository currently ships two standalone calculators:

* **Distance Calculator** – time evolution of distances between specific atom pairs.
* **Cumulative RDF Calculator** – cumulative radial distribution functions for selected atoms.

Each script reads its parameters from a JSON configuration file so that no code changes are required between runs.

## Requirements

* Python 3.8+
* [NumPy](https://numpy.org/)

Install the dependency via `pip install numpy` if it is not already available.

## Distance Calculator

> Location: `Distance Calculator/DistanceDistances.py`

This script computes the distance between multiple pairs of atoms for every timestep in the trajectory. The list of pairs and the input file are taken from `distance_config.json`:

```json
{
  "filename": "All.xyz",
  "AtmA": [1, 2],
  "AtmB": [3, 4],
  "multiply": 1
}
```

Run the program with:

```bash
python Distance\ Calculator/DistanceDistances.py distance_config.json
```

Results are written to `ORDERED_LIST.txt` with one column per pair.

## Cumulative RDF Calculator

> Location: `Cumulative RDF Calculator/atomRDF.py`

This tool calculates cumulative radial distribution functions over the entire trajectory. Runtime options are provided through `rdf_config.json`:

```json
{
  "filename": "Distances5.xyz",
  "AtmA": [1, 2],
  "AtmB": ["O"],
  "Range": 8.0,
  "Limit": 0.05,
  "Type": ["All", "O"],
  "multiply": 25
}
```

Execute with:

```bash
python Cumulative\ RDF\ Calculator/atomRDF.py rdf_config.json
```

The script writes `<index> <element>.txt` files containing the histogram for each `AtmA` index and `COMBINED_VALUES_<element>.txt` which aggregates the results.

## Preparing Trajectories with Ovito

Both calculators expect trajectories exported from Ovito with atom indices intact and an orthogonal simulation cell. Visualise your system in Ovito to determine the correct atom IDs for the `AtmA` and `AtmB` lists.

## Citation

If you use this software for published research, please cite it using the information in `CITATION.cff`.

## License

This project is distributed under the terms of the GNU General Public License v3. See [`LICENSE`](LICENSE) for details.
