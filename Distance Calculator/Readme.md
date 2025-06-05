# Distance Calculator

This script computes distances between selected pairs of atoms for every step of a trajectory. The input trajectory must be exported from Ovito in `.xyz` format with atom indices in the first column:

```
AtomID  Element  X  Y  Z
```

All runtime options are defined in `distance_config.json` (trajectory filename and lists of atom indices). Run the program as:

```bash
python DistanceDistances.py distance_config.json
```

The distances for each pair are written to `ORDERED_LIST.txt`.

### Performance note

Distance computations use NumPy for vectorised operations. A synthetic benchmark (1000 steps, 60 pairs) completed in **~0.005 s** on Python 3.12 versus **~0.49 s** for the previous nested-loop version.
