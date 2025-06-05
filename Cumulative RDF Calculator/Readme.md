The scripts in this folder compute cumulative radial distribution functions (RDF) from
extended XYZ trajectories.

`atomRDF.py` is the original implementation which operates on explicit atom index lists. A newer script `group_rdf.py` provides a more flexible interface allowing atom groups to be defined for each molecule. Groups can be selectively included or excluded and specific group-to-group interactions may be enabled or disabled via pair filters. The results can optionally be plotted individually.

All runtime parameters are read from a JSON configuration file (default
`rdf_config.json`). See the example below for the new options used by
`group_rdf.py`.

Pair filters can restrict which neighbour groups contribute to the RDF for each
reference group.

```json
{
  "filename": "trajectory.xyz",
  "Range": 8.0,
  "Limit": 0.05,
  "Type": ["All"],
  "groups": {
    "mol1": [1,2,3],
    "mol2": [4,5,6]
  },
  "include_groups": ["mol1"],
  "exclude_groups": [],
  "pair_include": {
    "mol1": ["mol2"]
  },
  "pair_exclude": {},
  "plot_individual": true
}
```
