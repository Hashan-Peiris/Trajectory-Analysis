# Cumulative RDF Calculator

This program computes the cumulative radial distribution function for a set of reference atoms across the entire trajectory. The trajectory must be provided in Ovito-generated `.xyz` format with atom indices.

Configuration parameters such as the trajectory file, atom indices and bin settings are read from `rdf_config.json`.

Run the script with:

```bash
python atomRDF.py rdf_config.json
```

Output files `<index> <element>.txt` contain the histogram for each reference atom and `COMBINED_VALUES_<element>.txt` aggregates them for easy plotting.
