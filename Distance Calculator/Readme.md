These codes are used to analyze the AIMD/CMD trajectories stored in the .xyz format. This calculates distances between specified pairs of atoms.
All xyz files need to be tabulated in the following format which can be easily done using the Ovito Basic program.
------------------------------------------------------------------------------------
Atom ID    |    Element Name   |    X   |    Y   |    Z    |

The script now reads its input parameters from a JSON file named `distance_config.json`.
Provide the filename and atom index pairs there before running the program.

### Performance note
Distance computations now use NumPy array operations. A small synthetic
benchmark (1000 steps, 60 pairs) completed in **~0.005 s** versus **~0.49 s**
with the previous double loops on Python 3.12.
