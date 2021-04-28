# paagen: Poseidon artificial ancestry generator

This is a highly experimental proof of concept CLI application to generate artificial genotype data for ancient human individuals. It requires input data in the [Poseidon format](https://poseidon-framework.github.io/#/). A jupyter notebook in `/playground` illustrates an example run, but that requires [jupyter notebook](https://jupyter.org/), an [R kernel](https://github.com/IRkernel/IRkernel) for it, [trident](https://poseidon-framework.github.io/#/trident), [plink1.9](https://www.cog-genomics.org/plink), and some R packages.

### Install

To instead install the latest development version you can follow these steps:

1. Install the Haskell build tool [Stack](https://docs.haskellstack.org/en/stable/README/)
2. Clone the repository
3. Execute `stack install` inside the repository to build the tool and automatically copy the executables to `~/.local/bin` (which you may want to add to your path). This will install the compiler and all dependencies into folders that won't interfere with any installation you might already have.

### Commands

#### spacetime

A nondeterministic generator based on spatiotemporal position.

<p align="left">
  <img src="figures/spacetime_schema.png" width = 400>
</p>

```
Usage: paagen spacetime (-d|--baseDir DIR) [-p|--positionString ARG] 
                        [--positionFile ARG] [--neighbors ARG] 
                        [--temporalDistanceScaling ARG] [--outFormat ARG]
                        (-o|--outPath ARG)
  Genotype profile generation based on spatiotemporal position

Available options:
  -h,--help                Show this help text
  -d,--baseDir DIR         A base directory to search for Poseidon Packages
                           (could be a Poseidon repository)
  -p,--positionString ARG  Spatiotemporal positions of interest: each position
                           is a string of the form
                           "[id:group](time,latitude,longitude)". Multiple
                           positions can be listed separated by ;. id and group
                           are simple strings, latitude and longitude must be
                           given in decimal degrees and time in calBC/AD, so
                           3245BC would be -3245 and 1148AD just 1148
  --positionFile ARG       A file with a list of spatiotemporal positions. Works
                           just as -p, but multiple values can also be separated
                           by newline, not just by ;. -p and --positionFile can
                           be combined.
  --neighbors ARG          Number of nearest neighbors to consider for the
                           calculation (default: 50)
  --temporalDistanceScaling ARG
                           Scaling factor for temporal distance in relation to
                           spatial distance. Default is 1, so 1year =
                           1km. (default: 1.0)
  --outFormat ARG          The format of the output genotype data: EIGENSTRAT or
                           PLINK
  -o,--outPath ARG         The output directory path

```
