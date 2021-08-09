[![GitHub Workflow Status](https://github.com/nevrome/paagen/actions/workflows/normalCheck.yml/badge.svg)](https://github.com/nevrome/paagen/actions/workflows/normalCheck.yml)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/nevrome/paagen?include_prereleases)
![GitHub all releases](https://img.shields.io/github/downloads/nevrome/paagen/total)](https://github.com/nevrome/paagen/releases)

# paagen

The **P**oseidon **a**rtificial **a**ncestry **gen**erator (paagen) is an experimental CLI application written in [Haskell](https://www.haskell.org/) to generate artificial genotype data for (ancient) human individuals. It requires input data in the [Poseidon format](https://poseidon-framework.github.io/#/).

### Install

For stable release versions we automatically prepare binaries that can be downloaded and run.

You can download them here: [ [Linux 📥](https://github.com/nevrome/paagen/releases/latest/download/paagen-Linux) | [macOS 📥](https://github.com/nevrome/paagen/releases/latest/download/paagen-macOS) | [Windows 📥](https://github.com/nevrome/paagen/releases/latest/download/paagen-Windows.exe) ]. Older release versions are available [here](https://github.com/nevrome/paagen/releases).

### Subcommands

- admixpops: Genotype profile generation based on admixture proportions
- spacetime: Genotype profile generation based on spatiotemporal position

## For developers

To install the latest development version you can follow these steps:

1. Install the Haskell build tool [Stack](https://docs.haskellstack.org/en/stable/README/)
2. Clone the repository
3. Execute `stack install` inside the repository to build the tool and automatically copy the executables to `~/.local/bin` (which you may want to add to your path). This will install the compiler and all dependencies into folders that won't interfere with any installation you might already have.

### Preparing a new stable release

The Github Actions script in `.github/workflows/release.yml` registers a new draft release and automatically builds and uploads currycarbon binaries when a new Git tag with the prefix `v*` is pushed. 

```bash
# locally register a new tag (e.g. 0.3.1)
git tag -a v0.3.1 -m "see CHANGELOG.md"
# push tag
git push origin v0.3.1
```

In case of a failing build delete the tag and the release draft on Github and then delete the tag locally with

```bash
git tag -d v0.3.1
```

before rerunning the procedure above.
