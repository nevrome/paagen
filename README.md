[![GitHub Workflow Status](https://github.com/nevrome/paagen/actions/workflows/normalCheck.yml/badge.svg)](https://github.com/nevrome/paagen/actions/workflows/normalCheck.yml)

# paagen

The **P**oseidon **a**rtificial **a**ncestry **gen**erator (paagen) is an experimental CLI application written in [Haskell](https://www.haskell.org/) to generate artificial genotype data for (ancient) human individuals. It requires input data in the [Poseidon format](https://poseidon-framework.github.io/#/).

### Subcommands

- admixpops: Genotype profile generation based on admixture proportions
- spacetime: Genotype profile generation based on spatiotemporal position

### Install

To install the latest development version you can follow these steps:

1. Install the Haskell build tool [Stack](https://docs.haskellstack.org/en/stable/README/)
2. Clone the repository
3. Execute `stack install` inside the repository to build the tool and automatically copy the executables to `~/.local/bin` (which you may want to add to your path). This will install the compiler and all dependencies into folders that won't interfere with any installation you might already have.