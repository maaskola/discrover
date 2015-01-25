[![Build Status](https://travis-ci.org/maaskola/discrover.svg?branch=master)](https://travis-ci.org/maaskola/discrover)

# Discrover
### Discriminative discovery of sequence motifs with hidden Markov models

Copyright 2011, Jonas Maaskola.
This is free software under the GPL version 3, or later.
See the file [COPYING](COPYING) for detailed conditions of distribution.




This software package also contains code coming from an R library:<br>
Mathlib : A C Library of Special Functions<br>
Copyright (C) 2005-6 Morten Welinder <terra@gnome.org><br>
Copyright (C) 2005-10 The R Foundation<br>
Copyright (C) 2006-10 The R Core Development Team<br>
Mathlib is free software and distributed under the GNU GPL version 2.<br>
This software package uses routines from Mathlib to compute Chi-Square distribution probabilities by means of the incomplete Gamma function.


## Publication

The tool is described in the following [open-access article](http://nar.oxfordjournals.org/content/42/21/12995.full):<br/>
Jonas Maaskola and Nikolaus Rajewsky.
Binding site discovery from nucleic acid sequences by discriminative learning of hidden Markov models<br/>
Nucleic Acid Research, 42(21):12995-13011, Dec 2014. doi:10.1093/nar/gku1083


## Documentation

The package includes UNIX man pages for the programs.

The sub-directory doc contains a manual for this package, written in LaTeX.
A [PDF version](doc/discrover-manual.pdf) of the manual will be generated during the build process of this package.


## Galaxy front-end

There's a module in development to use Discrover inside the bioinformatics web framework Galaxy.
You can find it [here](https://github.com/maaskola/discrover-galaxy).


## Obtaining Discover

[Binary packages of Discrover](https://github.com/maaskola/discrover/releases) are available for select Linux distributions.
Notice that for Ubuntu a [PPA](https://launchpad.net/~maaskola/+archive/ubuntu/discrover) has been set up that can also be used for installing Discrover.
Instructions on how to [install](INSTALL.md) the packages (and how to use the PPA) or how to [manually build](BUILDING.md) Discrover are available in separate files.


## Sequence data

The synthetic sequence data used in the [publication](http://nar.oxfordjournals.org/content/42/21/12995) for motif discover performance evaluation are available [here](http://dorina.mdc-berlin.de/public/rajewsky/discrover/).

## Usage

Below is a minimal description on how to use this package.
Please refer to the UNIX man pages, the [manual](doc/discrover-manual.pdf), and the command line help for more information.

The package contains two main programs: ```plasma``` and ```discrover```.

```plasma``` is used to find IUPAC regular expression type motifs, and ```discrover``` learns HMMs.
Both use discriminative objective functions.
If no seeds are specified for ```discrover```, ```plasma``` will be used to find seeds automatically.

Command line help is available with ```discrover -h``` or ```discrover --help``` and, similarly, ```plasma -h``` or ```plasma --help```.

Note that some infrequently used options are hidden by default, and may be shown with the verbose switch: ```discrover -hv```

Even more obscure options are available by adding the very verbose switch: ```discrover -hV```

### Sequence logos
Both ```plasma``` and ```discrover``` can generate sequence logos (and ```discrover``` does so by default).
The same sequence logo creation routines are also available in the separate program ```discrover-logo```.

### Sequence shuffling
When ```plasma``` and ```discrover``` are given just a single FASTA file for analysis, they will automatically shuffle the sequences to create control sequences.
You can use the same sequence shuffling routines via the separate program ```discrover-shuffle```.
