
# discrover - discriminative discovery of sequence motifs with hidden Markov models

## Binding site pattern discovery from nucleic acid sequences by discriminative learning of hidden Markov models

Copyright 2011, Jonas Maaskola


This is free software under the GPL version 3, or later. See the file [COPYING](COPYING)
for detailed conditions of distribution.




This software packages also contains code coming from an R library:
Mathlib : A C Library of Special Functions
Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
Copyright (C) 2005-10 The R Foundation
Copyright (C) 2006-10 The R Core Development Team
Mathlib is free software and distributed under the GNU GPL version 2.
This software package uses routines from Mathlib to compute Chi-Square
distribution probabilities by means of the incomplete Gamma function.





## Documentation

The sub-directory doc contains a manual for this package, written in LaTeX.
A [PDF version](doc/discrover-manual.pdf) of the manual will be generated during the build process of this
package. 





## Dependencies

This software depends on other software. During the configuration phase it is
checked whether these packages are found on your system.
The dependencies are:


### A C++11 supporting compiler
This software package is written in C++11, so it is necessary to use an
up-to-date version of your compiler. The GNU compiler collection supports all
necessary features to compile this project as of version 4.6.


### CMake
We use [CMake](http://www.cmake.org/) to construct Makefiles to build discrover.


### Boost C++ libraries
We make use of code in the [Boost C++ libraries](http://www.boost.org/), so please install them.
These libraries are available under the [Boost License](http://www.boost.org/users/license.html).
The version of Boost must be recent enough to include the V3 filesystem library.
Version 1.48 and more recent versions are known to work.


### OpenMP
We use OpenMP to support parallelization. As of version 4.2 the GCC supports
OpenMP out of the box, so we suggest building with a recent compiler version.


### R library
We are using code from the R library to compute the logarithm of the chi-square
distribution function.  During configuration it is checked if the R library is
found. If that is the case then we link to it, otherwise code extracted from it
is built and used.


### DREME from the MEME suite
Optionally, the program DREME from the [MEME suite](http://meme.nbcr.net/meme/) of motif analysis tools can be
integrated for seeding of HMM motifs. During configuration it is checked if any
programs named ```dreme``` or ```meme-dreme``` can be found. If that is the case, you
will be able to use DREME for motif seeding by specifying ```--algo dreme```.




## Building

The code contained in this package is built in four steps.

1. [Edit the CMake build script (optional)](#step1)
2. [Execute the CMake build script](#step2)
3. [Compile the source code](#step3)
4. [Link and install the binary file](#step4)



### <a name="step1"></a>Step 1: Edit the CMake build script (optional)

The file [CMakeLists.txt](CMakeLists.txt) in the root directory of the source code tree contains
the instructions to configure and build discrover.
Among other things, you can set the target installation directory in it.
To do this, find the following lines, and modify them to your liking.

    SET(LOCAL_PREFIX "~/local")
    SET(CMAKE_INSTALL_PREFIX ${LOCAL_PREFIX})
    SET(CMAKE_PREFIX_PATH ${LOCAL_PREFIX})

Explanation:

The variable ```CMAKE_INSTALL_PREFIX``` determines where the software will be
installed after building. In particular, programs will be installed in
```${CMAKE_INSTALL_PREFIX}/bin``` and libraries in ```${CMAKE_INSTALL_PREFIX}/lib```.
The variable ```CMAKE_PREFIX_PATH``` has to be set such that the required dependencies
may be found. In particular, the headers of Boost need to be found in
```${CMAKE_PREFIX_PATH}/include```, if not installed in a standard system directory
like ```/usr/include```.
Similarly, the directory ```${CMAKE_PREFIX_PATH}/lib``` tells CMake where to find
the Boost libraries if they are not installed in a standard system directory.

By default, both ```CMAKE_INSTALL_PREFIX``` and ```CMAKE_PREFIX_PATH``` are constructed from
```LOCAL_PREFIX```, hence setting only ```LOCAL_PREFIX``` to match your system is the only
change required in most cases.


After this variable has been adapted to your system, you may proceed with
the next step.





### <a name="step2"></a> Step 2: Execute the CMake build script

Change to the root directory of the package and execute

```cmake .```

Alternatively, you can execute

```cmake -DCMAKE_INSTALL_PREFIX:PATH=/desired/installation/path .```

where you would replace ```/desired/installation/path``` by the path to which you
to install the package. Note that you would have to issue this every time you
configure the package, so it may be preferable to set this permanently as
described in [step 1](#step1) above.


Explanation:
This will search for the paths to the required headers and libraries, create a
directory called build and prepare everything for the subsequent compilation.
It will also check if your compiler supports the required features (C++11
and OpenMP support).

If anything fails at this step please have a look at [CMakeLists.txt](CMakeLists.txt) and see if
some of the commented-out statements may help you. If difficulties persists,
please contact the author of this software.





### <a name="step3"></a> Step 3: Compile the source code


Execute from the package's root directory

```make```


Explanation:
This will compile the source code. You may make use of parallel building by
running

```make -j N```

where N is the number of CPUs that you want to use.




### <a name="step4"></a> Step 4: Link and install the binary file


While still in the root directory of this package, execute

```make install```


Explanation:
This will copy the libraries and binaries into the default installation path,
which is ```/usr/local/bin```, ```/usr/local/lib```, and ```/usr/local/share/doc``` if you did
not specify otherwise via [CMakeLists.txt](CMakeLists.txt) in [step 1](#step1), or via the command line in
[step 2](#step2).

If these locations are part of your ```$PATH``` and ```$LD_LIBRARY_PATH``` environment
variables you can then simply run the ```discrover``` program from anywhere.

Otherwise, you might still have to add the directories you installed to to ```$PATH```
and ```$LD_LIBRARY_PATH```. This can be done with a command like

    export PATH=HERE/bin:$PATH
    export LD_LIBRARY_PATH=HERE/lib:$LD_LIBRARY_PATH

where ```HERE``` will have to be the path that you installed to.
You might consider putting these commands into your ```~/.bashrc``` file or some place
similar such that they are executed every time you log into your machine.






## How to run

Below is a minimal description on how use this package. Please refer to the [manual](doc/discrover-manual.pdf) for more information.

The package contains two binaries, ```plasma``` and ```discrover```.

```plasma``` is used to find IUPAC regular expression type motifs, and ```discrover``` learns
HMMs. Both use discriminative objective functions.
If no seeds are specified for ```discrover```, ```plasma``` will be used to find seeds automatically.

Please refer to the command line help, which is available after installing with

```discrover -h```
  or
```discrover --help```

Similarly:

```plasma -h```
  or
```plasma --help```


Note that some infrequently used options are hidden by default, and may be shown
with the verbose switch:

```discrover -hv```

Even more obscure options are available by adding the very verbose switch:

```discrover -hV```

