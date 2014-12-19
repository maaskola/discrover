# Discrover
### Discriminative discovery of sequence motifs with hidden Markov models

Copyright 2011, Jonas Maaskola.
This is free software under the GPL version 3, or later.
See the file [COPYING](COPYING) for detailed conditions of distribution.

# Manual building

## Dependencies for building Discrover

To build Discrover you need to have some software.
During the configuration phase it is checked whether these packages are found on your system.
Some functionality of Discrover is enabled if you have certain, optional dependencies.
The required and optional dependencies include:


## A C++11 supporting compiler
Discrover is written in C++11, so it is necessary to use an up-to-date version of your compiler.
The GNU compiler collection supports all necessary features to compile this project as of version 4.6.


## CMake
We use [CMake](http://www.cmake.org/) to construct Makefiles for building Discrover.


### Boost C++ libraries
We make use of code in the [Boost C++ libraries](http://www.boost.org/), so please install them.
These libraries are available under the [Boost License](http://www.boost.org/users/license.html).
The version of Boost must be recent enough to include the V3 filesystem library.
Version 1.48 and more recent versions are known to work.


### OpenMP
We use OpenMP to support parallelization.
As of version 4.2 the GCC supports OpenMP out of the box, so we suggest building with a recent compiler version.


### R library
We are using code from the R library to compute the logarithm of the chi-square distribution function.
During configuration it is checked if the R library is found.
If that is the case then we link to it, otherwise code extracted from it is built and used.


### DREME from the MEME suite
Optionally, the program DREME from the [MEME suite](http://meme.nbcr.net/meme/) of motif analysis tools can be integrated for seeding of HMM motifs.
During configuration it is checked if any programs named ```dreme``` or ```meme-dreme``` can be found.
If that is the case, you will be able to use DREME for motif seeding by specifying ```--algo dreme```.


### LaTeX
LaTeX is another optional dependency.
In particular, ```pdflatex``` is used to compile the manual.
During configuration it is checked if your system has ```pdflatex```.
If it has, then the manual is built, otherwise building of the manual is skipped.

Note that some LaTeX-packages are used that may not be part of your core LaTeX installation.
Unfortunately, while our build-system will detect the presence of ```pdflatex```, it is currently not sophisticated enough to figure out if these LaTeX packages are present.
Thus, if later during compilation you see errors while building the manual, please ensure that you have the requisite LaTeX packages available on your system.
For example, on Debian, they are part of the ```texlive-latex-extra``` package.
Similarly, on Gentoo, they are part of the ```dev-texlive/texlive-latexextra``` package.


## Installing build-time dependencies
Note that the above-given list of dependencies is required only for BUILDING, not for running.

### Debian and Ubuntu
On Debian and Ubuntu, you can install all necessary and optional software to build Discrover with the following command:

```sh
apt-get install git cmake g++ imagemagick libboost-all-dev texlive texlive-latex-base latex-xcolor texlive-latex-extra pgf ruby ruby-dev
```

### Gentoo
Similarly, on Gentoo you can use:

```sh
emerge -av dev-vcs/git dev-util/cmake sys-devel/gcc media-gfx/imagemagick dev-libs/boost dev-texlive/texlive-latexextra dev-lang/ruby
```

### Arch
The corresponding command for Arch linux:

```sh
pacman -S git cmake make gcc boost texlive-core texlive-latexextra imagemagick ruby
```

### Fedora
On Fedora 20 the following command will install all dependencies required for building:

```sh
yum install gcc-c++ cmake git ImageMagick boost boost-devel texlive-latex-bin texlive-pgf texlive-xcolor texlive-collection-latexextra ruby ruby-devel
```

### Mac OS X
On Mac OS X, git is provided with XCode.
Using (brew)[http://brew.sh], you can install CMake and ImageMagick like this:

```sh
brew install cmake imagemagick
```

While brew also provides a binary package for GCC, we cannot use it, as it does not support OpenMP.
For this reason you need to rebuild GCC on your system:

```sh
brew install gcc --without-multilib
```

The boost brew package can also not be used because it relies on the GCC brew package and does not work with the manualy built GCC.
Therefore, please build boost according to [these instructions](http://qiita.com/misho/items/0c0b3ca25bb8f62aa681).
Finally, adapt the ``SET(BOOST_ROOT "...")`` statement in ``CMakeLists.txt`` to point to the place where you installed boost.

For TeX, please install [MacTex](https://tug.org/mactex/).

## Building

The code contained in this package is built in four steps.

1. [Edit the CMake build script (optional)](#step1)
2. [Execute the CMake build script](#step2)
3. [Compile the source code](#step3)
4. [Link and install the binary file](#step4)



### <a name="step1"></a>Step 1: Edit the CMake build script (optional)

The file [CMakeLists.txt](CMakeLists.txt) in the root directory of the source code tree contains the instructions to configure and build Discrover.
Among other things, you can set the installation target directory in it.
By default, it will be installed into the directory rooted at ```/usr/local```.
To change the installation target directory, find the following lines, and modify them to your liking.

```cmake
SET(LOCAL_PREFIX "/usr/local")
SET(CMAKE_INSTALL_PREFIX ${LOCAL_PREFIX})
SET(CMAKE_PREFIX_PATH ${LOCAL_PREFIX})
```

Explanation:
The variable ```CMAKE_INSTALL_PREFIX``` determines where the software will be installed after building.
In particular, programs will be installed in ```${CMAKE_INSTALL_PREFIX}/bin```, libraries in ```${CMAKE_INSTALL_PREFIX}/lib```, and documentation into ```${CMAKE_INSTALL_PREFIX}/share/doc/discrover```.
The variable ```CMAKE_PREFIX_PATH``` has to be set such that the required dependencies may be found.
In particular, the headers of Boost need to be found in ```${CMAKE_PREFIX_PATH}/include```, if not installed in a standard system directory like ```/usr/include```.
Similarly, the directory ```${CMAKE_PREFIX_PATH}/lib``` tells CMake where to find the Boost libraries if they are not installed in a standard system directory.

By default, both ```CMAKE_INSTALL_PREFIX``` and ```CMAKE_PREFIX_PATH``` are constructed from ```LOCAL_PREFIX```.
Hence, setting ```LOCAL_PREFIX``` to match your system might be the only change required in most cases, if necessary at all.


After this variable has been adapted to your system, you may proceed with the next step.





### <a name="step2"></a> Step 2: Execute the CMake build script

Change to the root directory of the package and execute

```sh
cmake .
```

Alternatively, you can execute

```sh
cmake -DCMAKE_INSTALL_PREFIX:PATH=/desired/installation/path .
```

where you would replace ```/desired/installation/path``` by the path to which you to install the package.
Note that you would have to issue this every time you configure the package, so it may be preferable to set this permanently as described in [step 1](#step1) above.


Explanation:
This will search for the paths to the required headers and libraries, create a directory called build and prepare everything for the subsequent compilation.
It will also check if your compiler supports the required features (C++11 and OpenMP support).

If anything fails at this step please have a look at [CMakeLists.txt](CMakeLists.txt) and see if some of the commented-out statements may help you.
If difficulties persists, please contact the author of this software.





### <a name="step3"></a> Step 3: Compile the source code


Execute from the package's root directory

```sh
make
```

Explanation:
This will compile the source code.
You may make use of parallel building by running

```sh
make -j N
```

where ```N``` is the number of CPUs that you want to use.




### <a name="step4"></a> Step 4: Link and install the binary file


While still in the root directory of this package, execute

```sh
make install
```

Explanation:
This will copy the libraries and binaries into the default installation path.
If you did not specify otherwise via [CMakeLists.txt](CMakeLists.txt) in [step 1](#step1), or via the command line in [step 2](#step2), then the installation will go to ```/usr/local/bin```, ```/usr/local/lib```, and ```/usr/local/share/doc```.
If these locations are part of your ```$PATH``` and ```$LD_LIBRARY_PATH``` environment variables you can then simply run the ```discrover``` program from anywhere.

Otherwise, you might still have to add the directories you installed to ```$PATH``` and ```$LD_LIBRARY_PATH```.
This can be done with a command like

```sh
export PATH=HERE/bin:$PATH
export LD_LIBRARY_PATH=HERE/lib:$LD_LIBRARY_PATH
```

where ```HERE``` will have to be the path that you installed to.
You might consider putting these commands into your ```~/.bashrc``` file or some place similar such that they are executed every time you log into your machine.
