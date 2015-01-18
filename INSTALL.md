# Discrover
### Discriminative discovery of sequence motifs with hidden Markov models

Copyright 2011, Jonas Maaskola.
This is free software under the GPL version 3, or later.
See the file [COPYING](COPYING) for detailed conditions of distribution.

There are multiple ways how you can install Discrover.

1. We provide [binary packages](#binary-packages) for various Linux distributions.
   If your Linux distribution is supported, this is the easiest way to install Discrover.
2. For other Linux distributions we provide [source-based packages](#source-based-packages) scripts that automatically build Discrover for you.
3. Should none of the above steps work for, we also have instructions on [how to build](BUILDING.md) Discrover.

Note that below, we frequently write ``X.Y.Z`` to indicate the version number (and sometimes ``X.Y.Z-R`` for version number and package release number).
You have to replace it with whatever version of Discrover that you want to install.

# Binary packages
Binary packages of Discrover are [available](https://github.com/maaskola/discrover/releases) for select Linux distributions.
Generally, these packages specify dependent software that is required for full functionality of Discrover.
Your distribution's package management  will install dependencies alongside Discrover.

The packages were built for 64-bit architectures and are available for:

* Debian 7.7 (Wheezy)
* Fedora 20
* Ubuntu 14.04 (Trusty)

The binary packages were generated with the specifications found in [another repository](https://github.com/maaskola/discrover-packages).

Once downloaded, these packages can be installed as follows.

## Debian

Just use ``dpkg`` to install the package:

```sh
dpkg -i discrover-X.Y.Z_amd64_debian_wheezy.deb
```

## Fedora

You can install with ``yum``:

```sh
yum install discrover-X.Y.Z-R.fc20.x86_64.rpm
```

To update just use

```sh
yum update discrover-X.Y.Z-R.fc20.x86_64.rpm
```

## Ubuntu

A Ubuntu [PPA](https://launchpad.net/~maaskola/+archive/ubuntu/discrover) has been set up that can be used for installing Discrover.
You can add this PPA to your system and install Discrover with the following commands (you might have to replace `trusty` with your Ubuntu system's release name):

```sh
deb http://ppa.launchpad.net/maaskola/discrover/ubuntu trusty main
deb-src http://ppa.launchpad.net/maaskola/discrover/ubuntu trusty main
apt-get install discrover
```

Aside from using this PPA, you can also directly download and install a Discrover release package.
Ubuntu being a Debian-derivative, you can use ``dpkg`` to install the package:

```sh
dpkg -i discrover-X.Y.Z_amd64_ubuntu_trusty.deb
```

Note the following advantage of using the PPA: when future Discrover releases happen your system's Discrover installation will be automatically updated, while you manually have to retreive and install packages if you opt for the latter solution.
Further, the PPA also includes packages for 32-bit architectures, as well as older Ubuntu releases.
It is thus recommended to use the PPA-based installation for Ubuntu.


# Source-based packages
For Gentoo and Arch Linux, there are source-based packages available.
You can use these to automatically compile and install Discrover.

## Gentoo

Ebuilds for Discrover are available [here](https://github.com/maaskola/discrover-packages/tree/master/gentoo/), but are also included in the gentoo-science overlay.
You can install this overlay, accept unstable versions of Discrover, and then simply emerge discrover:

```sh
overlay -a science
echo "=sci-biology/discrover-X.Y.Z ~amd64" >> /etc/portage/package.keywords
emerge -av discrover
```

## Arch Linux

We have written a PKGBUILD file for Arch Linux that specifies dependencies and how to build Discrover.
You can [get it](https://aur.archlinux.org/packages/discrover/) via the [Arch User Repository (AUR)](https://aur.archlinux.org/).
The most current version is typically available [here](https://github.com/maaskola/discrover-packages/blob/master/arch/PKGBUILD).
