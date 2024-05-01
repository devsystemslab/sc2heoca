## Installation

sc2heoca runs on Linux and Windows 10 using the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10). 
Mac OSX should work as well, but we don't use it ourselves, so unexpected issues might pop up. 
If you have any issues let us know, so we can try to fix it!

### Creat conda environment

The recommended way to install sc2heoca is by using [conda](https://docs.continuum.io/anaconda). 
Activate the [bioconda](https://bioconda.github.io/) channel if you haven't used bioconda before.
You only have to do this once.

``` bash
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

You can install sc2heoca by creating a specific environment.
You only have to do this once.

``` bash
$ conda create -n sc2heoca python==3.9.16
```

Don't forget to activate the environment whenever you want to use sc2heoca.

``` bash
# Activate the environment before you use sc2heoca
$ conda activate sc2heoca
```

### Install `sc2heoca` from GitHub

You can also install sc2heoca with `pip`. 

``` bash
$ pip install git+https://github.com/devsystemslab/sc2heoca
``` 
