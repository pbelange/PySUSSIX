# PySUSSIX

A Python wrapper for the frequency analysis tool `SUSSIX`, cf. https://cds.cern.ch/record/702438/ .

## Installation on Linux

`PySussix` itself is only a wrapper around a Fortran implementation of `SUSSIX`, which first needs to be built.
For this, first clone the repository with:
```bash
git clone https://github.com/PyCOMPLETE/PySUSSIX
```

In order to build the package, ensure you have access to `numpy`, `gcc` and potentially `gfortran` in your desired Python environment.
In a `conda` environment, this can be done with:
```bash
conda install -c conda-forge numpy gcc gfortran
```

Then, run the `install` installation script:
```bash
bash install
```

This will move into the `src` directory and run the build command.
If you experience problems with compiling the Fortran code, make sure you have access to `gfortran` and try to change the `Makefile`'s compiler flag to `--fcompiler=gfortran` (instead of the default `--fcompiler=gnu95`).

The compilation is successful if it produced a `pysussix.cpython-<python_version>-x86_64-linux-gnu.so` "shared object" (`.so`) file, where `<python_version>` here would be `38` if your environment's Python version is `3.8.x`.

In order to guarantee the package is in your `PYTHONPATH`, the simplest solution is to move or copy the `PySUSSIX` folder with this build into your `site-packages`.
One can find the location of this folder for the currently activated Python environment by looking for the `USER_SITE` entry in the output produced by running:
```bash
python -m site
```

Copy the `PySUSSIX` folder in there with the built shared object file into your `site-packages`:
```bash
cd ..
cp -rf PySUSSIX /path/to/site-packages
```

You're done!

## Usage

Once the aforementioned steps have been followed, one can simply import the module with:
```python
from PySUSSIX import PySussix
```

The `Sussix` class can also be directly imported and used with:
```python
from PySUSSIX.PySussix import Sussix

Sussix.sussix_inp(...)
Sussix.sussix(...)
```

#### Note: Compatibility Across Environments

Due to the build and linking process, it is most likely that the shared object file compiled with a specific Python version is incompatible with another version.
For instance, the example file above built with Python `3.8` would most likely not work if installed in the same way in a Python `3.6` environment.

This means the safest way to get `PySUSSIX` available in several environments is to make a build for each one.
Thankfully it is quite quick to compile.