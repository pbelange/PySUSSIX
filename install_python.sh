                                              

# installing conda
mkdir ./Executables
if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform
    if [ "$(uname -m)" == "x86_64" ]; then
        wget -O ./Executables/Miniforge3-latest.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh

    elif [ "$(uname -m)" == "arm64" ]; then
        wget -O ./Executables/Miniforge3-latest.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
    fi
elif [ "$(uname)" == "Linux" ]; then
    # Do something under Linux platform
    wget -O ./Executables/Miniforge3-latest.sh ./Executables/Miniconda3-latest.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
fi
bash ./Executables/Miniforge3-latest.sh -b  -p ./Executables/miniforge3 -f


# create your own virtual environment in a new folder
source ./Executables/miniforge3/bin/activate

conda update -n base -c conda-forge conda
conda create -n py-sussix python=3.10
conda activate py-sussix


# conda install -c conda-forge numpy gfortran
conda install -c conda-forge compilers

# makefile of newton.f
cd PySUSSIX/ducksussix
make
cd ../..


# Install generic python packages
#========================================
pip install jupyterlab
pip install ipywidgets
pip install pandas
# pip install ipython
pip install numpy
pip install matplotlib
pip install scipy
pip install ipympl
pip install ruamel.yaml
pip install rich
pip install lfm
pip install pynaff
pip install NAFFlib
pip install pyarrow

# Adding the jupyter kernel to the list of kernels
python -m ipykernel install --user --name py-sussix --display-name "py-sussix"
# ========================================

# Install project package
pip install -e ./


