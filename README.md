# ORCA-tools

Our collection of python tools for [ORCA](https://www.kofo.mpg.de/de/forschung/services/orca).

## Installation

Clone this repository and run `pip install .` inside the main directory. If you want to always use the latest content of the repo you can use the 'developement' install of pip by running `pip install -e .`. Just doing `git pull` to get the latest content of the repo will then automatically result in the usage of the latest code without need to reinstall.

## Usage

You can either run `python3 -m orcatools` to use the command line interface, or import the needed features in your own scripts.

E.g. in a Jupyter notebook you can do

```python
from orcatools.go import plot_orca_go
path = "path to your ORCA calculation"
plot_orca_go(filename="myjobname", presentation=True, path=path, show=True)
```

to also see your plots inside the notebook.
