# flutype-analysis
Analysis scripts for FluType microarray and microtiter data


## Installation
To run the jupyter notebook scripts setup a virtual environment and register
the respective kernel for the notebook via

```
mkvirtualenv flutype_analysis
(flutype_analysis) pip install -r requirements.txt
cd ~/git/flutype_analysis
(flutype_analysis) python -m ipykernel install --user --name=flutype-analysis
```

Installation of package for development
```
(flutype_analysis) cd ~/git/flutype_analysis
(flutype_analysis) pip install -e .
```

Run tox tests against different python versions:
```
sudo apt-get install tox
tox -e py27
tox -e py35
tox -e py36
```
