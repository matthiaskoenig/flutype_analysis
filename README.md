[![Build Status](https://travis-ci.org/matthiaskoenig/flutype_analysis.svg?branch=develop)](https://travis-ci.org/matthiaskoenig/flutype_analysis)
[![License (LGPL version 3)](https://img.shields.io/badge/license-LGPLv3.0-blue.svg?style=flat-square)](http://opensource.org/licenses/LGPL-3.0)


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

&copy; 2017 Janek Grzegorzewski & Matthias König.