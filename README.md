# flutype-analysis
Analysis scripts for FluType microarray and microtiter data


## Installation
To run the jupyter notebook scripts setup a virtual environment and register
the respective kernel for the notebook via

```
mkvirtualenv flutype-analysis
pip install -r requirements.txt
cd ~/git/flutype-analysis
python -m ipykernel install --user --name=flutype-analysis
```

Installation of package for development
```
pip install -e .
```