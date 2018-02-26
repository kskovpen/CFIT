CFIT is a tool to perform a template fit with including
systematic uncertainties in a correlation matrix. Systematic
correlations between templates are taken into account. The tool also
allows one to estimate statistical and systematic uncertainties
on the fit parameters and to perform the measurement of the scale factors.

```
To build a python wrapper:

git clone https://github.com/swig/swig.git
cd swig
./autogen.sh && ./configure  --prefix=<your build folder> && make && make install

./buildPyLib.sh

export LD_LIBRARY_PATH=/path/to/PyCFIT/:$LD_LIBRARY_PATH

Make sure that you specify the correct path to your python installation.
The PyCFIT library is _cfit.so. Documention of the methods is in cfit.py
```

Documentation is [available](./doc/doc.pdf)

