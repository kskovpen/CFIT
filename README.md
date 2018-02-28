CFIT is a tool to perform a template fit with including
systematic uncertainties in a correlation matrix. Systematic
correlations between templates are taken into account. The tool also
allows one to estimate statistical and systematic uncertainties
on the fit parameters and to perform the measurement of the scale factors.

```
To build a python wrapper:

# Build swig:
git clone https://github.com/swig/swig.git
cd swig
./autogen.sh && ./configure --prefix=<your build folder> && make && make install

# Add swig to path
export PATH=/path/to/swig_build_dir/bin:$PATH

# Build wrapper
./buildPyLib.sh

# Add library to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/PyCFIT/:$LD_LIBRARY_PATH

The PyCFIT library is _cfit.so. Documentation of the methods is in cfit.py
Pythonic examples are included in examples/ 
	test.py
	testSF.py
```

Documentation is [available](./doc/doc.pdf)

