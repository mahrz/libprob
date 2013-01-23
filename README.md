libprob
=======

Libprob is a discrete probability distribution library for C++11 based on [Eigen](http://eigen.tuxfamily.org/). The library uses variadic templates to encapsulate the matrix storage of multivariate probability distributions. Additional support for basic algebraic operations like joining, marginalization or Bayes' rules is provided, as well as information theoretic measures like entropy and mutual information.

Installation
------------

You will need to have [Eigen](http://eigen.tuxfamily.org/) and gcc 4.7 (std=c++11 is needed) installed on your computer to use the library. Cmake, gtest as well as doxygen are required to build the tests and documentation. Check out the repository or download the repository as a zip and extract it to a directory of your preference. Then download [gtest-1.6](http://code.google.com/p/googletest/downloads/list) and unpack the `gtest-1.6` folder to the libprob folder. Lets suppose the repository resides in your home directory in a folder named `libprob`, then create a build directory on the same level as the repository and execute the following commands:

```sh
~> mkdir libprob_build
~> cd libprob_build
~/libprob_build> cmake ../libprob
~/libprob_build> make test
~/libprob_build> make doc
~/libprob_build> make install
```

This will compile and execute the tests, build the documentation in the repository folder and installs the header only library into your include directory (which can be set via `CMAKE_INSTALL_PREFIX`). If you want to use the experimental information decomposition features (which are not tested yet and possibly broken) consult the documentation. These features require the [nlopt](http://ab-initio.mit.edu/wiki/index.php/NLopt) optimization library to be present on your system.

Examples
--------

Here is a simple example of the look & feel of the library, for more examples have a look at the testcases and the documentation.

```c++
using namespace prob;

RVAR_STATIC(W, 5)
RVAR_STATIC(A, 3)
RVAR_STATIC(S, 2)

[...]

distribution<double, W, A, given, S, W> p;

W w(3), wn(1);
A a(2);
S s(0);

init::random(p,gen);

p(w, a | s, wn) = 0.5;

p.normalize();

p.each_index([] (W w, A a, prob::given, S s, W wp)
{ 
  std::cout << "Index " << w << " " << a << "|"
    	       	      	<< s << " " << wp << std::endl;
});

```

Documentation
-------------

The documentation can be built locally using doxygen, for a documentation of the latest release look [here](http://mahrz.github.com/libprob/docs/).

