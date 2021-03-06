# Examples of "plugins" for [fwdpy](http://molpopgen.github.io/fwdpy/).

This repo contains examples of extending fwdpy using custom Cython code.

The gist of the trick is to:

1. make a Cython source file (foo.pyx) defining your cool new function.
2. Make foo.pyxbld that tells pyximport that this will be a C++ module requiring the flag -std=c++11.  Note: the contents of this file will be the same for any fwdpy extension.  The rule, though, is that module.pyx must have a module.pyxbld associated with it.
3. We use pximport to build our module on the fly and import it into our script using it.  This is analgous to Rcpp's "cppfunction".

For this to work:

1. fwdpy must be installed.
2. Cython must be installed.
3. cythongsl must be installed.

The latter two may be "pip installed".

To compile a module, you need to know where the fwdpy headers are.  This can be tricky, and I'm trying to figure out a way to automate it.  For right now, I'll assume you know the location.  You need to set CPPFLAGS to include that location.  For example:

~~~{sh}
python run_change_selection_coefficient.py
~~~

That command will:

1. compile our custom module if needed.
2. import it
3. do some evolution, apply the new custom function, and print some output

Please note that the "pyxbld" files attempt to deduce the install location of the fwdpy headers.  If that does not work, please post details to the fwdpy-users Google Group.