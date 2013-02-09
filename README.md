mcohconstraints
===============

Matlab/Octave functions for generating coherence and avoiding sure loss constraints for lower previsions


License:
    https://creativecommons.org/licenses/by-sa/3.0/


The functions are documented in the m-files, and this documentation can be
accessed from the Matlab/Octave console with ```help <function_name>```.

They have been tested with Octave 3.4.3 and Matlab R2012a
(i.e., 7.14.0.739) on 64bit Linux.

Some functions call some non-standard supporting functions, i.e., you need
cddmex and/or bensolve:

* cddmex: 
    http://control.ee.ethz.ch/~hybrid/cdd.php
(I last used version version 1.00)

* bensolve:
    http://ito.mathematik.uni-halle.de/~loehne/index_en_dl.php
(I last used version 1.2)


Notes about cddmex and bensolve:

* cddmex comes with precompiled mexfiles linked to cddlib; when using it on a
recent linux distribution, my advice is to install cddlib using your package
manager and then compile the included c-file into a mexfile yourself.
For me, it sufficed to type the following on the Matlab/Octave command prompt:
```mex -v cddmex.c -lcdd```

* bensolve's code contains Matlabisms that Octave chokes on; to get it running
in Octave sufficiently for mcohconstraints, just replace all occurences of
`display(` by `disp(` in bensolve.m.
