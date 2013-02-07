mcohconstraints
===============

Matlab/Octave scripts for generating coherence and avoiding sure loss constraints for lower previsions


License:
    https://creativecommons.org/licenses/by-sa/3.0/


These files have been tested with Octave 3.4.3 and Matlab R2012a
(i.e., 7.14.0.739) on 64bit Linux.

To use them, you need cddmex and/or bensolve:

* cddmex: (I last used version version 1.00)
    http://control.ee.ethz.ch/~hybrid/cdd.php

* bensolve: (I last used version 1.2)
    http://ito.mathematik.uni-halle.de/~loehne/index_en_dl.php


Notes about cddmex and bensolve:

* cddmex comes with precompiled mexfiles linked to cddlib; when using it on a
recent linux distribution, my advice is to install cddlib using your package
manager and then compile the included c-file into a mexfile yourself.
For me, it sufficed to type the following on the Matlab/Octave command prompt:

    mex -v cddmex.c -lcdd

* bensolve's code contains Matlabisms that Octave chokes on; to get it running
in Octave sufficiently for mcohconstraints, just replace all occurences of
"display(" by "disp(" in bensolve.m (without the quotes, of course).
