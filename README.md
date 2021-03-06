mcohconstraints
===============

Matlab/Octave functions for generating coherence and avoiding sure loss constraints for lower previsions and for calculating downward corrections of incoherent lower previsions

Code
----
**Usage** The functions are documented in the m-files, and this documentation can be
accessed from the Matlab/Octave console with ```help <function_name>```.

They have been tested with Octave 3.4.3 and a little bit with Matlab R2012a
(i.e., 7.14.0.739) on 64bit Linux.

**Licensing** CC BY-SA 3.0, see https://creativecommons.org/licenses/by-sa/3.0/

External code
-------------
Some functions call some non-standard supporting functions, i.e., you need
cddmex and/or bensolve:

* cddmex: 
    http://control.ee.ethz.ch/~hybrid/cdd.php
(I last used version version 1.00)

* bensolve:
    http://ito.mathematik.uni-halle.de/~loehne/index_en_dl.php
(I last used version 1.2)

Also, some functions use Octave's GLPK interface for linear programming.
These are not compatible with Matlab, but should be easy to convert to
linprog if needed. Another option would be to use glpkmex
(http://glpkmex.sourceforge.net/).

Notes about cddmex and bensolve:

* cddmex comes with precompiled mexfiles linked to cddlib; when using it on a
recent linux distribution, my advice is to install cddlib using your package
manager and then compile the included c-file into a mexfile yourself.
For me, it sufficed to type the following on the Matlab/Octave command prompt:
```mex -v cddmex.c -lcdd```

* bensolve's code contains Matlabisms that Octave chokes on; to get it running
in Octave sufficiently for mcohconstraints, just replace all occurences of
`display(` by `disp(` in bensolve.m.

References & citing
-------------------
This code has been used for my research.
This research has been presented at a conference:

> Quaeghebeur, Erik. 2010.
> “Characterizing coherence, correcting incoherence”
> Accepted for ISIPTA '13.

Please cite this paper if you use any material in this repository.