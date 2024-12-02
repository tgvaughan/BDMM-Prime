BDMM-Prime
==========
[![Unit/integration tests](https://github.com/tgvaughan/BDMM-Prime/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/tgvaughan/BDMM-Prime/actions/workflows/main.yml)

The BDMM-Prime project provides a [BEAST 2](http://www.beast2.org/) package for
performing phylodynamic inference under both structured and unstructured
birth-death models.

The BDMM-Prime project is a fork of the original
[BDMM project](https://github.com/denisekuehnert/bdmm).  The intention is to
extend the functionality of the original package, while improving its
flexibility and ease of use.  It incorporates the following enhancements:
- an improved BEAUti interface that allows a much more diverse range of analyses to be configured,
- automatic fall-back to analytical solutions for unstructured (single
  type) analyses (meaning BDMM-Prime includes much of the
  functionality of [BDSKY](https://github.com/BEAST2-Dev/bdsky)),
- use of stochastic mapping for sampling ancestral states,
- a particle filtering algorithm allowing joint sampling of population trajectories,
- a heavily refactored code base intended to make the package easier to use,
  extend and maintain.
  
As a result of the many changes that were required in making this transition,
BDMM-Prime is completely incompatible with BDMM itself.  Thus the original
package will be maintained separately to ensure that BEAST 2 XMLs and packages
that depend on it remain usable.

This repository is primarily of interest to people keen on building BDMM-Prime
from the source or contributing to its development.  If this doesn't include
you, please instead visit the project web page at:

https://tgvaughan.github.io/BDMM-Prime

There you'll find all relevant usage information including installation,
a complete tutorial, as well as reference guides for using BDMM-Prime via
the BEAUti interface together with detailed instructions on its XML interface.


Building from Source
--------------------

To build BDMM-Prime from source you'll need the following to be installed:
- OpenJDK version 17 or greater
- A recent version of OpenJFX
- the Apache Ant build system

Once these are installed and in your execution path, issue the following
command from the root directory of this repository:

```sh
JAVA_FX_HOME=/path/to/openjfx/ ant
```
The package archive will be left in the `dist/` subdirectory.

Note that unless you already have a local copy of the latest
[BEAST 2 source](https://github.com/CompEvol/beast2)
in the directory `../beast2` and the latest
[BeastFX source](https://github.com/CompEvol/beastfx)
in the directory `../beastfx` relative to the BDMM-Prime root, the build
script will attempt to download them automatically.  Thus, most builds
will require a network connection.


Acknowledgements and Citations
------------------------------

As this is a fork of [BDMM](https://github.com/denisekuehnert/bdmm),
BDMM-Prime owes its existence to the authors and contributors of
that project, in particular [Denise Kühnert](https://github.com/denisekuehnert/)
and [Jérémie Scire](https://github.com/jscire).

If you use this package as part of your research, please cite these papers:

* Vaughan and Stadler, "Bayesian phylodynamic inference of multi-type population trajectories using genomic data", [doi:10.1101/2024.11.26.625381](https://doi.org/10.1101/2024.11.26.625381) (preprint)
* Scire et al., "Robust Phylodynamic Analysis of Genetic Sequencing
  Data from Structured Populations", Viruses, 14:8, 1648 (2022),
  [doi:10.3390/v14081648](https://doi.org/10.3390/v14081648).
* Kühnert, et al., "Phylodynamics with Migration: A
  ComputationalFramework to Quantify Population Structure from Genomic
  Data", MBE, 33:8, 2102-2116 (2016),
  [doi:10.1093/molbev/msw064](https://doi.org/10.1093/molbev/msw064).
  

License
-------

BDMM-Prime is free software.  It is distributed under the terms of version 3
of the GNU General Public License.  A copy of this license should
be found in the file COPYING located in the root directory of this repository.
If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.
