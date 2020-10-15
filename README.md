BDMM-Prime
==========

[![Build Status](https://travis-ci.org/tgvaughan/BDMM-Prime.svg?branch=master)](https://travis-ci.org/tgvaughan/BDMM-Prime)

The BDMM-Prime project provides a [BEAST 2](http://www.beast2.org/) package for
performing phylodynamic inference under both structured and unstructured
birth-death models.

The BDMM-Prime project is a fork of the original
[BDMM project](https://github.com/denisekuehnert/bdmm).  The intention is to
extend the functionality of the original package, while improving its
flexibility and ease of use.  It incorporates the following enhancements:
- an improved BEAUti interface that allows a much more diverse range of analyses to be configured,
- automatic fall-back to analytical solutions for unstructured (single type) analyses,
- use of stochastic mapping for sampling ancestral states,
- a particle filtering algorithm allowing joint sampling of population trajectories,
- a heavily refactored code base intended to make the package easier to use,
  extend and maintain.

As a result of the many changes that were required in making this transition,
BDMM-Prime is completely incompatible with BDMM itself.  Thus the original
package will be maintained separately to ensure that BEAST 2 XMLs and packages
that depend on it remain usable.

BDMM-Prime is currently still in development and is not yet ready for general use.
Instead, refer to the original [BDMM project](https://github.com/denisekuehnert/bdmm).

However, if you are interested in playing around with the new version, you can install
it by adding https://tgvaughan.github.io/BDMM-Prime/package.xml as a third party
BEAST package repository and installing the package that appears.

About the Name
--------------

"BDMM-Prime" is named according to the
[mathematical convention](https://en.wikipedia.org/wiki/Prime_(symbol)) of using,
for example, x' (read "x prime") to name an alternative variable which is
in some sense related to x, but is nevertheless independent.

Importantly, it is not our intention to suggest that BDMM-Prime is superior
to BDMM - just different.

Acknowledgements
----------------

As this is a fork of [BDMM](https://github.com/denisekuehnert/bdmm),
BDMM-Prime owes its existance to all of the authors and contributors of
that project, in particular [Denise Kühnert](https://github.com/denisekuehnert/)
and [Jérémie Scire](https://github.com/jscire).

If you use this package as part of your research, please cite the
original BDMM paper:

* Kühnert, et al., "Phylodynamics with Migration: A
  ComputationalFramework to Quantify Population Structure from Genomic
  Data", MBE, 33:2102 (2016),
  [doi:10.1098/rsif.2013.1106](http://dx.doi.org/10.1098/rsif.2013.1106).

License
-------

BDMM-Prime is free software.  It is distributed under the terms of version 3
of the GNU General Public License.  A copy of this license should
be found in the file COPYING located in the root directory of this repository.
If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.
