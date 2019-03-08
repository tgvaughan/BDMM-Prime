BDMM-Prime
==========

The BDMM-Prime project provides a [BEAST 2](http://www.beast2.org/) package for
performing phylodynamic inference under a structured birth-death model.

The BDMM-Prime project is a fork of the original [BDMM project](https://github.com/denisekuehnert/bdmm).
The intention is to maintain exactly the same the functionaility of the original package,
while improving its flexibility and ease of use. To this end, the fork
employs an exact stochastic mapping algorithm rather than MCMC in
order to sample state changes along tree lineages.  As a result of the
many changes that were required in making this transition, it is
completely incompatible with BDMM itself.

BDMM-Pime is currently still in development and is not yet ready for general use.
Instead, refer to the original [BDMM project](https://github.com/denisekuehnert/bdmm).

BDMM-Prime is free software.  It is distributed under the terms of version 3
of the GNU General Public License.  A copy of this license should
be found in the file COPYING located in the root directory of this repository.
If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.
