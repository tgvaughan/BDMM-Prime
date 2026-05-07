/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.util;

import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;

import java.util.stream.Collectors;

/**
 * TraitSet in which value input is optional and the values
 * are initialized to a stand-in value.  Used by the BDMM-Prime BEAUti template,
 * where the trait set must be specified before any value (or indeed
 * the taxa themselves) can be known.
 *
 * Use of this trait set class in the BEAUti template is also important as
 * it causes BDMM-Prime's own traitset input editor to be used rather than
 * some other input editor (such as the one that MTT provdes).
 */
public class InitializedTraitSet extends TraitSet {

    public InitializedTraitSet() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        if (traitsInput.get() == null) {
            String value = taxaInput.get().getTaxaNames().stream()
                    .map(n -> n + "=NOT_SET")
                    .collect(Collectors.joining(","));

            traitsInput.setValue(value, this);
        }

        super.initAndValidate();
    }
}
