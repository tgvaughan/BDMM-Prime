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

package bdmmprime.testclasses;

import beast.base.spec.domain.Domain;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;

public class RealVectorParamFromString<T extends Real> extends RealVectorParam<T> {

    public RealVectorParamFromString(String initString, T domain) {
        String[] spl = initString.split(" ");
        double[] valArray = new double[spl.length];
        for (int i=0; i<valArray.length; i++)
            valArray[i] = Double.parseDouble(spl[i]);

        super(valArray, domain);
    }
}
