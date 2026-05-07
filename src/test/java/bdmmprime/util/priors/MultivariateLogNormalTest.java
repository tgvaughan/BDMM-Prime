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

package bdmmprime.util.priors;

import bdmmprime.testclasses.RealVectorParamFromString;
import beast.base.inference.Distribution;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;


public class MultivariateLogNormalTest {

    @Test
    public void testValue() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealVectorParamFromString<>("1", PositiveReal.INSTANCE),
                "x", new RealVectorParamFromString<>("2", PositiveReal.INSTANCE),
                "M", new RealVectorParamFromString<>("0 0", Real.INSTANCE),
                "S", new RealVectorParamFromString<>("1 0.1 0.1 1", PositiveReal.INSTANCE));

        assertEquals(-2.8181831410741953, dist.calculateLogP(), 1e-10);
    }

    @Test
    public void testNonSymmetric() {
        assertThrows(IllegalArgumentException.class, () -> {
            Distribution dist = new MultivariateLogNormal();
            dist.initByName("x", new RealVectorParamFromString<>("1", PositiveReal.INSTANCE),
                    "x", new RealVectorParamFromString<>("2", PositiveReal.INSTANCE),
                    "M", new RealVectorParamFromString<>("0 0", Real.INSTANCE),
                    "S", new RealVectorParamFromString<>("1 1.1 1.2 1", PositiveReal.INSTANCE));

            System.out.println(dist.calculateLogP());
        });
    }

    @Test
    public void testNonPositiveDefinite() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealVectorParamFromString<>("1", PositiveReal.INSTANCE),
                "x", new RealVectorParamFromString<>("2", PositiveReal.INSTANCE),
                "M", new RealVectorParamFromString<>("0 0", Real.INSTANCE),
                "S", new RealVectorParamFromString<>("1 1.1 1.1 1", PositiveReal.INSTANCE));

        assertEquals(Double.NEGATIVE_INFINITY, dist.calculateLogP(), 0.0);
    }
}
