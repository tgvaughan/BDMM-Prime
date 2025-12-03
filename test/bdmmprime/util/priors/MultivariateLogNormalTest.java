/*
 * Copyright (C) 2019-2025 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.util.priors;

import beast.base.inference.Distribution;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class MultivariateLogNormalTest {

    @Test
    public void testValue() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealParameter("1"),
                "x", new RealParameter("2"),
                "M", new RealParameter("0 0"),
                "S", new RealParameter("1 0.1 0.1 1"));

        assertEquals(-2.8181831410741953, dist.calculateLogP(), 1e-10);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testNonSymmetric() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealParameter("1"),
                "x", new RealParameter("2"),
                "M", new RealParameter("0 0"),
                "S", new RealParameter("1 1.1 1.2 1"));

        System.out.println(dist.calculateLogP());
    }

    @Test
    public void testNonPositiveDefinite() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealParameter("1"),
                "x", new RealParameter("2"),
                "M", new RealParameter("0 0"),
                "S", new RealParameter("1 1.1 1.1 1"));

        assertEquals(Double.NEGATIVE_INFINITY, dist.calculateLogP(), 0.0);
    }
}
