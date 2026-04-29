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

package bdmmprime.util.operators;

import beast.base.inference.parameter.RealParameter;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class SmartScaleOperatorTest {

    @Test
    public void test() {

        RealParameter startParam = new RealParameter("1 1 2 2 2");
        RealParameter param = new RealParameter("1 1 2 2 2");
        param.setLower(0.0);
        param.setUpper(Double.POSITIVE_INFINITY);
        SmartScaleOperator operator = new SmartScaleOperator();
        operator.initByName("parameter", param,
                "weight", 1.0);

        assertEquals(2, operator.nClasses);

        double hr = operator.proposal();
        double f = param.getArrayValue(0) != startParam.getArrayValue(0)
                ? param.getArrayValue(0)/startParam.getArrayValue(0)
                : param.getArrayValue(2)/startParam.getArrayValue(2);
        assertEquals(-Math.log(f), hr, 1e-10);
    }

    @Test
    public void testScaleAll() {

        RealParameter startParam = new RealParameter("1 1 2 2 2");
        RealParameter param = new RealParameter("1 1 2 2 2");
        param.setLower(0.0);
        param.setUpper(Double.POSITIVE_INFINITY);
        SmartScaleOperator operator = new SmartScaleOperator();
        operator.initByName("parameter", param,
                "scaleAll", true,
                "weight", 1.0);

        assertEquals(2, operator.nClasses);

        double hr = operator.proposal();
        double f = param.getArrayValue(0)/startParam.getArrayValue(0);
        assertEquals(-(operator.nClasses-2)*Math.log(f), hr, 1e-10);
    }
}
