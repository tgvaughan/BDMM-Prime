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

import bdmmprime.testclasses.RealVectorParamFromString;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SmartScaleOperatorTest {

    @Test
    public void test() {

        RealVectorParam<?> startParam = new RealVectorParamFromString<>("1 1 2 2 2", NonNegativeReal.INSTANCE);
        RealVectorParam<?> param = new RealVectorParamFromString<>("1 1 2 2 2", NonNegativeReal.INSTANCE);
        SmartScaleOperator operator = new SmartScaleOperator();
        operator.initByName("parameter", param,
                "weight", 1.0);

        assertEquals(2, operator.nClasses);

        double hr = operator.proposal();
        double f = param.get(0) != startParam.get(0)
                ? param.get(0)/startParam.get(0)
                : param.get(2)/startParam.get(2);
        assertEquals(-Math.log(f), hr, 1e-10);
    }

    @Test
    public void testScaleAll() {

        RealVectorParam<?> startParam = new RealVectorParamFromString<>("1 1 2 2 2", NonNegativeReal.INSTANCE);
        RealVectorParam<?> param = new RealVectorParamFromString<>("1 1 2 2 2", NonNegativeReal.INSTANCE);
        SmartScaleOperator operator = new SmartScaleOperator();
        operator.initByName("parameter", param,
                "scaleAll", true,
                "weight", 1.0);

        assertEquals(2, operator.nClasses);

        double hr = operator.proposal();
        double f = param.get(0)/startParam.get(0);
        assertEquals(-(operator.nClasses-2)*Math.log(f), hr, 1e-10);
    }
}
