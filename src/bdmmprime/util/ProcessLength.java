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

package bdmmprime.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import feast.function.LoggableFunction;

@Description("This class abstracts the birth-death process length, which can" +
        " be provided by either a parameter representing the \"origin\" of the" +
        " process, or by a tree in which case the process length is assumed to" +
        " be determined by the tree root height and a final sample offset." +
        " Used primarily by the BDMM-Prime BEAUti template.")
public class ProcessLength extends LoggableFunction {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree whose age might define the process length.");

    public Input<Function> originInput = new Input<>("origin",
            "Parameter whose value might define the process length.",
            Input.Validate.XOR, treeInput);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "Final sample offset.  Only used in combination with the tree.",
            new RealParameter("0"));

    @Override
    public void initAndValidate() { }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue(int dim) {
        if (treeInput.get() != null)
            return treeInput.get().getRoot().getHeight() +
                    finalSampleOffsetInput.get().getArrayValue();
        else
            return originInput.get().getArrayValue();
    }
}
