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

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.io.PrintStream;

@Description("This class abstracts the birth-death process length, which can" +
        " be provided by either a parameter representing the \"origin\" of the" +
        " process, or by a tree in which case the process length is assumed to" +
        " be determined by the tree root height and a final sample offset." +
        " Used by the BDMM-Prime BEAUti template.")
public class ProcessLength extends CalculationNode implements Loggable, RealScalar<NonNegativeReal> {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree whose age might define the process length.");

    public Input<RealScalar<? extends NonNegativeReal>> originInput = new Input<>("origin",
            "Parameter whose value might define the process length.",
            Input.Validate.XOR, treeInput);

    public Input<RealScalar<? extends NonNegativeReal>> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "Final sample offset.  Only used in combination with the tree.",
            new RealScalarParam<>(0, NonNegativeReal.INSTANCE));

    public Input<Boolean> isEstimatedInput = new Input<>("estimate",
            "Indicates to BEAUti whether a prior on this parameter is " +
                    "needed.", true);

    public ProcessLength() {}

    public ProcessLength(Tree tree) {
        this.treeInput.setValue(tree, this);
    }

    public ProcessLength(RealScalar<? extends NonNegativeReal> origin) {
        this.originInput.setValue(origin, this);
    }

    @Override
    public void initAndValidate() { }

    @Override
    public NonNegativeReal getDomain() {
        return NonNegativeReal.INSTANCE;
    }

    public boolean isRoot() {
        return treeInput.get() != null;
    }

    @Override
    public double get() {
        if (treeInput.get() != null)
            return treeInput.get().getRoot().getHeight() +
                    finalSampleOffsetInput.get().get();
        else
            return originInput.get().get();
    }

    @Override
    public void init(PrintStream out) {
        out.print(getID() + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(get() + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
