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

package bdmmprime.mapping;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

@Description("Logger for recording the duration of stochastic mapping calculations.")
public class TypeMappingDurationLogger extends CalculationNode implements Loggable {

    public Input<TypeMappedTree> typeMappedTreeInput = new Input<>(
            "typeMappedTree",
            "TypeMappedTree object used to stochastically map types onto " +
                    "a tip-typed tree.",
            Input.Validate.REQUIRED);

    TypeMappedTree tree;

    @Override
    public void initAndValidate() {
        tree = typeMappedTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        String prefix = tree.getID() != null
                ? tree.getID() + "."
                : "";

        out.print(prefix + "treeMappingCalcDuration" + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(tree.getPrevCalculationTime() + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
