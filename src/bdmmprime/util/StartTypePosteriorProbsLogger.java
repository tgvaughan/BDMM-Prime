/*
 * Copyright (C) 2019-2024 ETH Zurich
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

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.TypeSet;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

public class StartTypePosteriorProbsLogger extends CalculationNode implements Loggable {

    public Input<BirthDeathMigrationDistribution> treePriorInput = new Input<>(
            "bdmmTreePrior",
            "Instance of BirthDeathMigrationModel which records the " +
                    "initial type posterior probabilities",
            Input.Validate.REQUIRED);

    BirthDeathMigrationDistribution treePrior;
    TypeSet typeSet;
    String prefix;

    @Override
    public void initAndValidate() {
        treePrior = treePriorInput.get();
        typeSet = treePrior.parameterizationInput.get().getTypeSet();

        String loggerID;
        if (getID() != null)
            loggerID = getID() + ".";
        else if (treePrior.getID() != null)
            loggerID = treePrior.getID() + ".";
        else loggerID = "";

        prefix = loggerID + "startTypePosteriorProbs.";
    }

    @Override
    public void init(PrintStream out) {
        if (typeSet.getNTypes()==1)
            return;

        for (int i=0; i<typeSet.getNTypes(); i++)
            out.print(prefix + typeSet.getTypeName(i) + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (typeSet.getNTypes()==1)
            return;

        for (double startTypeProb : treePrior.getStartTypePosteriorProbs())
            out.print(startTypeProb + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
