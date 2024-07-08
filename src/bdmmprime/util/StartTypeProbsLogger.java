/*
 * Copyright (C) 2019-2024 Tim Vaughan
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

import bdmmprime.parameterization.TypeSet;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;

public class StartTypeProbsLogger extends BEASTObject implements Loggable {

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set used for analysis.", Input.Validate.REQUIRED);

    public Input<RealParameter> startTypeProbsInput = new Input<>("startTypeProbs",
            "Probabilities of different start types.",
            Input.Validate.REQUIRED);

    TypeSet typeSet;
    RealParameter startTypeProbs;
    String prefix;

    @Override
    public void initAndValidate() {
        typeSet = typeSetInput.get();
        startTypeProbs = startTypeProbsInput.get();

        prefix = (startTypeProbs.getID() != null ? startTypeProbs.getID()
                : "startTypeProbs") + ".";
    }

    @Override
    public void init(PrintStream out) {
        if (typeSet.getNTypes()==1)
            return;

        for (int i = 0; i < typeSet.getNTypes(); i++) {
            out.print(prefix + typeSet.getTypeName(i) + "\t");
        }

    }

    @Override
    public void log(long sample, PrintStream out) {
        if (typeSet.getNTypes()==1)
            return;

        startTypeProbs.log(sample, out);
    }

    @Override
    public void close(PrintStream out) {

    }
}
