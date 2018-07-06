/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package beast.app.beauti;

import beast.core.BEASTInterface;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import multitypetree.distributions.ExcludablePrior;

import java.util.List;
import java.util.stream.Collectors;


/**
 * Class containing a static method used as a "custom connector" in the
 * MultiTypeBirthDeath BEAUti template.  This connector ensures that the
 * trait set is initialized properly.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class TraitSetInitConnector {

    public static boolean customConnector(BeautiDoc doc) {

        for (BEASTInterface p : doc.getPartitions("Tree")) {
            TreeLikelihood treeLikelihood = (TreeLikelihood) p;
            Tree tree =  (Tree) treeLikelihood.treeInput.get();

            String pID = BeautiDoc.parsePartition(tree.getID());

            TraitSet typeTraitSet = (TraitSet) doc.pluginmap.get(
                    "typeTraitSet.t:" + pID);

            if (typeTraitSet.traitsInput.get() == null
                    || typeTraitSet.traitsInput.get().trim().isEmpty()) {

                try {
                    String newValue = typeTraitSet.taxaInput.get().getTaxaNames().stream()
                            .map(n -> n + "=NOT_SET").collect(Collectors.joining(","));
                    typeTraitSet.traitsInput.setValue(newValue, typeTraitSet);
                    typeTraitSet.initAndValidate();
                } catch (Exception ex) {
                    System.err.println("Error configuration initial migration model.");
                }
            }
        }

        return false;
    }

    private static String getParameterString(Parameter.Base param) {

        String str = "";
        for (Object value : (List<Object>) param.valuesInput.get()) {
            if (str.length()>0)
                str += " ";
            str += value;
        }

        return str;
    }
}
