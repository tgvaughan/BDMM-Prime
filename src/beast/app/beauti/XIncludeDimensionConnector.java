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
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.SCMigrationModel;
import beast.evolution.tree.StructuredCoalescentMultiTypeTree;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import multitypetree.distributions.ExcludablePrior;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;


/**
 * Class containing a static method used as a "custom connector" in the
 * MultiTypeBirthDeath BEAUti template.  This connector ensures that the
 * xInclude input for the samplingProportion priors exclude elements that
 * are zero when computing the prior.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class XIncludeDimensionConnector {

    public static boolean customConnector(BeautiDoc doc) {

        for (BEASTInterface p : doc.getPartitions("Tree")) {
            TreeLikelihood treeLikelihood = (TreeLikelihood) p;
            Tree tree =  (Tree) treeLikelihood.treeInput.get();

            String pID = BeautiDoc.parsePartition(tree.getID());

            ExcludablePrior excludablePrior = (ExcludablePrior) doc.pluginmap.get(
                    "samplingProportionPrior.t:" + pID);


            RealParameter samplingPropParam =  (RealParameter) excludablePrior.m_x.get();
            int samplingPropDim = samplingPropParam.getDimension();

            StringBuilder xIncludeStrBuilder = new StringBuilder();
            for (int i=0; i<samplingPropDim; i++) {
                if (samplingPropParam.getValue(i) > 0.0)
                    xIncludeStrBuilder.append("true ");
                else
                    xIncludeStrBuilder.append("false ");
            }

            excludablePrior.xIncludeInput.get().setDimension(samplingPropDim);
            excludablePrior.xIncludeInput.get().valuesInput.setValue(
                    xIncludeStrBuilder.toString(),
                    excludablePrior.xIncludeInput.get());

            try {
                excludablePrior.xIncludeInput.get().initAndValidate();
                excludablePrior.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error configuring initial migration model.");
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
