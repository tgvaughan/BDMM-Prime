/*
 * Copyright (C) 2019-2026 ETH Zurich
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
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.parameterization.SkylineVectorParameter;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

@Description("Logs a typed node tree, annotating each node with the values of "
        + "one or more skyline vector parameters evaluated at the node time and type.")
public class SkylineNodeTreeLogger extends TypedNodeTreeLogger {

    public Input<Parameterization> parameterizationInput = new Input<>(
            "parameterization",
            "Parameterization used to convert node heights to times at which skyline parameters are evaluated.",
            Input.Validate.REQUIRED);

    public Input<List<SkylineVectorParameter>> skylineParametersInput = new Input<>(
            "skylineParameter",
            "One or more skyline vector parameters whose values are annotated at each node of the tree.",
            new ArrayList<>(), Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "Optional time offset between the final sample and the end of the birth–death process, "
                    + "used when converting node heights to absolute times. Default is 0.",
            new RealParameter("0.0"), Input.Validate.OPTIONAL);

    public Input<Integer> precisionInput = new Input<>(
            "precision",
            "Number of decimal places used when logging skyline parameter values. Default is 6.",
            6, Input.Validate.OPTIONAL);

    private List<SkylineVectorParameter> skylineParameters;
    private Parameterization parameterization;
    private double finalSampleOffset;
    private int precision;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        skylineParameters = skylineParametersInput.get();
        parameterization = parameterizationInput.get();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue();
        precision = precisionInput.get() ;
    }

    @Override
    public String getStrippedNewick(Node node) {
        // Iterate over all nodes including the root
        // and annotate each node with skyline parameter values
        // evaluated at the node time and type
        for (Node n : node.getAllChildNodesAndSelf()) {
            double nodeTime = parameterization.getNodeTime(n, finalSampleOffset);

            Object nodeType = n.getMetaData("type");
            int typeIndex = nodeType != null ? (Integer) nodeType : 0;

            for (int i = 0; i < skylineParameters.size(); i++) {
                SkylineVectorParameter param = skylineParameters.get(i);
                String paramName = param.getID() != null ? param.getID() : "param" + i;

                double value = param.getValuesAtTime(nodeTime)[typeIndex];
                BigDecimal rounded = new BigDecimal(value)
                        .setScale(precision, RoundingMode.HALF_UP);
                n.setMetaData(paramName, rounded);
            }

            n.metaDataString = n.getMetaDataNames().stream()
                    .map(k -> k + "=" + n.getMetaData(k))
                    .collect(Collectors.joining(","));
        }

        return super.getStrippedNewick(node);
    }
}
