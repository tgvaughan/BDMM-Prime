package bdmmprime.util.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;

import java.util.*;

public abstract class SmartRealOperator extends Operator {

    public Input<List<RealParameter>> parametersInput = new Input("parameter",
            "One or more parameters to operate on", new ArrayList<>());

    public Input<List<Double>> classesToExcludeInput = new Input<>("classToExclude",
            "Elements with this value will not be operated on.",
            new ArrayList<>());

    protected List<RealParameter> parameters;
    protected Map<RealParameter, Integer[]> groups;

    protected int nClasses;

    @Override
    public void initAndValidate() {

        parameters = parametersInput.get();

        SortedSet<Double> seenValuesSet = new TreeSet<>();

        for (RealParameter param : parameters) {
            for (int i=0; i<param.getDimension(); i++) {
                if (param.getValue(i) != 0.0)
                    seenValuesSet.add(param.getValue(i));
            }
        }

        // Explicitly exclude certain classes (identified by the element value)
        classesToExcludeInput.get().forEach(seenValuesSet::remove);

        List<Double> seenValues = new ArrayList<>(seenValuesSet);
        nClasses = seenValues.size();

        groups = new HashMap<>();

        for (RealParameter param : parameters) {
            Integer[] groupIDs = new Integer[param.getDimension()];

            for (int i = 0; i < param.getDimension(); i++)
                groupIDs[i] = seenValues.indexOf(param.getValue(i));

            groups.put(param, groupIDs);
        }
    }
}
