package bdmmprime.util.operators;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;

import java.util.*;

public abstract class SmartRealOperator extends Operator {

    public Input<List<RealVectorParam<? extends Real>>> parametersInput = new Input<>("parameter",
            "One or more parameters to operate on", new ArrayList<>());

    public Input<RealVector<? extends Real>> classesToExcludeInput = new Input<>("classesToExclude",
            "Elements with these value will not be operated on.");

    protected List<RealVectorParam<? extends Real>> parameters;
    protected Map<RealVectorParam<? extends Real>, Integer[]> groups;

    protected int nClasses;

    @Override
    public void initAndValidate() {

        parameters = parametersInput.get();

        SortedSet<Double> seenValuesSet = new TreeSet<>();

        for (RealVectorParam<? extends Real> param : parameters) {
            for (int i=0; i<param.size(); i++) {
                if (param.get(i) != 0.0)
                    seenValuesSet.add(param.get(i));
            }
        }

        // Explicitly exclude certain classes (identified by the element value)
        if (classesToExcludeInput.get() != null) {
            for (double value : classesToExcludeInput.get().getElements())
                seenValuesSet.remove(value);
        }

        List<Double> seenValues = new ArrayList<>(seenValuesSet);
        nClasses = seenValues.size();

        groups = new HashMap<>();

        for (RealVectorParam<? extends Real> param : parameters) {
            Integer[] groupIDs = new Integer[param.size()];

            for (int i = 0; i < param.size(); i++)
                groupIDs[i] = seenValues.indexOf(param.get(i));

            groups.put(param, groupIDs);
        }
    }

    @Override
    public List<StateNode> listStateNodes() {
        return new ArrayList<>(parametersInput.get());
    }
}
