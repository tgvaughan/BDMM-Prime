package bdmmprime.mapping;

import bdmmprime.parameterization.TypeSet;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public abstract class AbstractTypeTreeStatsLogger extends CalculationNode implements Loggable {

    public Input<Tree> typedTreeInput = new Input<>("typedTree",
            "Tree with type changes mapped.",
            Input.Validate.REQUIRED);

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set specifying model types.",
            Input.Validate.REQUIRED);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used to store type information on tree.",
            Input.Validate.REQUIRED);

    public Input<Boolean> includeRootEdgeInput = new Input<>(
            "includeRootEdge",
            "If true, include root edge in summary stats calculations.",
            false);

    Tree tree;
    TypeSet typeSet;
    int nTypes;
    String typeLabel;
    boolean includeRootEdge;

    @Override
    public void initAndValidate() {
        tree = typedTreeInput.get();
        typeSet = typeSetInput.get();
        nTypes = typeSet.getNTypes();
        typeLabel = typeLabelInput.get();
        includeRootEdge = includeRootEdgeInput.get();

    }

    protected int getType(Node node) {
        Object typeObj = node.getMetaData(typeLabel);
        if (typeObj == null)
            throw new RuntimeException("Tree does not have type metadata with the label '" + typeLabel + "'");

        if (typeObj instanceof Integer)
            return (int)typeObj;

        if (typeObj instanceof String)
            return typeSet.getTypeIndex((String)typeObj);

        throw new RuntimeException("Tree does not contain valid type metadata.");
    }
}
