package bdmmprime.mapping;

import bdmmprime.parameterization.TypeSet;
import beast.core.BEASTObject;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;

/**
 * Logger for generating statistics from type mapped trees.
 */
public class TypedTreeStatsLogger extends CalculationNode implements Loggable {

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
            "If true, include changes occurring on root edge in count.",
            false);

    int[][] countMatrix;
    double[] lengthVector;

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

        countMatrix = new int[nTypes][nTypes];
        lengthVector = new double[nTypes];
    }


    private void update() {
        // Zero length and count arrays
        for (int i=0; i<nTypes; i++) {
            lengthVector[i] = 0.0;
            for (int j=0; j<nTypes; j++) {
                countMatrix[i][j] = 0;
            }
        }

        Node rootNode = tree.getRoot();
        if (!includeRootEdge) {
            while (rootNode.getChildCount() == 1)
                rootNode = rootNode.getChild(0);
        }

        computeMetricsOnSubtree(rootNode);
    }

    private int getType(Node node) {
        Object typeObj = node.getMetaData(typeLabel);
        if (typeObj == null)
            throw new RuntimeException("Tree does not have type metadata with the label '" + typeLabel + "'");

        if (typeObj instanceof Integer)
            return (int)typeObj;

        if (typeObj instanceof String)
            return typeSet.getTypeIndex((String)typeObj);

        throw new RuntimeException("Tree does not contain valid type metadata.");
    }

    private void computeMetricsOnSubtree(Node subtreeRoot) {

        Node node = subtreeRoot;
        int type = getType(subtreeRoot);

        while (node.getChildCount()==1) {
            Node childNode = node.getChild(0);
            int childType = getType(childNode);
            countMatrix[type][childType] += 1;
            lengthVector[childType] += node.getHeight() - childNode.getHeight();

            type = childType;
            node = childNode;
        }

        for (Node childNode : node.getChildren()) {
            int childType = getType(childNode);

            if (childType != type)
                countMatrix[type][childType] += 1;

            lengthVector[childType] += node.getHeight() - childNode.getHeight();

            computeMetricsOnSubtree(childNode);
        }
    }


    @Override
    public void init(PrintStream out) {

        String prefix = tree.getID() != null
                ? tree.getID() + "."
                : "";

        for (int type=0; type<nTypes; type++) {
            out.print(prefix + "length_" + typeSet.getTypeName(type) + "\t");

            for (int typeP=0; typeP<nTypes; typeP++) {
                if (type == typeP)
                    continue;

                out.print(prefix + "count_" + typeSet.getTypeName(type)
                        + "_to_" + typeSet.getTypeName(typeP) + "\t");

            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (tree instanceof TypeMappedTree)
            ((TypeMappedTree)tree).remapForLog(sample);

        update();

        for (int type=0; type<nTypes; type++) {
            out.print(lengthVector[type] + "\t");

            for (int typeP = 0; typeP < nTypes; typeP++) {
                if (type == typeP)
                    continue;

                out.print(countMatrix[type][typeP] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}
