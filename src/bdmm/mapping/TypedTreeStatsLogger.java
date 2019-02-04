package bdmm.mapping;

import bdmm.parameterization.TypeSet;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;

/**
 * Logger for generating statistics from type mapped trees.
 */
public class TypedTreeStatsLogger extends BEASTObject implements Loggable {

    public Input<Tree> typedTreeInput = new Input<>("typedTree",
            "Tree with type changes mapped.",
            Input.Validate.REQUIRED);

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set specifying model types.",
            Input.Validate.REQUIRED);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used to store type information on tree.",
            Input.Validate.REQUIRED);

    public Input<Boolean> includeRootEdgeChangesInput = new Input<>(
            "includeRootEdgeChanges",
            "If true, include changes occuring on root edge in count.",
            false);

    int[][] countMatrix;

    Tree tree;
    TypeSet typeSet;
    int nTypes;
    String typeLabel;
    boolean includeRootEdgeChanges;

    @Override
    public void initAndValidate() {
        tree = typedTreeInput.get();
        typeSet = typeSetInput.get();
        nTypes = typeSet.getNTypes();
        typeLabel = typeLabelInput.get();
        includeRootEdgeChanges = includeRootEdgeChangesInput.get();

        countMatrix = new int[nTypes][nTypes];
    }


    private void update() {
        // Zero count matrix
        for (int i=0; i<nTypes; i++) {
            for (int j=0; j<nTypes; j++) {
                countMatrix[i][j] = 0;
            }
        }

        Node rootNode = tree.getRoot();
        if (!includeRootEdgeChanges) {
            while (rootNode.getChildCount() == 1)
                rootNode = rootNode.getChild(0);
        }

        countChangesInSubtree(rootNode);
    }

    private int getType(Node node) {
        return (int)node.getMetaData(typeLabel);
    }

    private void countChangesInSubtree(Node subtreeRoot) {

        Node node = subtreeRoot;
        int type = getType(subtreeRoot);

        while (node.getChildCount()==1) {
            Node childNode = node.getChild(0);
            int childType = getType(childNode);
            countMatrix[type][childType] += 1;

            type = childType;
            node = childNode;
        }

        for (Node childeNode : node.getChildren()) {
            int childType = getType(childeNode);

            if (childType != type)
                countMatrix[type][childType] += 1;

            countChangesInSubtree(childeNode);
        }
    }


    @Override
    public void init(PrintStream out) {

        String prefix = tree.getID() != null
                ? tree.getID() + "."
                : "";

        for (int type=0; type<nTypes; type++) {
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
        update();

        for (int type=0; type<nTypes; type++) {
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
