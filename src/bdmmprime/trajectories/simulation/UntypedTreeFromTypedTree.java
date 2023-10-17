package bdmmprime.trajectories.simulation;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

public class UntypedTreeFromTypedTree extends Tree {

    public Input<Tree> typedTreeInput = new Input<>("typedTree",
            "Typed tree used to initialize untyped tree.",
            Input.Validate.REQUIRED);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used in typed tree.",
            "type");

    Tree typedTree;
    int nextIntNr;

    String typeLabel;

    @Override
    public void initAndValidate() {

        typedTree = typedTreeInput.get();
        typeLabel = typeLabelInput.get();

        nextIntNr = typedTree.getLeafNodeCount();

        Node newRoot = getUntypedTree(typedTree.getRoot());

        assignFromWithoutID(new Tree(newRoot));

        super.initAndValidate();
    }

    Node getUntypedTree(Node cladeRoot) {

        while (cladeRoot.getChildCount() == 1)
            cladeRoot = cladeRoot.getChild(0);

        Node newCladeRoot = new Node();
        newCladeRoot.setHeight(cladeRoot.getHeight());

        if (cladeRoot.isLeaf()) {
            newCladeRoot.setNr(cladeRoot.getNr());
            newCladeRoot.setID(cladeRoot.getID());
            newCladeRoot.setMetaData(typeLabel,
                    cladeRoot.getMetaData(typeLabel));
            newCladeRoot.metaDataString = cladeRoot.metaDataString;
        } else {
            for (Node child : cladeRoot.getChildren())
                newCladeRoot.addChild(getUntypedTree(child));

            newCladeRoot.setNr(nextIntNr++);
        }

        return newCladeRoot;
    }
}


