package bdmm.loggers;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import bdmm.distributions.BirthDeathMigrationDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.io.PrintStream;

/**
 * Adapted from TypedNodeTreeLogger class (written by Tim Vaughan) in MultiTypeTree package
 */

public class TipTypedTreeLogger extends BEASTObject implements Loggable {
    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Sample-typed tree to log.",
            Input.Validate.REQUIRED);
    public Input<BirthDeathMigrationDistribution> bdmucInput = new Input<>(
            "bdmuc",
            "Birth-Death migration uncoloured model",
            Input.Validate.REQUIRED);
    public Input<String> typeLabelInput = new Input<>(
            "typeLabel",
            "string for type label",
            "type");

    Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        tree.init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {

        // Set up metadata string
        for (Node node : tree.getNodesAsArray()) {
            if(node.isLeaf() ) {
                node.metaDataString = typeLabelInput.get()
                        + "=\""
                        + bdmucInput.get().getLeafStateForLogging(node, nSample)
                        + "\"";
            }
        }

        out.print("tree STATE_" + nSample + " = ");
        out.print(tree.getRoot().toSortedNewick(new int[1], true));
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        tree.close(out);
    }
}
