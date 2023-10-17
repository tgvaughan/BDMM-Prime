package bdmmprime.util;

import beast.base.core.BEASTObject;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class NodeAgeLogger extends BEASTObject implements Loggable {

    public Input<List<Tree>> treesInput = new Input<>("tree",
            "Tree whose node ages to log.", new ArrayList<>());

    public Input<List<Function>> fsosInput = new Input<>("finalSampleOffset",
            "Final sample offsets for each tree.", new ArrayList<>());

    List<Tree> trees;
    List<Function> fsos;

    int nNodes;

    @Override
    public void initAndValidate() {
        trees = treesInput.get();
        fsos = fsosInput.get();

        if (!fsos.isEmpty() && fsos.size() != trees.size())
            throw new IllegalArgumentException("If finalSampleOffset is provided," +
                    "there must be as many as there are trees.");

        nNodes = 0;
        for (Tree tree : trees)
            nNodes += tree.getNodeCount();
    }

    @Override
    public void init(PrintStream out) {
        int i=0;
        for (Tree tree : trees) {
            for (int j=0; j<tree.getNodeCount(); j++) {
                out.print((i++) + "\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int treeIdx = 0; treeIdx<trees.size(); treeIdx++) {
            Tree tree = trees.get(treeIdx);
            double fso = fsos.isEmpty() ? 0.0 : fsos.get(treeIdx).getArrayValue();
            for (Node node : tree.getNodesAsArray())
                out.print((node.getHeight() + fso) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
