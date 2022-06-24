package bdmmprime.mapping;

import beast.core.Input;
import beast.evolution.tree.Node;

import java.io.PrintStream;

public class TransitionTimeLogger extends AbstractTypeTreeStatsLogger {

    public Input<String> sourceTypeInput = new Input<>(
            "sourceType",
            "Source (oldest) type",
            Input.Validate.REQUIRED);

    public Input<String> destTypeInput = new Input<>(
            "destType",
            "Destination (younguest) type",
            Input.Validate.REQUIRED);

    private Double earliest, latest;
    private int count;

    int sourceTypeIdx, destTypeIdx;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        sourceTypeIdx = typeSet.getTypeIndex(sourceTypeInput.get());
        destTypeIdx = typeSet.getTypeIndex(destTypeInput.get());
    }

    public void update() {
        earliest = null;
        latest = null;
        count = 0;

        Node rootNode = tree.getRoot();
        if (!includeRootEdge) {
            while (rootNode.getChildCount() == 1)
                rootNode = rootNode.getChild(0);
        }
        updateMetricsOnSubtree(rootNode);
    }

    boolean isTransition(Node node) {
        if (node.isLeaf() || getType(node) != sourceTypeIdx)
            return false;

        for (Node child : node.getChildren())
            if (getType(child) == destTypeIdx)
                return true;

        return false;
    }

    public void updateMetricsOnSubtree(Node subtreeRoot) {
        if (isTransition(subtreeRoot)) {
            if (earliest == null)
                earliest = subtreeRoot.getHeight();

            if (latest == null || subtreeRoot.getHeight()<latest)
                latest = subtreeRoot.getHeight();

            count += 1;
        }

        for (Node child : subtreeRoot.getChildren())
            updateMetricsOnSubtree(child);
    }

    @Override
    public void init(PrintStream out) {
        String prefix = (tree.getID() != null ? tree.getID() + "." : "")
                + sourceTypeInput.get() + "_to_" + destTypeInput.get();

        out.print(prefix + "earliestTransitionAge" + "\t");
        out.print(prefix + "latestTransitionAge" + "\t");
        out.print(prefix + "transitionCount" + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (tree instanceof TypeMappedTree)
            ((TypeMappedTree)tree).remapForLog(sample);

        update();

        out.print( (earliest == null ? "NA" : earliest) + "\t"
                + (latest == null ? "NA" : latest) + "\t"
                + count + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
