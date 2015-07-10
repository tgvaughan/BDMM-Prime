package beast.evolution.tree;

import beast.core.*;

import java.io.PrintStream;

/**
 * User: Denise
 * Date: 12.07.14
 * Time: 18:45
 */


@Description("Logger to report type changes along root branch")
public class OriginTypeChangeLogger extends CalculationNode implements Loggable, Function {

    public Input<MultiTypeRootBranch> multiTypeRootBranchInput =
            new Input<MultiTypeRootBranch>("multiTypeRootBranch", "MultiTypeRootBranch for origin coloring", Input.Validate.REQUIRED);



    private MultiTypeRootBranch multiTypeRootBranch;

    @Override
    public void initAndValidate() {

        multiTypeRootBranch = multiTypeRootBranchInput.get();
    }

    @Override
    public void init(PrintStream out) throws Exception {
            out.print(multiTypeRootBranch.getID() + "\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(multiTypeRootBranch.toString());
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        return multiTypeRootBranch.getChangeCount();
    }

    @Override
    public double getArrayValue(int iDim) {
        return multiTypeRootBranch.getChangeCount();
    }
}

