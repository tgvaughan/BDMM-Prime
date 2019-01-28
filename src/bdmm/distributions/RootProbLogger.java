package bdmm.distributions;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

public class RootProbLogger extends BEASTObject implements Loggable {

    public Input<BirthDeathMigrationDistribution> treePriorInput = new Input<>(
            "bdmmTreePrior",
            "Instance of BirthDeathMigrationModel which records the " +
                    "root type probabilities",
            Input.Validate.REQUIRED);

    BirthDeathMigrationDistribution treePrior;

    @Override
    public void initAndValidate() {
        treePrior = treePriorInput.get();
    }

    @Override
    public void init(PrintStream out) {
        String loggerID;
        if (getID() != null)
            loggerID = getID() + ".";
        else if (treePrior.getID() != null)
            loggerID = treePrior.getID() + ".";
        else loggerID = "";

        double[] rootTypeProbs = treePrior.getRootTypeProbs();

        for (int i=0; i<rootTypeProbs.length; i++)
            out.print(loggerID + "probForRootType" + i + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (double rootTypeProb : treePrior.getRootTypeProbs())
            out.print(rootTypeProb + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
