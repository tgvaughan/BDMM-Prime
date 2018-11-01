package bdmm.loggers;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.speciation.BirthDeathMigrationModelUncoloured;

import java.io.PrintStream;

public class RootProbLogger extends BEASTObject implements Loggable {

    public Input<BirthDeathMigrationModelUncoloured> bdmmucInput = new Input<>(
            "bdmmuc",
            "Instance of BirthDeathMigrationModelUncoloured which records the root state probabilities",
            Input.Validate.REQUIRED);

    BirthDeathMigrationModelUncoloured bdmmuc;

    @Override
    public void initAndValidate() {
        bdmmuc = bdmmucInput.get();
    }

    @Override
    public void init(PrintStream out) {
        String loggerID;
        if (getID() != null)
            loggerID = getID() + ".";
        else if (bdmmuc.getID() != null)
            loggerID = bdmmuc.getID() + ".";
        else loggerID = "";

        double[] rootTypeProbs = bdmmuc.getRootTypeProbs();

        for (int i=0; i<rootTypeProbs.length; i++)
            out.print(loggerID + "probForRootType" + i + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (double rootTypeProb : bdmmuc.getRootTypeProbs())
            out.print(rootTypeProb + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
