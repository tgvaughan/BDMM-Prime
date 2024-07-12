package bdmmprime.distribution;

import bdmmprime.parameterization.TypeSet;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

public class StartTypePosteriorProbsLogger extends CalculationNode implements Loggable {

    public Input<BirthDeathMigrationDistribution> treePriorInput = new Input<>(
            "bdmmTreePrior",
            "Instance of BirthDeathMigrationModel which records the " +
                    "initial type posterior probabilities",
            Input.Validate.REQUIRED);

    BirthDeathMigrationDistribution treePrior;
    TypeSet typeSet;
    String prefix;

    @Override
    public void initAndValidate() {
        treePrior = treePriorInput.get();
        typeSet = treePrior.parameterizationInput.get().getTypeSet();

        String loggerID;
        if (getID() != null)
            loggerID = getID() + ".";
        else if (treePrior.getID() != null)
            loggerID = treePrior.getID() + ".";
        else loggerID = "";

        prefix = loggerID + "startTypePosteriorProbs.";
    }

    @Override
    public void init(PrintStream out) {
        if (typeSet.getNTypes()==1)
            return;

        for (int i=0; i<typeSet.getNTypes(); i++)
            out.print(prefix + typeSet.getTypeName(i) + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (typeSet.getNTypes()==1)
            return;

        for (double startTypeProb : treePrior.getStartTypePosteriorProbs())
            out.print(startTypeProb + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
