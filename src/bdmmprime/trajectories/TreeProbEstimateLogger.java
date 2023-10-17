package bdmmprime.trajectories;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

public class TreeProbEstimateLogger extends CalculationNode implements Loggable {

    public Input<SampledTrajectory> sampledTrajectoryInput = new Input<>("sampledTrajectory",
            "Sampled trajectory object.", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() { }

    @Override
    public void init(PrintStream out) {
        if (getID() != null)
            out.print(getID() + "\t");
        else
            out.print("logProbEst\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(sampledTrajectoryInput.get().getLogTreeProbEstimate() + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
