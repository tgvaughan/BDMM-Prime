package bdmmprime.trajectories;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;

import java.util.List;
import java.util.Random;

public class TreeProbEstimateLogger extends Distribution {

    public Input<SampledTrajectory> sampledTrajectoryInput = new Input<>("sampledTrajectory",
            "Sampled trajectory logger.", Input.Validate.REQUIRED);

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }

    @Override
    public double calculateLogP() {
        logP = sampledTrajectoryInput.get().getTreeProbEstimate();
        return logP;
    }
}
