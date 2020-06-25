package bdmmprime.simulation;

import bdmmprime.parameterization.Parameterization;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

public class SimulatedTree extends Tree {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The equilibrium frequencies for each type",
            Input.Validate.REQUIRED);

    public Input<RealParameter> maxSimulationTimeInput = new Input<>("maxSimulationTime",
            Input.Validate.REQUIRED);

    Parameterization param;
    RealParameter frequencies;
    RealParameter maxSimulationTime;


    @Override
    public void initAndValidate() {
        param = parameterizationInput.get();
        frequencies = frequenciesInput.get();
        maxSimulationTime = maxSimulationTimeInput.get();

        super.initAndValidate();
    }

    void simulateTrajectory() {
        double t = 0;

        while (true) {

        }
    }

    void simulate() {



    }
}
