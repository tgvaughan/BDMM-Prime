package bdmmprime.trajectories.simulation;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

public class SimulatedTrajectoryLogger extends BEASTObject implements Loggable {

    public Input<SimulatedTree> simulatedTreeInput = new Input<>("simulatedTree",
            "Simulated tree whose trajectory you want to log.",
            Input.Validate.REQUIRED);

    SimulatedTree simulatedTree;

    @Override
    public void initAndValidate() {
        simulatedTree = simulatedTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        if (getID() == null)
            out.print("trajectory\t");
        else
            out.print(getID() + "\t");
    }


    @Override
    public void log(long sample, PrintStream out) {
        if (simulatedTree.traj == null)
            out.print("NA");
        else
            out.print(simulatedTree.traj);

        out.print("\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
