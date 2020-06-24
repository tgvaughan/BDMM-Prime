package bdmmprime.trajectories;

import bdmmprime.parameterization.Parameterization;
import beast.util.Transform;

import java.util.ArrayList;
import java.util.List;

public class Trajectory {
    List<TrajectoryEvent> events = new ArrayList<>();
    double[] currentState;

    public Trajectory(double[] currentState) {
        this.currentState = currentState.clone();
    }

    public Trajectory(double[] currentState, List<TrajectoryEvent> events) {
        this.currentState = currentState.clone();
        this.events = new ArrayList<>(events);
    }

    public Trajectory copy() {
        return new Trajectory(currentState, events);
    }

    public void addEvent(TrajectoryEvent event) {
        events.add(event);
        event.updateState(currentState);
    }
}
