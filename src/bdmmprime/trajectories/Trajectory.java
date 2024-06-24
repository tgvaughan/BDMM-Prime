package bdmmprime.trajectories;

import bdmmprime.trajectories.trajevents.TrajectoryEvent;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A simulated birth-death trajectory.
 *
 * At bare minimum, objects of this class have a list of events and the state following the final event.
 * This is sufficient to (a) continue simulating, or (b) reconstruct all previous states.
 */
public class Trajectory {
    public List<TrajectoryEvent> events = new ArrayList<>();
    public double[] currentState;

    public Trajectory(double[] currentState) {
        this.currentState = currentState.clone();
    }

    public Trajectory(double[] currentState, List<TrajectoryEvent> events) {
        this.currentState = currentState.clone();
        this.events = new ArrayList<>(events);
    }

    public void addEvent(TrajectoryEvent event) {
        events.add(event);
        event.updateState(currentState);
    }

    /**
     * Replace current trajectory events and state with those from
     * another trajectory.
     *
     * @param other Trajectory whose events and state will replace the current values.
     */
    public void assignFrom(Trajectory other) {
        events.clear();
        events.addAll(other.events);

        System.arraycopy(other.currentState, 0, currentState, 0, currentState.length);
    }

    public List<double[]> getStateList() {
        List<double[]> states = new ArrayList<>();

        double[] state = currentState.clone();
        states.add(state.clone());
        for (int i=events.size()-1; i>=0; i--) {
            events.get(i).reverseUpdateState(state);
            states.add(state.clone());
        }

        Collections.reverse(states);

        return states;
    }

    public boolean currentStateValid() {
        for (double v : currentState) {
            if (v < 0)
                return false;
        }

        return true;
    }

    public boolean currentStateValid(int[] lineages) {
        for (int s=0; s<currentState.length; s++) {
            if (Math.round(currentState[s]) < lineages[s])
                return false;
        }

        return true;
    }

    public boolean currentStateEmpty() {
        for (int s=0; s<currentState.length; s++) {
            if (Math.round(currentState[s]) > 0)
                return false;
        }

        return true;
    }

    public List<Double> getEventTimes() {
        return this.events.stream().map(e -> e.time).collect(Collectors.toList());
    }

    public int getSampleCount() {
        return (int) this.events.stream()
                .filter(TrajectoryEvent::isSamplingEvent)
                .mapToInt(e -> e.multiplicity)
                .sum();
    }

    public double getFinalSampleTime() {
        return this.events.stream()
                .filter(TrajectoryEvent::isSamplingEvent)
                .mapToDouble(e -> e.time)
                .max().getAsDouble();
    }

    public double getFirstSampleTime() {
        return this.events.stream()
                .filter(TrajectoryEvent::isSamplingEvent)
                .mapToDouble(e -> e.time)
                .min().getAsDouble();
    }

    public void dump(PrintStream out, boolean includeEvents) {
        List<double[]> states = getStateList();
        List<Double> eventTimes = getEventTimes();

        out.print("t");
        for (int s=0; s<currentState.length; s++)
            out.print("\tn" + s);

        if (includeEvents)
            out.print("\tevent");

        out.println();

        for (int i=0; i<states.size(); i++ ) {
            if (i == 0)
                out.print(0.0);
            else
                out.print(eventTimes.get(i - 1));

            for (int s = 0; s < currentState.length; s++) {
                out.print("\t" + states.get(i)[s]);
            }

            if (includeEvents) {
                if (i == 0)
                    out.print("\tSTART");
                else
                    out.print("\t" + events.get(i - 1));
            }

            out.println();
        }
    }


    /**
     * Method used to construct trajectory log files which can be directly
     * read into R without the need for custom parsing code.
     *
     * @param ps
     * @param sample
     * @param time
     * @param isFirst
     */
    public void addToLog(PrintStream ps, long sample, double time,
                         double[] state, TrajectoryEvent event,
                         boolean isFirst) {


        for (int s=0; s<state.length; s++) {
            if (s>0 || !isFirst ) {
                ps.print("\n" + sample + "\t");
            }

            ps.print(time + "\t");
            ps.print("N\t" + s + "\tNA\t" + state[s]);
        }

        ps.print("\n" + sample + "\t" + time + "\t");

        if (event == null) {
            ps.print("O\tNA\tNA\tNA");
        } else {
            ps.print(event.getEventCode());
        }
    }

    public void log (PrintStream out, long sample) {
        List<double[]> states = getStateList();
        List<Double> eventTimes = getEventTimes();

        addToLog(out, sample, 0, states.get(0), null, true);

        for (int i = 1; i < states.size(); i++) {
            addToLog(out, sample, eventTimes.get(i - 1),
                    states.get(i), events.get(i - 1), false);

        }
        out.print("\t");
    }
}
