/*
 * Copyright (C) 2019-2024 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.trajectories;

import bdmmprime.trajectories.trajevents.SamplingEvent;
import bdmmprime.trajectories.trajevents.TrajectoryEvent;

import java.io.PrintStream;
import java.util.*;
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
     * This method constructs a discretized copy of the current trajectory,
     * in which trajectory events are binned to composite events occuring
     * at regularly spaced intervals.
     *
     * @param timeStep time between composite events
     * @return discretized trajectory
     */
    public Trajectory getDiscretized(double timeStep) {

        List<TrajectoryEvent> discretizedEvents = new ArrayList<>();

        Map<String, TrajectoryEvent> eventMap = new HashMap<>();

        double gridTime = 0;

        for (TrajectoryEvent thisEvent : events) {

            if (thisEvent.time - gridTime > timeStep) {
                discretizedEvents.addAll(eventMap.values());
                eventMap.clear();
                gridTime += timeStep;
            }

            String thisEventFingerprint = thisEvent.getEventFingerprint();
            if (eventMap.keySet().contains(thisEventFingerprint)) {
                if (thisEvent instanceof SamplingEvent thisSamplingEvent) {
                    SamplingEvent discretizedSamplingEvent = (SamplingEvent) eventMap.get(thisEventFingerprint);
                    discretizedSamplingEvent.nRemoveSamp += thisSamplingEvent.nRemoveSamp;
                    discretizedSamplingEvent.nNoRemoveSamp += thisSamplingEvent.nNoRemoveSamp;
                }
                eventMap.get(thisEventFingerprint).multiplicity += thisEvent.multiplicity;
            } else {
                TrajectoryEvent eventCopy = thisEvent.copy();
                eventCopy.time = gridTime;
                eventMap.put(thisEventFingerprint, eventCopy);
            }

        }

        // Flush any remaining events
        if (!eventMap.isEmpty())
            discretizedEvents.addAll(eventMap.values());

        return new Trajectory(currentState.clone(), discretizedEvents);
    }

    /**
     * Method used to construct trajectory log files which can be directly
     * read into R without the need for custom parsing code.
     *
     * @param ps
     * @param sample
     * @param states
     * @param stateIdx
     * @param isFirst
     */
    public static void addToLog(PrintStream ps, long sample,
                                List<double[]> states,
                                List<TrajectoryEvent> events,
                                int stateIdx, double processLength,
                                boolean isFirst) {

        double[] state = states.get(stateIdx);
        double eventTime = stateIdx > 0 ? events.get(stateIdx-1).time : 0.0;
        double eventAge = processLength - eventTime;

        // Only include population size following last event recorded at a given time.
        boolean includePopSize = stateIdx==states.size()-1 || events.get(stateIdx).time>eventTime;

        for (int s=0; s<state.length; s++) {
            if (s>0 || !isFirst ) {
                ps.print("\n" + sample + "\t");
            }

            ps.print(eventTime + "\t" + eventAge + "\t");

            if (includePopSize)
                ps.print("N\t" + s + "\tNA\t" + state[s]);
        }

        ps.print("\n" + sample + "\t" + eventTime + "\t" + eventAge + "\t");

        if (stateIdx==0) {
            ps.print("O\tNA\tNA\tNA");
        } else {
            ps.print(events.get(stateIdx-1).getEventCode());
        }
    }

    public static void init(PrintStream out) {
        out.print("t\tage\tvariable\ttype\ttype2\tvalue");
    }

    public static void logEmpty(PrintStream out) {
        out.print("NA\tNA\tNA\tNA\tNA\tNA");
    }

    public static void log(long sample, List<double[]> states,
                            List<TrajectoryEvent> events,
                            double processLength,
                            PrintStream out) {

        addToLog(out, sample, states, events, 0, processLength, true);

        for (int i = 1; i < states.size(); i++)
            addToLog(out, sample, states, events, i, processLength, false);

        out.print("\t");
    }
}
