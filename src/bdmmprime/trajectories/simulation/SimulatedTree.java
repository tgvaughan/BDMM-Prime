package bdmmprime.trajectories.simulation;

import bdmmprime.parameterization.*;
import bdmmprime.trajectories.*;
import bdmmprime.trajectories.trajevents.*;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.util.Randomizer;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static bdmmprime.util.Utils.nextBinomial;
import static org.apache.commons.math3.stat.StatUtils.sum;

/**
 * Simulates a tree from a multi-type birth-death skyline process.
 *
 * Note that the time of origin is also sampled, as this is a random
 * variable that depends on the time of the most recent sample in the tree.
 */
public class SimulatedTree extends Tree {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The equilibrium frequencies for each type",
            Input.Validate.REQUIRED);


    public Input<Integer> minSamplesInput = new Input<>("minSamples",
            "Minimum number of samples to accept in simulated trajectory.", 1);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Used to label tree nodes with corresponding types.",
            "type");

    public Input<String> trajFileNameInput = new Input<>("trajFileName",
            "Name of file to write simulated trajectory to.");

    public Input<String> treeFileNameInput = new Input<>("treeFileName",
            "Name of file to write simulated tree to.");

    public Input<RealParameter> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "Will be set to the time between the final sample and the end of the simulation.",
            Input.Validate.REQUIRED);

    Parameterization param;
    RealParameter frequencies;
    double simulationTime;

    double[] a_birth, a_death, a_sampling;
    double[][] a_migration, a_crossbirth;

    int nTypes;
    String typeLabel;

    int minSamples;

    public Trajectory traj;

    @Override
    public void initAndValidate() {
        param = parameterizationInput.get();
        frequencies = frequenciesInput.get();
        simulationTime = param.originInput.get().getValue();

        minSamples = minSamplesInput.get();

        nTypes = param.getNTypes();
        typeLabel = typeLabelInput.get();

        a_birth = new double[nTypes];
        a_death = new double[nTypes];
        a_sampling = new double[nTypes];
        a_migration = new double[nTypes][nTypes];
        a_crossbirth = new double[nTypes][nTypes];

        traj = null;
        do {
            traj = simulateTrajectory();
        } while (traj.getSampleCount() < Math.max(minSamples,1));

        finalSampleOffsetInput.get().setValue(param.originInput.get().getValue() - traj.getFinalSampleTime());

        if (trajFileNameInput.get() != null) {
            try (PrintStream out = new PrintStream(trajFileNameInput.get())) {

                traj.dump(out, false);

            } catch (FileNotFoundException e) {
                System.err.println("Error writing trajectory to file.");
            }
        }

        Tree tree = simulateTree();
        assignFromWithoutID(tree);

        if (treeFileNameInput.get() != null) {
            try (PrintStream out = new PrintStream(treeFileNameInput.get())) {

                out.println(tree);

            } catch (FileNotFoundException e) {
                System.err.println("Error writing tree to file.");
            }
        }

        super.initAndValidate();
    }

    Trajectory simulateTrajectory() {
        double t = 0;
        int interval = 0;

        double[] initialState = new double[nTypes];
        int startType;
        double u = Randomizer.nextDouble()*sum(frequencies.getDoubleValues());
        for (startType=0; startType<nTypes-1; startType++) {
            if (u < frequencies.getValue(startType))
                break;
            u -= frequencies.getValue(startType);
        }
        initialState[startType] = 1.0;

        Trajectory traj = new Trajectory(initialState);
        while (true) {

            double a_tot = 0.0;
            for (int s=0; s<nTypes; s++) {
                a_birth[s] = traj.currentState[s] * param.getBirthRates()[interval][s];
                a_death[s] = traj.currentState[s] * param.getDeathRates()[interval][s];
                a_sampling[s] = traj.currentState[s] * param.getSamplingRates()[interval][s];
                a_tot += a_birth[s] + a_death[s] + a_sampling[s];

                for (int sp = 0; sp < nTypes; sp++) {
                    if (sp == s)
                        continue;

                    a_migration[s][sp] = traj.currentState[s] * param.getMigRates()[interval][s][sp];
                    a_crossbirth[s][sp] = traj.currentState[s] * param.getCrossBirthRates()[interval][s][sp];
                    a_tot += a_migration[s][sp] + a_crossbirth[s][sp];
                }
            }

            if (a_tot > 0)
                t += Randomizer.nextExponential(a_tot);
            else
                t += Double.POSITIVE_INFINITY;

            if (param.getIntervalEndTimes()[interval] < simulationTime
                    && t > param.getIntervalEndTimes()[interval]) {
                t = param.getIntervalEndTimes()[interval];

                for (int s=0; s<nTypes; s++) {
                    double rho = param.getRhoValues()[interval][s];
                    if (rho > 0) {
                        int nRhoSamp = nextBinomial((int)Math.round(traj.currentState[s]), rho);
                        int nRemoveSamp = nextBinomial(nRhoSamp, param.getRemovalProbs()[interval][s]);
                        int nNoRemoveSamp = nRhoSamp - nRemoveSamp;

                        if (nRhoSamp > 0)
                            traj.addEvent(new SamplingEvent(t, s, nRemoveSamp, nNoRemoveSamp));
                    }
                }

                interval += 1;
                continue;
            }

            if (t > simulationTime)
                break;

            u = Randomizer.nextDouble()*a_tot;

            TrajectoryEvent event = null;
            for (int s=0; s<nTypes; s++) {

                if (u < a_birth[s]) {
                    event = new BirthEvent(t, s);
                    break;
                }
                u -= a_birth[s];

                if (u < a_death[s]) {
                    event = new DeathEvent(t, s);
                    break;
                }
                u -= a_death[s];

                if (u < a_sampling[s]) {
                    if (u < param.getRemovalProbs()[interval][s]*a_sampling[s])
                        event = new SamplingEvent(t, s, 1, 0);
                    else
                        event = new SamplingEvent(t, s, 0, 1);
                    break;
                }
                u -= a_sampling[s];

                for (int sp=0; sp<nTypes; sp++) {
                    if (sp == s)
                        continue;

                    if (u < a_migration[s][sp]) {
                        event = new MigrationEvent(t, s, sp);
                        break;
                    }
                    u -= a_migration[s][sp];

                    if (u < a_crossbirth[s][sp]) {
                        event = new CrossBirthEvent(t, s, sp);
                        break;
                    }
                    u -= a_crossbirth[s][sp];
                }

                if (event != null)
                    break;
            }

            if (event == null)
                throw new IllegalStateException("Event selection loop fell through.");

            traj.addEvent(event);
        }

        return traj;
    }

    public Tree simulateTree() {


        List<TrajectoryEvent> events = new ArrayList<>(traj.events);
        Collections.reverse(events);

        double[] state = traj.currentState.clone();

        List<List<Node>> activeLineages = new ArrayList<>();
        for (int s=0; s<nTypes; s++)
            activeLineages.add(new ArrayList<>());

        NodeFactory nodeFactory = new NodeFactory(traj.getFinalSampleTime(), traj.getSampleCount(),
                typeLabel, param.getTypeSet());

        for (TrajectoryEvent event : events) {
                event.simulateTreeEvent(state, activeLineages, nodeFactory);
                event.reverseUpdateState(state);
        }

        int nRemainingLineages = 0;
        for (int s=0; s<nTypes; s++)
            nRemainingLineages += activeLineages.get(s).size();

        if (nRemainingLineages != 1)
            throw new IllegalStateException("Number of remaining lineages not equal to 1.");

        Node root = null;
        for (int s=0; s<nTypes; s++) {
            if (!activeLineages.get(s).isEmpty()) {
                root = activeLineages.get(s).get(0);
                break;
            }
        }

        if (root == null)
            throw new IllegalStateException("Tree simulation failed.");

        return new Tree(root);
    }



    @Override
    public void log(long sample, PrintStream out) {
        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String newick = tree.getRoot().toSortedNewick(new int[1], true);
        out.print(newick);
        out.print(";");
    }

    /**
     * Debug method for testing
     * @param args unused
     */
    public static void main(String[] args) {

        Randomizer.setSeed(1);

        Parameterization param = new CanonicalParameterization();

        param.initByName(
                "typeSet", new TypeSet(2),
                "origin", "3.0",
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("4.0"),
                        new RealParameter("2.0 0.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "rhoSampling", new TimedParameter(
                        new RealParameter("2.0 3.0 4.0"),
                        new RealParameter("0.5 0.5 0.5"), 2));

        SimulatedTree sim = new SimulatedTree();
        sim.initByName("parameterization", param,
                "frequencies", new RealParameter("0.5 0.5"),
                "simulationTime", "5.0");

        System.out.println(sim.getRoot().toSortedNewick(new int[1], true));
    }
}
