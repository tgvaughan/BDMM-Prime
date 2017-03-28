package beast.evolution.speciation;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.TreeInterface;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;

import java.util.Arrays;

/**
 * Created by Denise Kühnert on 06.03.17.
 */
public class BirthDeathMigrationClusterModelUncoloured extends BirthDeathMigrationModelUncoloured {

    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches. Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.REQUIRED);

    final public Input<IntegerParameter> clusterNumbers = new Input<>("clusterNumbers", "the names of all clusters (need to be integers)");

    final public Input<Integer> currentCluster = new Input<>("currentCluster", "the number of the current cluster");

    final public Input<RealParameter> quantileInput = new Input<>("rateQuantiles", "the rate quantiles associated with clusters for sampling of individual rates among branches.", Input.Validate.REQUIRED);


    ParametricDistribution distribution;
    IntegerParameter categories;
    RealParameter quantiles;
    private int clusterCount;
    private Integer[] clusters;
    private Integer[] clusterIndices;

    private boolean recompute = true;

    @Override
    public void initAndValidate() {

        super.initAndValidate();

        clusters = clusterNumbers.get().getValues();
        clusterCount = clusters.length;

        Integer max = 0;
        for (int i = 0; i < clusters.length; ++i) {
            max = Math.max(max, clusters[i]);
        }

        clusterIndices = new Integer[max]; // todo: stop wasting so much memory!!
        Arrays.fill(clusterIndices,0);

        for (int i = 0; i < clusters.length; ++i) {
            clusterIndices[clusters[i]-1] = i;
        }

        quantiles = quantileInput.get();
        quantiles.setDimension(clusterCount);
        Double[] initialQuantiles = new Double[clusterCount];
        for (int i = 0; i < clusterCount; i++) {
            initialQuantiles[i] = Randomizer.nextDouble();
        }
        RealParameter other = new RealParameter(initialQuantiles);
        quantiles.assignFromWithoutID(other);
        quantiles.setLower(0.0);
        quantiles.setUpper(1.0);

        distribution = rateDistInput.get();

        TreeInterface tree = treeInput.get();

        updateRates(tree);

//        try {
//            double mean = rateDistInput.get().getMean();
//            if (Math.abs(mean - 1.0) > 1e-6) {
//                Log.warning.println("WARNING: mean of distribution for BirthDeathMigrationClusterModelUncoloured is not 1.0.");
//            }
//        } catch (RuntimeException e) {
//            // ignore
//        }
    }

    protected Double updateRates(TreeInterface tree) {

        birth = new Double[n*totalIntervals];
        death = new Double[n*totalIntervals];
        psi = new Double[n*totalIntervals];
        b_ij = new Double[totalIntervals*(n*(n-1))];
        M = new Double[totalIntervals*(n*(n-1))];
        if (SAModel) r =  new Double[n * totalIntervals];

        if (transform) {
            transformParameters();
        }
        else {

            Double[] birthAmongDemesRates = new Double[1];

            if (birthAmongDemes) birthAmongDemesRates = birthRateAmongDemes.get().getValues();

            updateBirthDeathPsiParams();

            if (birthAmongDemes) {

                updateAmongParameter(b_ij, birthAmongDemesRates, b_ij_Changes, b_ijChangeTimes);
            }
        }

        Double[] migRates = migrationMatrix.get().getValues();

        updateAmongParameter(M, migRates, migChanges, migChangeTimes);

        updateRho();

        for (int i = 0; i < totalIntervals; i++) birth[i]*=getRateForCluster(clusterIndices[currentCluster.get()-1]);

        freq = frequencies.get().getValues();

        setupIntegrators();

        return 0.;
    }



    public double getRateForCluster(int cluster) {
        if (recompute) {
            // this must be synchronized to avoid being called simultaneously by
            // two different likelihood threads
//            synchronized (this) {
                distribution = rateDistInput.get(); //prepare();
                recompute = false;
//            }
        }

        return getRawRateForQuantile(cluster) * birth[0]; // todo: make this work for birth dimension > 1?
    }


    private double getRawRateForQuantile(int cluster) {

        try {
            return distribution.inverseCumulativeProbability(quantiles.getValue(cluster));
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability!");
        }
    }


//    @Override
//    public void store() {
//        storedScaleFactor = scaleFactor;
//        super.store();
//    }
//
//    @Override
//    public void restore() {
//        scaleFactor = storedScaleFactor;
//        super.restore();
//    }

}
