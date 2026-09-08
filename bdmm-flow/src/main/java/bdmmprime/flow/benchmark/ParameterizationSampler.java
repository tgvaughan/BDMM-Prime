package bdmmflow.benchmark;

import bdmmprime.parameterization.*;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;
import java.util.Random;
import java.util.StringJoiner;

public class ParameterizationSampler {

    Random random = new Random();

    public Parameterization sampleParameterization() {
        Parameterization parameterization = new CanonicalParameterization();

        int numTypes = this.random.nextInt(2, 11);
        double processLength = this.random.nextDouble(1, 5);

        double[] birthRates = sampleUniformDoubles(numTypes, 1, 3);
        double[] deathRates = Arrays.stream(birthRates).map(x -> x * this.random.nextDouble()).toArray();
        double[] samplingRates = sampleUniformDoubles(numTypes, 0.05, 0.5);
        double[] removalProbabilities = sampleUniformDoubles(numTypes, 0.0, 1.0);
        double[] migrationRates = sampleUniformDoubles(numTypes*(numTypes-1), 0.0, 0.5);

        new TypeSet(numTypes);

        parameterization.initByName(
                "typeSet", new TypeSet(numTypes),
                "processLength", Double.toString(processLength),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter(buildParameterString(birthRates)),
                        numTypes
                ),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter(buildParameterString(deathRates)),
                        numTypes
                ),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter(buildParameterString(samplingRates)),
                        numTypes
                ),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter(buildParameterString(removalProbabilities)),
                        numTypes
                ),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter(buildParameterString(migrationRates)),
                        numTypes
                )
        );

        return parameterization;
    }

    public RealParameter sampleStartTypePriorProbs(Parameterization parameterization) {
        int numTypes = parameterization.getNTypes();

        double[] frequencies = new double[numTypes];
        Arrays.fill(frequencies, 1.0 / numTypes);

        return new RealParameter(buildParameterString(frequencies));
    }

    int sampleNumTypes() {
        return this.random.nextInt(2, 11);
    }

    double[] sampleUniformDoubles(int numValues, double lower, double upper) {
        return this.random.doubles(numValues, lower, upper).toArray();
    }

    int sampleMinIntervals() {
        return 1 ; //Runtime.getRuntime().availableProcessors();
    }

    String buildParameterString(double[] parameters) {
        StringJoiner joiner = new StringJoiner(" ");
        for (double value : parameters) {
            joiner.add(Double.toString(value));
        }
        return joiner.toString();
    }
}
