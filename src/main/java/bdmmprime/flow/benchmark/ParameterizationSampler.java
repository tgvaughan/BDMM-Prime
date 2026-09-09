package bdmmprime.flow.benchmark;

import bdmmprime.parameterization.*;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.Simplex;

import java.util.Arrays;
import java.util.Random;

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
                        new RealVectorParam<>(birthRates, Real.INSTANCE),
                        numTypes
                ),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(deathRates, Real.INSTANCE),
                        numTypes
                ),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(samplingRates, Real.INSTANCE),
                        numTypes
                ),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(removalProbabilities, Real.INSTANCE),
                        numTypes
                ),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(migrationRates, Real.INSTANCE),
                        numTypes
                )
        );

        return parameterization;
    }

    public Simplex sampleStartTypePriorProbs(Parameterization parameterization) {
        int numTypes = parameterization.getNTypes();

        double[] frequencies = new double[numTypes];
        Arrays.fill(frequencies, 1.0 / numTypes);

        return new SimplexParam(frequencies);
    }

    double[] sampleUniformDoubles(int numValues, double lower, double upper) {
        return this.random.doubles(numValues, lower, upper).toArray();
    }

    int sampleMinIntervals() {
        return 1 ; //Runtime.getRuntime().availableProcessors();
    }

}
