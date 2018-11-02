package bdmm.distributions;

import beast.core.Input;

import java.util.List;

public class CanonicalParameterization extends Parameterization {

    public Input<Skyline> migRateInput = new Input<>("migrationRate",
            "Migration rate skyline.", Input.Validate.REQUIRED);

    public Input<Skyline> birthRateInput = new Input<>("birthRate",
            "Birth rate skyline.", Input.Validate.REQUIRED);

    public Input<Skyline> crossBirthRateInput = new Input<>("birthRateAmongDemes",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<Skyline> deathRateInput = new Input<>("deathRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<Skyline> samplingRateInput = new Input<>("samplingRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<Skyline> removalRateInput = new Input<>("removalRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<Skyline> rhoSamplingInput = new Input<>("rhoSampling",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);


    @Override
    public List<Double> getMigRateChangeTimes(double maxTime) {
        return migRateInput.get().getSanitizedChangeTimes(maxTime);
    }

    @Override
    public List<Double> getBirthRateChangeTimes(double maxTime) {
        return birthRateInput.get().getSanitizedChangeTimes(maxTime);
    }

    @Override
    public List<Double> getCrossBirthRateChangeTimes(double maxTime) {
        return crossBirthRateInput.get().getSanitizedChangeTimes(maxTime);
    }

    @Override
    public List<Double> getDeathRateChangeTimes(double maxTime) {
        return deathRateInput.get().getSanitizedChangeTimes(maxTime);
    }

    @Override
    public List<Double> getSamplingRateChangeTimes(double maxTime) {
        return samplingRateInput.get().getSanitizedChangeTimes(maxTime);
    }

    @Override
    public List<Double> getRemovalProbChangeTimes(double maxTime) {
        return removalRateInput.get().getSanitizedChangeTimes(maxTime);
    }

    @Override
    public List<Double> getRhoSamplingTimes(double maxTime) {
        return rhoSamplingInput.get().getSanitizedChangeTimes(maxTime);
    }
}
