package beast.evolution.speciation;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Denise
 * Date: 22.08.14
 * Time: 14:05
 */
@Description("Piece-wise constant rates are assumed to be ordered by state and time. First k entries of an array give " +
        "values belonging to type 1, for intervals 1 to k, second k intervals for type 2 etc.")
public abstract class PiecewiseBirthDeathSamplingDistribution extends SpeciesTreeDistribution {


    // the interval times for the migration rates
    public Input<RealParameter> migChangeTimesInput =
            new Input<>("migChangeTimes", "The times t_i specifying when migration rate changes occur", (RealParameter) null);

    // the interval times for the birth rate
    public Input<RealParameter> birthRateChangeTimesInput =
            new Input<>("birthRateChangeTimes", "The times t_i specifying when birth/R rate changes occur", (RealParameter) null);

    // the interval times for the birth rate among demes
    public Input<RealParameter> b_ijChangeTimesInput =
            new Input<>("birthRateAmongDemesChangeTimes", "The times t_i specifying when birth/R among demes changes occur", (RealParameter) null);

    // the interval times for the death rate
    public Input<RealParameter> deathRateChangeTimesInput =
            new Input<>("deathRateChangeTimes", "The times t_i specifying when death/becomeUninfectious rate changes occur", (RealParameter) null);

    // the interval times for sampling rate
    public Input<RealParameter> samplingRateChangeTimesInput =
            new Input<>("samplingRateChangeTimes", "The times t_i specifying when sampling rate or sampling proportion changes occur", (RealParameter) null);

    // the interval times for removal probability
    public Input<RealParameter> removalProbabilityChangeTimesInput =
            new Input<RealParameter>("removalProbabilityChangeTimes", "The times t_i specifying when removal probability changes occur", (RealParameter) null);

    public Input<RealParameter> intervalTimes =
            new Input<>("intervalTimes", "The time t_i for all parameters if they are the same", (RealParameter) null);

    public Input<Boolean> migTimesRelativeInput =
            new Input<>("migTimesRelative", "True if migration rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> b_ijChangeTimesRelativeInput =
            new Input<>("birthRateTimesRelative", "True if birth rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> birthRateChangeTimesRelativeInput =
            new Input<>("birthRateAmongDemesTimesRelative", "True if birth rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> deathRateChangeTimesRelativeInput =
            new Input<>("deathRateTimesRelative", "True if death rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> samplingRateChangeTimesRelativeInput =
            new Input<>("samplingRateTimesRelative", "True if sampling rate times specified relative to tree height? Default false", false);

    Input<Boolean> removalProbabilityChangeTimesRelativeInput =
            new Input<Boolean>("removalProbabilityTimesRelative", "True if removal probability change times specified relative to tree height? Default false", false);

    public Input<BooleanParameter> reverseTimeArraysInput =
            new Input<>("reverseTimeArrays", "True if the time arrays are given in backwards time (from the present back to root). Order: 1) birth 2) death 3) sampling 4) rho 5) r 6) migration. Default false." +
                    "Careful, rate array must still be given in FORWARD time (root to tips).");

    // the times for rho sampling
    public Input<RealParameter> rhoSamplingTimes =
            new Input<>("rhoSamplingTimes", "The times t_i specifying when rho-sampling occurs", (RealParameter) null);
    public Input<Boolean> contemp =
            new Input<>("contemp", "Only contemporaneous sampling (i.e. all tips are from same sampling time, default false)", false);

    public Input<RealParameter> birthRate =
            new Input<>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time");
    public Input<RealParameter> deathRate =
            new Input<>("deathRate", "The deathRate vector with birthRates between times");
    public Input<RealParameter> samplingRate =
            new Input<>("samplingRate", "The sampling rate per individual");      // psi

    public Input<RealParameter> m_rho =
            new Input<>("rho", "The proportion of lineages sampled at rho-sampling times (default 0.)");


    public Input<RealParameter> R0 =
            new Input<>("R0", "The basic reproduction number", Input.Validate.XOR, birthRate);
    public Input<RealParameter> becomeUninfectiousRate =
            new Input<>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (throuch recovery or sampling)", Input.Validate.XOR, deathRate);
    public Input<RealParameter> samplingProportion =
            new Input<>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRate);


    public Input<RealParameter> migrationMatrix =
            new Input<>("migrationMatrix", "Flattened migration matrix, can be asymmetric, diagnonal entries omitted",  Input.Validate.REQUIRED);

    public Input<RealParameter> birthRateAmongDemes =
             new Input<>("birthRateAmongDemes", "birth rate vector with rate at which transmissions occur among locations");

    public Input<RealParameter> R0AmongDemes =
            new Input<>("R0AmongDemes", "The basic reproduction number determining transmissions occur among locations");


    public Input<RealParameter> removalProbability =
            new Input<RealParameter>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");


    //coupling R0 changes: assume there are 2 types with R0 = [r1,r2], R0AmongDemes = [r12,r21] and coupledR0Changes=[c1a,c1b,c2a,c2b],
    //                     i.e. there are dim(coupledR0Changes)/(#types) = 2 rate changes through time for the R0's
    // this translates to   R0 = [r1,r1*c1a,r1*c1b,r2,r2*c2a,r2*c2b]
    //    and               R0AmongDemes = [r12,r12*c2a,r12*c2b,r21,r21*c1a,r21*c1b]  // note that the scale of change is determined by the "receiving" deme here
    public Input<RealParameter> coupledR0Changes =
            new Input<>("coupledR0Changes", "The scale of change in R0 and R0AmongDemes per interval, when they are assumed to be equal for both");


    public Input<Integer> stateNumber =
            new Input<>("stateNumber", "The number of states or locations", Input.Validate.REQUIRED);

    public Input<RealParameter> adjustTimesInput =
            new Input<>("adjustTimes", "Origin of MASTER sims which has to be deducted from the change time arrays");
   // <!-- HACK ALERT for reestimation from MASTER sims: adjustTimes is used to correct the forward changetimes such that they don't include orig-root (when we're not estimating the origin) -->

    public Input<Boolean> useRKInput =
            new Input<>("useRK", "Use fixed step size Runge-Kutta with 1000 steps. Default true", true);


    // these four arrays are totalIntervals in length
    protected Double[] birth;
    Double[] death;
    Double[] psi;
    Double[] rho;
    Double[] r;

    /**
     * The number of change points in the birth rate, b_ij, death rate, sampling rate, rho, r
     */
    int migChanges;
    int birthChanges;
    int b_ij_Changes;
    int deathChanges;
    int samplingChanges;
    int rhoChanges;
    int rChanges;


    public boolean SAModel;

    /**
     * The number of times rho-sampling occurs
     */
    int rhoSamplingCount;
    Boolean constantRho;
    Boolean[] isRhoTip;

    /**
     * Total interval count
     */
    int totalIntervals;
    int n;  // number of states / locations

    protected List<Double> migChangeTimes = new ArrayList<>();
    protected List<Double> birthRateChangeTimes = new ArrayList<>();
    protected List<Double> b_ijChangeTimes = new ArrayList<>();
    protected List<Double> deathRateChangeTimes = new ArrayList<>();
    protected List<Double> samplingRateChangeTimes = new ArrayList<>();
    protected List<Double> rhoSamplingChangeTimes = new ArrayList<>();
    protected List<Double> rChangeTimes = new ArrayList<Double>();

    Boolean contempData;
    SortedSet<Double> timesSet = new TreeSet<>();

    protected Double[] times = new Double[]{0.};

    protected Boolean transform;

    Boolean migTimesRelative = false;
    Boolean birthRateTimesRelative = false;
    Boolean b_ijTimesRelative = false;
    Boolean deathRateTimesRelative = false;
    Boolean samplingRateTimesRelative = false;
    Boolean rTimesRelative = false;
    Boolean[] reverseTimeArrays;

    Double[] M;
    Double[] b_ij;
    Boolean birthAmongDemes = false;


    @Override
    public void initAndValidate() {

        if (removalProbability.get() != null) SAModel = true;

        birth = null;
        b_ij = null;
        death = null;
        psi = null;
        rho = null;
        r = null;
        birthRateChangeTimes.clear();
        deathRateChangeTimes.clear();
        samplingRateChangeTimes.clear();
        if (SAModel) rChangeTimes.clear();
        totalIntervals = 0;
        n = stateNumber.get();

        M = migrationMatrix.get().getValues();

        if (n>1 && M.length != n*(n-1))
            if (migChangeTimesInput.get()==null || M.length != n*(n-1) *migChangeTimesInput.get().getDimension())
                throw new RuntimeException("Migration matrix dimension is incorrect!");

        birthRateTimesRelative = birthRateChangeTimesRelativeInput.get();
        b_ijTimesRelative = b_ijChangeTimesRelativeInput.get();
        migTimesRelative = migTimesRelativeInput.get();
        deathRateTimesRelative = deathRateChangeTimesRelativeInput.get();
        samplingRateTimesRelative = samplingRateChangeTimesRelativeInput.get();
        if (SAModel) rTimesRelative = removalProbabilityChangeTimesRelativeInput.get();

        reverseTimeArrays = new Boolean[]{false, false, false, false, false, false};
        if (reverseTimeArraysInput.get()!= null )  {
            Boolean[] r = reverseTimeArraysInput.get().getValues();
            for (int i=0; i<r.length; i++)
                reverseTimeArrays[i] = r[i];
        }

        rhoSamplingCount = 0;
        contempData = contemp.get();


        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null) {

            transform = false;
            death = deathRate.get().getValues();
            psi = samplingRate.get().getValues();
            birth = birthRate.get().getValues();
            if (SAModel) r = removalProbability.get().getValues();

            if (birthRateAmongDemes.get()!=null ){

                birthAmongDemes = true;
                b_ij=birthRateAmongDemes.get().getValues();
            }

        } else if (R0.get() != null && becomeUninfectiousRate.get() != null && samplingProportion.get() != null) {

            transform = true;

        } else {
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }

        if (transform) {

            if (R0AmongDemes.get()!=null) {
                birthAmongDemes = true;
                b_ij_Changes = R0AmongDemes.get().getDimension()/Math.max(1,(n*(n-1))) - 1;
            }

            if (birthChanges < 1) birthChanges = R0.get().getDimension()/n - 1;
            samplingChanges = samplingProportion.get().getDimension()/n - 1;
            deathChanges = becomeUninfectiousRate.get().getDimension()/n - 1;

        } else {    //todo: b d s param doesn't work yet with rate changes (unless all parameters have equally many)

            if (birthChanges < 1) birthChanges = birthRate.get().getDimension()/n - 1;
            if (birthAmongDemes) b_ij_Changes = birthRateAmongDemes.get().getDimension()/(n*(n-1)) - 1;
            deathChanges = deathRate.get().getDimension()/n - 1;
            samplingChanges = samplingRate.get().getDimension()/n - 1;
        }

        migChanges = migrationMatrix.get().getDimension()/Math.max(1,(n*(n-1))) - 1;

        if (SAModel) rChanges = removalProbability.get().getDimension()/n -1;

        if (m_rho.get()!=null) {
            rho = m_rho.get().getValues();
            rhoChanges = m_rho.get().getDimension()/n - 1;
        }

        if (coupledR0Changes.get()!=null){

            if (this instanceof BirthDeathMigrationModel)  throw new RuntimeException("Error: Coupled R0 changes are not implemented for multitype trees!");

            if ((birthChanges>0 && birthChanges!=n) ||  (b_ij_Changes>0 && b_ij_Changes!=n*(n-1))) throw new RuntimeException("if coupledR0Changes!=null R0 and R0AmongDemes must be of dimension 1");

            birthChanges = coupledR0Changes.get().getDimension()/n;
            b_ij_Changes = coupledR0Changes.get().getDimension()/n;
        }
//          collectTimes(T); //  these need to be called from implementing initAndValidate
//          setRho();

    }

    void setRho(){

        isRhoTip = new Boolean[ treeInput.get().getLeafNodeCount()];
        Arrays.fill(isRhoTip,false);

        if (m_rho.get() != null) {

            constantRho = !(m_rho.get().getDimension() > n);

            if (m_rho.get().getDimension() <= n && (rhoSamplingTimes.get()==null || rhoSamplingTimes.get().getDimension() < 2)) {
                if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() <= n && samplingProportion.get().getValue() == 0.) || // todo:  instead of samplingProportion.get().getValue() == 0. need checked that samplingProportion[i]==0 for all i=0..n-1
                        (samplingRate.get() != null && samplingRate.get().getDimension() <= 2 && samplingRate.get().getValue() == 0.))) {                              // todo:  instead of samplingRate.get().getValue() == 0. need checked that samplingRate[i]==0 for all i=0..n-1
                    contempData = true;
                }
            }

            if (contempData) {
                if (m_rho.get().getDimension() != 1 && m_rho.get().getDimension() != n)
                    throw new RuntimeException("when contemp=true, rho must have dimension 1 (or equal to the stateNumber)");

                else {
                    rho = new Double[n*totalIntervals];
                    Arrays.fill(rho, 0.);
                    Arrays.fill(isRhoTip, true);
                    for (int i=1; i<=n; i++)  rho[i*totalIntervals - 1] = m_rho.get().getValue(i-1);

                    rhoSamplingCount = 1;
                }
            }
            else {
                Double[] rhos = m_rho.get().getValues();
                rho = new Double[n*totalIntervals];
                Arrays.fill(rho, 0.);
                for (int i = 0; i < totalIntervals; i++) {
                    for (int j=0;j<n;j++){
                       rho[j*totalIntervals+i]= rhoSamplingChangeTimes.contains(times[i]) ? (rhos[constantRho? j : j*(1+rhoChanges)+rhoSamplingChangeTimes.indexOf(times[i])]) : 0.;
                    }
                }
                computeRhoTips();
            }

        } else {
            rho = new Double[n*totalIntervals];
            Arrays.fill(rho, 0.);
        }

    }

    abstract void computeRhoTips();

    /**
     * Collect all the times of parameter value changes and rho-sampling events
     */
    void collectTimes(double maxTime) {

        timesSet.clear();

        getChangeTimes(maxTime, migChangeTimes,
                migChangeTimesInput.get() != null ? migChangeTimesInput.get() : intervalTimes.get(),
                migChanges, migTimesRelative, reverseTimeArrays[5]);

        getChangeTimes(maxTime, birthRateChangeTimes,
                birthRateChangeTimesInput.get() != null ? birthRateChangeTimesInput.get() : intervalTimes.get(),
                birthChanges, birthRateTimesRelative, reverseTimeArrays[0]);

       getChangeTimes(maxTime, b_ijChangeTimes,
                b_ijChangeTimesInput.get() != null ? b_ijChangeTimesInput.get() : intervalTimes.get(),
                b_ij_Changes, b_ijTimesRelative, reverseTimeArrays[0]);

        getChangeTimes(maxTime, deathRateChangeTimes,
                deathRateChangeTimesInput.get() != null ? deathRateChangeTimesInput.get() : intervalTimes.get(),
                deathChanges, deathRateTimesRelative, reverseTimeArrays[1]);

        getChangeTimes(maxTime, samplingRateChangeTimes,
                samplingRateChangeTimesInput.get() != null ? samplingRateChangeTimesInput.get() : intervalTimes.get(),
                samplingChanges, samplingRateTimesRelative, reverseTimeArrays[2]);

        getChangeTimes(maxTime, rhoSamplingChangeTimes,
                rhoSamplingTimes.get()!=null ? rhoSamplingTimes.get() : intervalTimes.get(),
                rhoChanges, false, reverseTimeArrays[3]);

        if (SAModel) getChangeTimes(maxTime, rChangeTimes,
                removalProbabilityChangeTimesInput.get() != null ? removalProbabilityChangeTimesInput.get() : intervalTimes.get(),
                rChanges, rTimesRelative, reverseTimeArrays[4]);

        for (Double time : migChangeTimes) {
            timesSet.add(time);
        }

        for (Double time : birthRateChangeTimes) {
            timesSet.add(time);
        }

        for (Double time : b_ijChangeTimes) {
             timesSet.add(time);
        }

         for (Double time : deathRateChangeTimes) {
            timesSet.add(time);
        }

        for (Double time : samplingRateChangeTimes) {
            timesSet.add(time);
        }

        for (Double time : rhoSamplingChangeTimes) {
            timesSet.add(time);
        }

        if (SAModel) {
            for (Double time : rChangeTimes) {
                timesSet.add(time);
            }
        }


        times = timesSet.toArray(new Double[timesSet.size()]);
        totalIntervals = times.length;

    }

    /**
     * set change times
     */
    public void getChangeTimes(double maxTime, List<Double> changeTimes, RealParameter intervalTimes, int numChanges, boolean relative, boolean reverse) {
        changeTimes.clear();

        if (intervalTimes == null) { //equidistant

            double intervalWidth = maxTime / (numChanges + 1);

            double end;
            for (int i = 1; i <= numChanges; i++) {
                end = (intervalWidth) * i;
                changeTimes.add(end);
            }
            end = maxTime;
            changeTimes.add(end);

        } else {

            if (!reverse && intervalTimes.getValue(0) != 0.0) {
                throw new RuntimeException("First time in interval times parameter should always be zero.");
            }

            if (numChanges > 0 && coupledR0Changes.get()==null && intervalTimes.getDimension() != numChanges + 1) {
                throw new RuntimeException("The time interval parameter should be numChanges + 1 long (" + (numChanges + 1) + ").");
            }

            int dim = intervalTimes.getDimension();

            double end;
            for (int i = (reverse?0:1); i < dim; i++) {
                end = reverse ? (maxTime - intervalTimes.getValue(dim - i - 1)) : intervalTimes.getValue(i);
                if (relative) end *= maxTime;
                if (end != maxTime) changeTimes.add(end);
            }
            end = maxTime;

            if (adjustTimesInput.get()!=null){

                double iTime;
                double aTime = adjustTimesInput.get().getValue();

                for (int i = 0 ; i < changeTimes.size(); i++){

                    iTime = intervalTimes.getArrayValue(i+1);

                    changeTimes.set(i, Math.abs(end-aTime+iTime) );
                }
            }

            changeTimes.add(end);
        }
    }

    /**
     * @param t the time in question
     * @return the index of the given time in the list of times, or if the time is not in the list, the index of the
     *         next smallest time
     *         This index function should only be used in transformParameters(), for likelihood calculations the times List needs to be used (e.g. with Untils.index(...))
     */
    public int index(double t, List<Double> times) {

        int epoch = Collections.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return epoch;
    }

    	// Interface requirements:

	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

    @Override
    public boolean requiresRecalculation(){
        return true;
    }

}
