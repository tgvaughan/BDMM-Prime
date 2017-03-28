package beast.evolution.speciation;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Utils;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.math.ScaledNumbers;
import beast.math.SmallNumber;
import beast.math.SmallNumberScaler;
import beast.math.p0_ODE;
import beast.math.p0ge_InitialConditions;
import beast.math.p0ge_ODE;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created with IntelliJ IDEA.
 * User: Denise
 * Date: 22.08.14
 * Time: 14:05
 */
@Citation("Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. 2016. " +
		"Phylodynamics with migration: \n\t" +
		"A computational framework to quantify population structure from genomic data. \n\t" +
		"Mol Biol Evol. 33(8):2102–2116.")

@Description("Piece-wise constant rates are assumed to be ordered by state and time. First k entries of an array give " +
		"values belonging to type 1, for intervals 1 to k, second k intervals for type 2 etc.")
public abstract class PiecewiseBirthDeathMigrationDistribution extends SpeciesTreeDistribution {


	public Input<RealParameter> frequencies =
			new Input<>("frequencies", "The frequencies for each type",  Input.Validate.REQUIRED);

	public Input<RealParameter> origin =
			new Input<>("origin", "The origin of infection x1");

	public Input<Boolean> originIsRootEdge =
			new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);

	public Input<Integer> maxEvaluations =
			new Input<>("maxEvaluations", "The maximum number of evaluations for ODE solver", 1000000);

	public Input<Boolean> conditionOnSurvival =
			new Input<>("conditionOnSurvival", "condition on at least one survival? Default true.", true);

	public Input<Double> relativeTolerance =
			new Input<>("relTolerance", "relative tolerance for numerical integration", 1e-7);

	public Input<Double> absoluteTolerance =
			new Input<>("absTolerance", "absolute tolerance for numerical integration", 1e-100);

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
			new Input<>("R0", "The basic reproduction number");
	public Input<RealParameter> becomeUninfectiousRate =
			new Input<>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (through recovery or sampling)", Input.Validate.XOR, deathRate);
	public Input<RealParameter> samplingProportion =
			new Input<>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRate);

	public Input<RealParameter> R0_base =
			new Input<>("R0_base",
					"The basic reproduction number for the base pathogen class, should have the same dimension as " +
					"the number of time intervals.");
	public Input<RealParameter> lambda_ratio =
			new Input<>("lambda_ratio",
					"The ratio of basic infection rates of all other classes when compared to the base lambda, " +
					"should have the dimension of the number of pathogens - 1, as it is kept constant over intervals.");

	public Input<RealParameter> migrationMatrix =
			new Input<>("migrationMatrix", "Flattened migration matrix, can be asymmetric, diagnonal entries omitted",  Input.Validate.REQUIRED);

	public Input<RealParameter> birthRateAmongDemes =
			new Input<>("birthRateAmongDemes", "birth rate vector with rate at which transmissions occur among locations");

	public Input<RealParameter> R0AmongDemes =
			new Input<>("R0AmongDemes", "The basic reproduction number determining transmissions occur among locations");


	public Input<RealParameter> removalProbability =
			new Input<RealParameter>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");


	public Input<Integer> stateNumber =
			new Input<>("stateNumber", "The number of states or locations", Input.Validate.REQUIRED);

	public Input<RealParameter> adjustTimesInput =
			new Input<>("adjustTimes", "Origin of MASTER sims which has to be deducted from the change time arrays");
	// <!-- HACK ALERT for reestimation from MASTER sims: adjustTimes is used to correct the forward changetimes such that they don't include orig-root (when we're not estimating the origin) -->

	//TO DO remove one of the inputs
	public Input<Boolean> useRKInput =
			new Input<>("useRK", "Use fixed step size Runge-Kutta with 1000 steps. Default false", false);

	public Input<Boolean> useSmallNumbers = new Input<>("useSN",
			"Use non-underflowing method (default: true)", true);

	public Input<Boolean> checkRho = new Input<>("checkRho", "check if rho is set if multiple tips are given at present (default true)", true);

	//  TO DO CHECKER QUE C'EST PAS POSSIBLE DE REMETTRE 1e-20
	public final static double globalPrecisionThreshold = 1e-10;

	double T;
	double orig;
	int ntaxa;

	p0_ODE P;
	p0ge_ODE PG;

	FirstOrderIntegrator pg_integrator;
	public int maxEvalsUsed;
	public Double minstep;
	public Double maxstep;

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
	Boolean[] isRhoInternalNode;

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

	Double[] freq;

	double[][] pInitialConditions;
	double[] sortedNodes;


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

		if (birthRate.get() == null && R0.get() == null && R0_base.get() == null && lambda_ratio.get() == null) {
			throw new RuntimeException("Either birthRate, R0, or R0_base and R0_ratio need to be specified!");
		} else if ((birthRate.get() != null && R0.get() != null)
				|| (R0.get() != null && (R0_base.get() != null || lambda_ratio.get() != null))
				|| (birthRate.get() != null && (R0_base.get() != null || lambda_ratio.get() != null))) {
			throw new RuntimeException("Only one of birthRate, or R0, or R0_base and lambda_ratio need to be specified!");
		} else if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null) {

			transform = false;
			death = deathRate.get().getValues();
			psi = samplingRate.get().getValues();
			birth = birthRate.get().getValues();
			if (SAModel) r = removalProbability.get().getValues();

			if (birthRateAmongDemes.get()!=null ){

				birthAmongDemes = true;
				b_ij=birthRateAmongDemes.get().getValues();
			}
		} else if ((R0.get() != null || (R0_base.get() != null && lambda_ratio.get() != null)) && becomeUninfectiousRate.get() != null && samplingProportion.get() != null) {
			transform = true;
		} else {
			throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0 (or R0_base AND R0_ratio), becomeUninfectiousRate and samplingProportion!");
		}

		if (transform) {

			if (R0AmongDemes.get()!=null) {
				birthAmongDemes = true;
				b_ij_Changes = R0AmongDemes.get().getDimension()/Math.max(1,(n*(n-1))) - 1;
			}

			if (birthChanges < 1) {
				if (R0.get()!=null) {
					birthChanges = R0.get().getDimension() / n - 1;
				} else {
					birthChanges = R0_base.get().getDimension() - 1;
				}
			}
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

		freq = frequencies.get().getValues();

		double freqSum = 0;
		for (double f : freq) freqSum+= f;
		if (freqSum!=1.)
			throw new RuntimeException("Error: frequencies must add up to 1 but currently add to " + freqSum + ".");

		//        // calculate equilibrium frequencies for 2 types:
		//        double LambMu = -b[0]-b[1]-(d[0]+s[0])+(d[1]+s[1]);
		//        double c = Math.sqrt(Math.pow(LambMu,2) +4*b_ij[0]*b_ij[1]);
		//        freq[0] = (c+LambMu)/(c+LambMu+2*b_ij[0]) ;
		//        freq[1] = 1 - freq[0];

	}

	/**
	 * Perform integration on ODEs g using a classical implementation (using double[] for the initial conditions and the output).
	 * WARNING: getG and getGSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param t
	 * @param PG0
	 * @param t0
	 * @return
	 */
	public double[] getG(double t, double[] PG0, double t0,
			FirstOrderIntegrator pg_integrator, p0ge_ODE PG, Double T, int maxEvalsUsed){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		try {

			if (Math.abs(T-t)<globalPrecisionThreshold || Math.abs(t0-t)<globalPrecisionThreshold ||  T < t) {
				return PG0;
			}

			double from = t;
			double to = t0;
			double oneMinusRho;

			int indexFrom = Utils.index(from, times, times.length);
			int index = Utils.index(to, times, times.length);

			int steps = index - indexFrom;
			if (Math.abs(from-times[indexFrom]) < globalPrecisionThreshold ) steps--;
			if (index>0 && Math.abs(to-times[index-1]) < globalPrecisionThreshold ) {
				steps--;
				index--;
			}

			index--;

			while (steps > 0){

				from = times[index];

				pg_integrator.integrate(PG, to, PG0, from, PG0); // solve PG , store solution in PG0

				if (rhoChanges>0){
					for (int i=0; i<n; i++){
						oneMinusRho = (1-rho[i*totalIntervals + index]);
						PG0[i] *= oneMinusRho;
						PG0[i+n] *= oneMinusRho;
						System.out.println("In getG, multiplying with oneMinusRho: " + oneMinusRho);
					}
				}

				to = times[index];

				steps--;
				index--;
			}

			pg_integrator.integrate(PG, to, PG0, t, PG0); // solve PG , store solution in PG0

		}catch(Exception e){

			throw new RuntimeException("couldn't calculate g");
		}

		if (pg_integrator.getEvaluations() > maxEvalsUsed) maxEvalsUsed = pg_integrator.getEvaluations();

		return PG0;
	}


	/**
	 * Implementation of getG with Small Number structure for the ge equations. Avoids underflowing of integration results.
	 * WARNING: getG and getGSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param t
	 * @param PG0
	 * @param t0
	 * @return
	 */
	public p0ge_InitialConditions getGSmallNumber(double t, p0ge_InitialConditions PG0, double t0,
			FirstOrderIntegrator pg_integrator, p0ge_ODE PG, Double T, int maxEvalsUsed){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		try {

			if (Math.abs(T-t) < globalPrecisionThreshold|| Math.abs(t0-t) < globalPrecisionThreshold ||  T < t) {
				return PG0;
			}

			double from = t;
			double to = t0;
			double oneMinusRho;

			double threshold  = T/10;

			int indexFrom = Utils.index(from, times, times.length);
			int index = Utils.index(to, times, times.length);

			int steps = index - indexFrom;
			if (Math.abs(from-times[indexFrom]) < globalPrecisionThreshold ) steps--;
			if (index>0 && Math.abs(to-times[index-1]) < globalPrecisionThreshold ) {
				steps--;
				index--;
			}
			index--;

			// pgScaled contains the set of initial conditions scaled made to fit the requirements on the values 'double' can represent. It also contains the factor by which the numbers were multiplied
			ScaledNumbers pgScaled = SmallNumberScaler.scale(PG0);
			// integrationResults will temporarily store the results of	 each integration step as 'doubles', before converting them back to 'SmallNumbers'
			double[] integrationResults = new double[2*n];

			while (steps > 0){

				from = times[index];

				if (useRKInput.get() || (to - from) < threshold) {
					pg_integrator.integrate(PG, to, pgScaled.getEquation(), from, integrationResults);
					PG0 = SmallNumberScaler.unscale(integrationResults, pgScaled.getScalingFactor());
				} else {
					pgScaled = safeIntegrate(pg_integrator, PG, to, pgScaled, from); // solve PG , store solution temporarily integrationResults
					// 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
					PG0 = SmallNumberScaler.unscale(pgScaled.getEquation(), pgScaled.getScalingFactor());
				}


				if (rhoChanges>0){
					for (int i=0; i<n; i++){
						oneMinusRho = (1-rho[i*totalIntervals + index]);
						PG0.conditionsOnP[i] *= oneMinusRho;
						PG0.conditionsOnG[i] = PG0.conditionsOnG[i].scalarMultiply(oneMinusRho);
						/*
						System.out.println("In getGSmallNumber, multiplying with oneMinusRho: " + oneMinusRho +", to = " +to);
						 */
					}
				}

				to = times[index];

				steps--;
				index--;

				// 'rescale' the results of the last integration to prepare for the next integration step
				pgScaled = SmallNumberScaler.scale(PG0);
			}

			if (useRKInput.get() || (to - t) < threshold) {
				pg_integrator.integrate(PG, to, pgScaled.getEquation(), t, integrationResults);
				PG0 = SmallNumberScaler.unscale(integrationResults, pgScaled.getScalingFactor());
			} else {
				pgScaled = safeIntegrate(pg_integrator, PG, to, pgScaled, t); // solve PG , store solution temporarily integrationResults
				// 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
				PG0 = SmallNumberScaler.unscale(pgScaled.getEquation(), pgScaled.getScalingFactor());
			}
		}catch(Exception e){

			throw new RuntimeException("couldn't calculate g");
		}

		if (pg_integrator.getEvaluations() > maxEvalsUsed) maxEvalsUsed = pg_integrator.getEvaluations();

		return PG0;
	}

	void setRho(){

		isRhoTip = new Boolean[ treeInput.get().getLeafNodeCount()];
		Arrays.fill(isRhoTip,false);

		isRhoInternalNode = new Boolean[ treeInput.get().getInternalNodeCount()];
		Arrays.fill(isRhoInternalNode,false);

		if (m_rho.get() != null) {

			constantRho = !(m_rho.get().getDimension() > n);

			if (m_rho.get().getDimension() <= n && (rhoSamplingTimes.get()==null || rhoSamplingTimes.get().getDimension() < 2)) {
				if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() <= n && samplingProportion.get().getValue() == 0.) || // todo:  instead of samplingProportion.get().getValue() == 0. need checked that samplingProportion[i]==0 for all i=0..n-1
						(samplingRate.get() != null && samplingRate.get().getDimension() <= 2 && samplingRate.get().getValue() == 0.))) {                              // todo:  instead of samplingRate.get().getValue() == 0. need checked that samplingRate[i]==0 for all i=0..n-1

					// check if data set is contemp!
					for (Node node : treeInput.get().getExternalNodes()){
						if (node.getHeight()>0.) throw new RuntimeException("Error in analysis setup: Parameters set for entirely contemporaneously sampled data, but some nodeheights are > 0!");
					}

					contempData = true;
					System.out.println("BDMM: setting contemp=true.");
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

			if (numChanges > 0 && intervalTimes.getDimension() != numChanges + 1) {
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


	void updateBirthDeathPsiParams(){

		Double[] birthRates = birthRate.get().getValues();
		Double[] deathRates = deathRate.get().getValues();
		Double[] samplingRates = samplingRate.get().getValues();
		Double[] removalProbabilities = new Double[1];

		if (SAModel) {
			removalProbabilities = removalProbability.get().getValues();
			r =  new Double[n*totalIntervals];
		}

		int state;

		for (int i = 0; i < n*totalIntervals; i++) {

			state =  i/totalIntervals;

			birth[i] = birthRates[birthRates.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state];
			death[i] = deathRates[deathRates.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state];
			psi[i] = samplingRates[samplingRates.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state];
			if (SAModel) r[i] = removalProbabilities[removalProbabilities.length > n ? (rChanges+1)*state+index(times[i%totalIntervals], rChangeTimes) : state];

		}

	}


	void updateAmongParameter(Double[] param, Double[] paramFrom, int nrChanges, List<Double> changeTimes){

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int dt = 0; dt < totalIntervals; dt++) {
					if (i != j) {
						param[(i * (n - 1) + (j < i ? j : j - 1)) * totalIntervals + dt]
								= paramFrom[(paramFrom.length > (n * (n - 1)))
								            ? (nrChanges + 1) * (n - 1) * i + index(times[dt], changeTimes)
								            : (i * (n - 1) + (j < i ? j : j - 1))];
					}
				}
			}
		}

	}


	void updateRho(){
		if (m_rho.get() != null && (m_rho.get().getDimension()==1 ||  rhoSamplingTimes.get() != null)) {

			Double[] rhos = m_rho.get().getValues();
			rho = new Double[n*totalIntervals];
			int state;

			for (int i = 0; i < totalIntervals*n; i++) {

				state =  i/totalIntervals;

				rho[i]= rhoChanges>0?
						rhoSamplingChangeTimes.contains(times[i]) ? rhos[rhos.length > n ? (rhoChanges+1)*state+index(times[i%totalIntervals], rhoSamplingChangeTimes) : state] : 0.
								: rhos[0];
			}
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


	public void transformWithinParameters(){

		Double[] p = samplingProportion.get().getValues();
		Double[] ds = becomeUninfectiousRate.get().getValues();
		Double[] R;
		if (R0.get() != null) {
			R = R0.get().getValues();
		} else {
			Double[] l_ratio = lambda_ratio.get().getValues();
			Double[] R_sens = R0_base.get().getValues();

			int totalIntervals = R_sens.length;
			int totalTypes = l_ratio.length + 1;
			R = new Double[totalIntervals * totalTypes];
			for (int i=0; i < totalIntervals; i++) {
				R[i] = R_sens[i];
				for (int j=1; j < totalTypes; j++) {
					double lambda = R_sens[i] * ds[ds.length > totalTypes ? index(times[i%totalIntervals], deathRateChangeTimes) : 0];
					R[i + totalIntervals * j] = (lambda * l_ratio[j - 1]) / ds[ds.length > totalTypes ? (deathChanges+1)*j+index(times[i%totalIntervals], deathRateChangeTimes) : j];
				}
			}
		}

		Double[] removalProbabilities = new Double[1];
		if (SAModel) removalProbabilities = removalProbability.get().getValues();

		int state;

		for (int i = 0; i < totalIntervals*n; i++){

			state =  i/totalIntervals;

			birth[i] = R[R.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state]
					* ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state];

			if (!SAModel) {
				psi[i] = p[p.length > n ? (samplingChanges + 1) * state + index(times[i % totalIntervals], samplingRateChangeTimes) : state]
						* ds[ds.length > n ? (deathChanges + 1) * state + index(times[i % totalIntervals], deathRateChangeTimes) : state];

				death[i] = ds[ds.length > n ? (deathChanges + 1) * state + index(times[i % totalIntervals], deathRateChangeTimes) : state] - psi[i];
			}

			else {
				r[i] = removalProbabilities[removalProbabilities.length > n ? (rChanges+1)*state+index(times[i%totalIntervals], rChangeTimes) : state];

				psi[i] = p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]
						* ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state]
								/ (1+(r[i]-1)*p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]);


				death[i] = ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] - psi[i]*r[i];
			}
		}

	}


	public void transformAmongParameters(){

		Double[] RaD = (birthAmongDemes) ? R0AmongDemes.get().getValues() : new Double[1];
		Double[] ds = becomeUninfectiousRate.get().getValues();

		if (birthAmongDemes)    {

			for (int i = 0; i < n; i++){

				for (int j=0; j<n ; j++){

					for (int dt=0; dt<totalIntervals; dt++){

						if (i!=j){
							b_ij[(i*(n-1)+(j<i?j:j-1))*totalIntervals+dt]
									= RaD[(RaD.length>(n*(n-1)))
									      ?  (b_ij_Changes+1)*(n-1)*i + index(times[dt], b_ijChangeTimes)
									      : (i*(n-1)+(j<i?j:j-1))]
									    		  * ds[ds.length > n ? (deathChanges+1)*i+index(times[dt], deathRateChangeTimes) : i];
						}
					}
				}

			}
		}
	}


	void checkOrigin(TreeInterface tree){

		if (origin.get()==null){
			T = tree.getRoot().getHeight();
		}
		else {

			updateOrigin(tree.getRoot());

			if (!Boolean.valueOf(System.getProperty("beast.resume")) && orig < 0)
				throw new RuntimeException("Error: origin("+T+") must be larger than tree height("+tree.getRoot().getHeight()+")!");
		}

	}


	void updateOrigin(Node root){

		T = origin.get().getValue();
		orig = T - root.getHeight();

		if (originIsRootEdge.get()) {

			orig = origin.get().getValue();
			T = orig + root.getHeight();
		}

	}


	void setupIntegrators(){   // set up ODE's and integrators

		if (minstep == null) minstep = T*1e-100;
		if (maxstep == null) maxstep = T/10;

		Boolean augmented = this instanceof BirthDeathMigrationModel;

		P = new p0_ODE(birth, ((!augmented && birthAmongDemes) ? b_ij : null), death,psi,M, n, totalIntervals, times);
		PG = new p0ge_ODE(birth, ((!augmented && birthAmongDemes) ? b_ij : null), death,psi,M, n, totalIntervals, T, times, P, maxEvaluations.get(), augmented);

		p0ge_ODE.globalPrecisionThreshold = globalPrecisionThreshold;

		if (!useRKInput.get() && useSmallNumbers.get()) {
			pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get()); //new HighamHall54Integrator(minstep, maxstep, absolutePrecision.get(), tolerance.get()); //new DormandPrince853Integrator(minstep, maxstep, absolutePrecision.get(), tolerance.get()); //new DormandPrince54Integrator(minstep, maxstep, absolutePrecision.get(), tolerance.get()); // 
			pg_integrator.setMaxEvaluations(maxEvaluations.get());

			PG.p_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get()); 
			PG.p_integrator.setMaxEvaluations(maxEvaluations.get());
		} else {
			pg_integrator = new ClassicalRungeKuttaIntegrator(T / 1000);
			PG.p_integrator = new ClassicalRungeKuttaIntegrator(T / 1000);

		}
	}



	/**
	 * Obtain element of rate matrix for migration model for use in likelihood
	 * calculation.
	 *
	 * @param i
	 * @param j
	 * @return Rate matrix element.
	 */
	public double getNbyNRate(int i, int j) {
		if (i==j)
			return 0;

		int offset = getArrayOffset(i, j);
		return migrationMatrix.get().getValue(offset);
	}


	/**
	 * Obtain offset into "rate matrix" and associated flag arrays.
	 *
	 * @param i
	 * @param j
	 * @return Offset (or -1 if i==j)
	 */
	protected int getArrayOffset(int i, int j) {

		if (i==j)
			throw new RuntimeException("Programmer error: requested migration "
					+ "rate array offset for diagonal element of "
					+ "migration rate matrix.");


		if (j>i)
			j -= 1;
		return i*(n-1)+j;   // todo: check if this is correct!!!
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

	/**
	 * If integration interval is too long to provide precise results, cuts it in half and starts integration again.
	 * @param integrator
	 * @param PG
	 * @param to
	 * @param pgScaled
	 * @param from
	 * @return
	 */
	public ScaledNumbers safeIntegrate(FirstOrderIntegrator integrator, p0ge_ODE PG, double to, ScaledNumbers pgScaled, double from){

		// if the integration interval is too small, nothing is done (to prevent infinite looping)
		if(Math.abs(from-to)< (T * 1e-10)) return pgScaled;

		// we test to see if the current interval size can produce a sufficiently-precise result
		double[] integrationResults = new double[pgScaled.getEquation().length];
		integrator.integrate(PG, to, pgScaled.getEquation(), from, integrationResults);
		// look for the minimal value
		double minRes = getMinIntegrationOnGe(integrationResults);

		if(absoluteTolerance.get() > relativeTolerance.get())
			throw new RuntimeException("Absolute tolerance higher than relative tolerance for the adaptive integrator. Change values for these inputs.");


		if(minRes != Double.MAX_VALUE && (minRes<0 || absoluteTolerance.get()/minRes > relativeTolerance.get())) {

			pgScaled = safeIntegrate(integrator, PG, to, pgScaled, from + (to-from)/2);
			pgScaled = safeIntegrate(integrator, PG, from + (to-from)/2, pgScaled, from);

		} else {

			int a = pgScaled.getScalingFactor();
			int n = integrationResults.length/2;
			double[] pConditions = new double[n];
			SmallNumber[] geConditions = new SmallNumber[n];

			for (int i = 0; i < n; i++) {
				pConditions[i] = integrationResults[i];
				geConditions[i] = new SmallNumber(integrationResults[i+n]);
			}
			pgScaled = SmallNumberScaler.scale(new p0ge_InitialConditions(pConditions, geConditions));
			pgScaled.augmentFactor(a);

		}

		return pgScaled;
	}

	/**
	 * Find the lowest non-zero value among the integration results, only on the ge part.
	 * @param values
	 * @return
	 */
	public double getMinIntegrationOnGe(double[] values) {
		double min=Double.MAX_VALUE;
		if (values.length < 2)
			throw new RuntimeException("Invalid inital-conditions array size");
		for (int i = (values.length -1); i> (values.length/2 - 1) ; i-- ) {
			if (values[i] < min && values[i]!=0 ) min=values[i];
		}
		return min;
	}


	public p0ge_ODE getPG() {
		return PG;
	}

	/*
	 * Find all initial conditions for all future integrations on p0 equations 
	 * @param tree
	 * @return an array of arrays storing the initial conditions values
	 */
	public double[][] getAllInitialConditionsForP(TreeInterface tree){
		int leafCount = tree.getLeafNodeCount();
		double[] leafHeights = new double[leafCount];
		int[] indicesSortedByLeafHeight  =new int[leafCount];
		for (int i=0; i<leafCount; i++){
			leafHeights[i] = T - tree.getNode(i).getHeight();
			// System.out.println(nodeHeight[i]);
			indicesSortedByLeafHeight[i] = i;
		}

		HeapSort.sort(leafHeights, indicesSortedByLeafHeight);
		//"sort" sorts in ascending order, so we have to be careful since the integration starts from the leaves at height T and goes up to the root at height 0 (or >0)

		double[][] pInitialCondsAtLeaves = new double[leafCount + 1][n]; 

		double t = leafHeights[indicesSortedByLeafHeight[leafCount-1]];

		boolean rhoSampling =  (m_rho.get()!=null);

		pInitialCondsAtLeaves[indicesSortedByLeafHeight[leafCount-1]] = PG.getP(t, rhoSampling, rho);
		double t0 = t;

		if (leafCount >1 ){
			for (int i = leafCount-2; i>-1; i--){
				t = leafHeights[indicesSortedByLeafHeight[i]];

				//If the next higher leaf is actually at the same height, store previous results and skip iteration
				if (Math.abs(t-t0) < globalPrecisionThreshold) {
					t0 = t;
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] = pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]];
					continue;
				} else {
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] = PG.getP(t, pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]], t0, rhoSampling, rho);
					t0 = t;
				}

			}
		}


		pInitialCondsAtLeaves[leafCount] = PG.getP(0, pInitialCondsAtLeaves[indicesSortedByLeafHeight[0]], t0, rhoSampling, rho);

		return pInitialCondsAtLeaves;
	}

}
