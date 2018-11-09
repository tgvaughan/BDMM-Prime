package bdmm.distributions;

import beast.core.*;
import beast.core.parameter.RealParameter;
import bdmm.util.Utils;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.HeapSort;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;


/**
 * Created with IntelliJ IDEA.
 * User: Denise
 * Date: 22.08.14
 * Time: 14:05
 */
@Citation(value="Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. (2016). " +
		"Phylodynamics with migration: \n" +
		"A computational framework to quantify population structure from genomic data. \n" +
		"Mol Biol Evol. 33(8):2102–2116."
		, DOI= "10.1093/molbev/msw064", year = 2016, firstAuthorSurname = "Kuehnert")

@Description("Piece-wise constant rates are assumed to be ordered by state and time. First k entries of an array give " +
		"values belonging to type 1, for intervals 1 to k, second k intervals for type 2 etc.")
public abstract class PiecewiseBirthDeathMigrationDistribution extends SpeciesTreeDistribution {


    public Input<Parameterization> parameterizationInput =
            new Input<>("parameterization", "BDMM parameterization", Input.Validate.REQUIRED);

	public Input<RealParameter> frequencies =
			new Input<>("frequencies", "The frequencies for each type",  Input.Validate.REQUIRED);

	public Input<Integer> maxEvaluations =
			new Input<>("maxEvaluations", "The maximum number of evaluations for ODE solver", 1000000);

	public Input<Boolean> conditionOnSurvival =
			new Input<>("conditionOnSurvival", "condition on at least one survival? Default true.", true);

	public Input<Double> relativeTolerance =
			new Input<>("relTolerance", "relative tolerance for numerical integration", 1e-7);

	public Input<Double> absoluteTolerance =
			new Input<>("absTolerance", "absolute tolerance for numerical integration", 1e-100 /*Double.MIN_VALUE*/);

	public Input<Boolean> checkRho = new Input<>("checkRho", "check if rho is set if multiple tips are given at present (default true)", true);


	public Input<Boolean> isParallelizedCalculationInput = new Input<>("parallelize", "is the calculation parallelized on sibling subtrees or not (default true)", true);

	//If a large number a cores is available (more than 8 or 10) the calculation speed can be increased by diminishing the parallelization factor
	//On the contrary, if only 2-4 cores are available, a slightly higher value (1/5 to 1/8) can be beneficial to the calculation speed.
	public Input<Double> minimalProportionForParallelizationInput = new Input<>("parallelizationFactor", "the minimal relative size the two children subtrees of a node" +
			" must have to start parallel calculations on the children. (default: 1/10). ", 1.0/10);

	Parameterization parameterization;

	public static boolean isParallelizedCalculation;

	public static double minimalProportionForParallelization;

	//  TODO check if it's possible to have 1e-20 there
	public final static double globalPrecisionThreshold = 1e-10;

	int ntaxa;

	p0_ODE P;
	p0ge_ODE PG;

	FirstOrderIntegrator pg_integrator;
	public static Double minstep;
	public static Double maxstep;

	/**
	 * Total interval count
	 */
	static int n;  // number of states / locations

	Double[] freq;

	static double[][] pInitialConditions;

	public double[] weightOfNodeSubTree;

	double parallelizationThreshold;

	static ExecutorService executor;
	static ThreadPoolExecutor pool;


	TreeInterface tree;

	@Override
	public void initAndValidate() {
	    parameterization = parameterizationInput.get();

		tree = treeInput.get();

		Double factor;

		freq = frequencies.get().getValues();

		double freqSum = 0;
		for (double f : freq) freqSum+= f;
		if (Math.abs(1.0-freqSum)>1e-10)
			throw new RuntimeException("Error: frequencies must add up to 1 but currently add to " + freqSum + ".");


		ntaxa = tree.getLeafNodeCount();

		int contempCount = 0;
		for (Node node : tree.getExternalNodes())
			if (node.getHeight()==0.)
				contempCount++;

		weightOfNodeSubTree = new double[ntaxa * 2];

		isParallelizedCalculation = isParallelizedCalculationInput.get();
		minimalProportionForParallelization = minimalProportionForParallelizationInput.get();

		if(isParallelizedCalculation) executorBootUp();

	}

	/**
	 *
	 * @param t
	 * @param PG0 initial conditions for p0 (0..n-1) and for ge (n..2n-1)
	 * @param t0
	 * @param PG
	 * @return
	 */
	public p0ge_InitialConditions getG(double t, p0ge_InitialConditions PG0, double t0, p0ge_ODE PG){

		try {

			if (Math.abs(PG.origin -t) < globalPrecisionThreshold|| Math.abs(t0-t) < globalPrecisionThreshold ||  PG.origin < t) {
				return PG0;
			}

			double from = t;
			double to = t0;
			double oneMinusRho;

			int indexFrom = Utils.index(from, PG.times, PG.nIntervals);
			int index = Utils.index(to, PG.times, PG.nIntervals);

			int steps = index - indexFrom;
			if (Math.abs(from-PG.times[indexFrom]) < globalPrecisionThreshold ) steps--;
			if (index>0 && Math.abs(to-PG.times[index-1]) < globalPrecisionThreshold ) {
				steps--;
				index--;
			}
			index--;

			// pgScaled contains the set of initial conditions scaled made to fit the requirements on the values 'double' can represent. It also contains the factor by which the numbers were multiplied
			ScaledNumbers pgScaled = SmallNumberScaler.scale(PG0);

			while (steps > 0){

				from = PG.times[index];

				pgScaled = safeIntegrate(PG, to, pgScaled, from); // solve PG , store solution temporarily integrationResults

				// 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
				PG0 = SmallNumberScaler.unscale(pgScaled.getEquation(), pgScaled.getScalingFactor());


                for (int i=0; i<n; i++){
                    oneMinusRho = 1-PG.rho[index][i];
                    PG0.conditionsOnP[i] *= oneMinusRho;
                    PG0.conditionsOnG[i] = PG0.conditionsOnG[i].scalarMultiply(oneMinusRho);
				}

				to = PG.times[index];

				steps--;
				index--;

				// 'rescale' the results of the last integration to prepare for the next integration step
				pgScaled = SmallNumberScaler.scale(PG0);
			}

			pgScaled = safeIntegrate(PG, to, pgScaled, t); // solve PG , store solution temporarily integrationResults

			// 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
			PG0 = SmallNumberScaler.unscale(pgScaled.getEquation(), pgScaled.getScalingFactor());

		}catch(Exception e){
			// e.printStackTrace(); // for debugging

			throw new RuntimeException("couldn't calculate g");
		}

		return PG0;
	}

	/**
	 * Perform an initial traversal of the tree to get the 'weights' (sum of all its edges lengths) of all sub-trees
	 * Useful for performing parallelized calculations on the tree.
	 * The weights of the subtrees tell us the depth at which parallelization should stop, so as to not parallelize on subtrees that are too small.
	 * Results are stored in 'weightOfNodeSubTree' array
	 * @param tree
	 */
	public void getAllSubTreesWeights(TreeInterface tree){
		Node root = tree.getRoot();
		double weight = 0;
		for(final Node child : root.getChildren()) {
			weight += getSubTreeWeight(child);
		}
		weightOfNodeSubTree[root.getNr()] = weight;
	}

	/**
	 * Perform an initial traversal of the subtree to get its 'weight': sum of all its edges.
	 * @param node
	 * @return
	 */
	public double getSubTreeWeight(Node node){

		// if leaf, stop recursion, get length of branch above and return
		if(node.isLeaf()) {
			weightOfNodeSubTree[node.getNr()] = node.getLength();
			return node.getLength();
		}

		// else, iterate over the children of the node
		double weight = 0;
		for(final Node child : node.getChildren()) {
			weight += getSubTreeWeight(child);
		}
		// add length of parental branch
		weight += node.getLength();
		// store the value
		weightOfNodeSubTree[node.getNr()] = weight;

		return weight;
	}


	void updateParallelizationThreshold(){
		if(isParallelizedCalculation) {
			getAllSubTreesWeights(tree);
			// set 'parallelizationThreshold' to a fraction of the whole tree weight.
			// The size of this fraction is determined by a tuning parameter. This parameter should be adjusted (increased) if more computation cores are available
			parallelizationThreshold = weightOfNodeSubTree[tree.getRoot().getNr()] * minimalProportionForParallelization;
		}
	}


	void setupIntegrators(){   // set up ODE's and integrators

		//TODO set these minstep and maxstep to be a class field
		if (minstep == null) minstep = parameterization.getOrigin()*1e-100;
		if (maxstep == null) maxstep = parameterization.getOrigin()/10;

		P = new p0_ODE(parameterization);
		PG = new p0ge_ODE(parameterization, P, maxEvaluations.get());

		p0ge_ODE.globalPrecisionThreshold = globalPrecisionThreshold;

        pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
        PG.p_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
	}

	/**
	 * Perform the integration of PG with initial conds in pgScaled between to and from
	 * Use an adaptive-step-size integrator
	 * "Safe" because it divides the integration interval in two
	 * if the interval is (arbitrarily) judged to be too big to give reliable results
	 * @param PG
	 * @param to
	 * @param pgScaled
	 * @param from
	 * @return
	 */
	public static ScaledNumbers safeIntegrate(p0ge_ODE PG, double to, ScaledNumbers pgScaled, double from){

		// if the integration interval is too small, nothing is done (to prevent infinite looping)
		if(Math.abs(from-to) < globalPrecisionThreshold /*(T * 1e-20)*/) return pgScaled;

		//TODO make threshold a class field
		if(PG.origin >0 && Math.abs(from-to)>PG.origin /6 ) {
			pgScaled = safeIntegrate(PG, to, pgScaled, from + (to-from)/2);
			pgScaled = safeIntegrate(PG, from + (to-from)/2, pgScaled, from);
		} else {

			//setup of the relativeTolerance and absoluteTolerance input of the adaptive integrator
			//TODO set these two as class fields
			double relativeToleranceConstant = 1e-7;
			double absoluteToleranceConstant = 1e-100;
			double[] absoluteToleranceVector = new double [2*n];
			double[] relativeToleranceVector = new double [2*n];

			for(int i = 0; i<n; i++) {
				absoluteToleranceVector[i] = absoluteToleranceConstant;
				if(pgScaled.getEquation()[i+n] > 0) { // adapt absoluteTolerance to the values stored in pgScaled
					absoluteToleranceVector[i+n] = Math.max(1e-310, pgScaled.getEquation()[i+n]*absoluteToleranceConstant);
				} else {
					absoluteToleranceVector[i+n] = absoluteToleranceConstant;
				}
				relativeToleranceVector[i] = relativeToleranceConstant;
				relativeToleranceVector[i+n] = relativeToleranceConstant;
			}

			double[] integrationResults = new double[pgScaled.getEquation().length];
			int a = pgScaled.getScalingFactor(); // store scaling factor
			int n = pgScaled.getEquation().length/2; // dimension of the ODE system


			FirstOrderIntegrator integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteToleranceVector, relativeToleranceVector);
			integrator.integrate(PG, to, pgScaled.getEquation(), from, integrationResults); // perform the integration step

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
	 * Find all initial conditions for all future integrations on p0 equations
	 * @param tree
	 * @return an array of arrays storing the initial conditions values
	 */
	public double[][] getAllInitialConditionsForP(TreeInterface tree){

		int leafCount = tree.getLeafNodeCount();
		double[] leafHeights = new double[leafCount];
		int[] indicesSortedByLeafHeight  =new int[leafCount];

		for (int i=0; i<leafCount; i++){ // get all leaf heights
			leafHeights[i] = parameterization.getOrigin() - tree.getNode(i).getHeight();
			// System.out.println(nodeHeight[i]);
			indicesSortedByLeafHeight[i] = i;
		}

		HeapSort.sort(leafHeights, indicesSortedByLeafHeight); // sort leafs in order their height in the tree
		//"sort" sorts in ascending order, so we have to be careful since the integration starts from the leaves at height T and goes up to the root at height 0 (or >0)

		double[][] pInitialCondsAtLeaves = new double[leafCount + 1][n];

		double t = leafHeights[indicesSortedByLeafHeight[leafCount-1]];

		pInitialCondsAtLeaves[indicesSortedByLeafHeight[leafCount-1]] = PG.getP(t);
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
					//TODO the integration performed in getP is done before all the other potentially-parallelized getG, so should not matter that it has its own integrator, but if it does (or to simplify the code), take care of passing an integrator as a local variable
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] = PG.getP(t, pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]], t0);
					t0 = t;
				}

			}
		}

		pInitialCondsAtLeaves[leafCount] = PG.getP(0, pInitialCondsAtLeaves[indicesSortedByLeafHeight[0]], t0);

		return pInitialCondsAtLeaves;
	}

	static void executorBootUp(){
		executor = Executors.newCachedThreadPool();
		pool = (ThreadPoolExecutor) executor;
	}

	static void executorShutdown(){
		pool.shutdown();
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

	abstract class TraversalService implements Callable<p0ge_InitialConditions> {

		protected Node rootSubtree;
		protected double from;
		protected double to;
		protected p0ge_ODE PG;
		protected FirstOrderIntegrator pg_integrator;

		public TraversalService(Node root, double from, double to) {
			this.rootSubtree = root;
			this.from = from;
			this.to = to;
			this.setupODEs();
		}

		private void setupODEs(){  // set up ODE's and integrators

			//TODO set minstep and maxstep to be PiecewiseBDDistr fields
			if (minstep == null) minstep = parameterization.getOrigin()*1e-100;
			if (maxstep == null) maxstep = parameterization.getOrigin()/10;

			PG = new p0ge_ODE(parameterization, P, maxEvaluations.get());

			p0ge_ODE.globalPrecisionThreshold = globalPrecisionThreshold;

			pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
		}

		abstract protected p0ge_InitialConditions calculateSubtreeLikelihoodInThread();

		@Override
		public p0ge_InitialConditions call() throws Exception {
			// traverse the tree in a potentially-parallelized way
			return calculateSubtreeLikelihoodInThread();
		}
	}

}
