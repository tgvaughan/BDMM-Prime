package bdmm.distributions;

import bdmm.parameterization.Parameterization;
import beast.core.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import bdmm.util.Utils;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.TreeInterface;
import beast.util.HeapSort;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.util.List;
import java.util.Random;
import java.util.concurrent.*;

/**
 * @author Denise Kuehnert
 * Date: Jul 2, 2013
 * Time: 10:28:16 AM
 *
 */

@Citation(value="Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. (2016). " +
		"Phylodynamics with migration: \n" +
		"A computational framework to quantify population structure from genomic data. \n" +
		"Mol Biol Evol. 33(8):2102–2116."
		, DOI= "10.1093/molbev/msw064", year = 2016, firstAuthorSurname = "Kuehnert")

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
        "This implementation also works with sampled ancestor trees.")
public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

	public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The frequencies for each type",
            Input.Validate.REQUIRED);

	public Input<TraitSet> tiptypes = new Input<>("tiptypes",
            "trait information for initializing traits " +
                    "(like node types/locations) in the tree");

	public Input<String> typeLabel = new Input<>("typeLabel",
            "type label in tree for initializing traits " +
                    "(like node types/locations) in the tree");

	public Input<IntegerParameter> tipTypeArray = new Input<>("tipTypeArray",
            "integer array of traits (like node types/locations) " +
                    "in the tree, index corresponds to node number in tree");

	public Input<Integer> maxEvaluations = new Input<>("maxEvaluations",
            "The maximum number of evaluations for ODE solver",
            1000000);

	public Input<Boolean> conditionOnSurvival = new Input<>("conditionOnSurvival",
            "condition on at least one survival? Default true.",
            true);

	public Input<Double> relativeTolerance = new Input<>("relTolerance",
            "relative tolerance for numerical integration",
            1e-7);

	public Input<Double> absoluteTolerance = new Input<>("absTolerance",
            "absolute tolerance for numerical integration",
            1e-100 /*Double.MIN_VALUE*/);

	public Input<Boolean> checkRho = new Input<>("checkRho",
            "check if rho is set if multiple tips are " +
                    "given at present (default true)",
            true);

	public Input<Boolean> parallelizeInput = new Input<>(
	        "parallelize",
            "is the calculation parallelized on sibling subtrees " +
                    "or not (default true)",
            true);

	/* If a large number a cores is available (more than 8 or 10) the
	calculation speed can be increased by diminishing the parallelization
	factor. On the contrary, if only 2-4 cores are available, a slightly
	higher value (1/5 to 1/8) can be beneficial to the calculation speed. */
	public Input<Double> minimalProportionForParallelizationInput = new Input<>(
	        "parallelizationFactor",
            "the minimal relative size the two children " +
                    "subtrees of a node must have to start parallel " +
                    "calculations on the children. (default: 1/10). ",
            1.0/10);

	public Input<Boolean> storeNodeTypes = new Input<>("storeNodeTypes",
            "store tip node types? this assumes that tip types cannot " +
                    "change (default false)", false);

	private int[] nodeStates;

//	private final boolean print = false;
    private final boolean print = true;

	private double[] rootTypeProbs, storedRootTypeProbs;
	private boolean[] isRhoTip;

    private Parameterization parameterization;

    private P0System P;
	private P0GeSystem PG;
    private static Double minstep, maxstep;


	private static double[][] pInitialConditions;

	private static boolean isParallelizedCalculation;
	private static double minimalProportionForParallelization;
	private final static double globalPrecisionThreshold = 1e-10;
    private double parallelizationThreshold;
    private static ThreadPoolExecutor pool;

	private double[] weightOfNodeSubTree;

	private TreeInterface tree;

	@Override
	public void initAndValidate() {
		if ((tiptypes.get()==null?0:1) + (typeLabel.get()==null?0:1) + (tipTypeArray.get()==null?0:1) != 1 )
			throw new RuntimeException("Tip types need to be specified exactly once using either tiptypes OR typeLabel OR tipTypeArray.");

	    parameterization = parameterizationInput.get();

		tree = treeInput.get();

		Double factor;

		double freqSum = 0;
		for (double f : frequenciesInput.get().getValues()) freqSum+= f;
		if (Math.abs(1.0-freqSum)>1e-10)
			throw new RuntimeException("Error: frequencies must add up to 1 but currently add to " + freqSum + ".");


        int nLeaves = tree.getLeafNodeCount();

		int contempCount = 0;
		for (Node node : tree.getExternalNodes())
			if (node.getHeight()==0.)
				contempCount++;

		weightOfNodeSubTree = new double[nLeaves * 2];

		isParallelizedCalculation = parallelizeInput.get();
		minimalProportionForParallelization = minimalProportionForParallelizationInput.get();

		if(isParallelizedCalculation) executorBootUp();

		setupIntegrators();

				if (storeNodeTypes.get()) {

			nodeStates = new int[nLeaves];

			for (Node node : tree.getExternalNodes()){
				nodeStates[node.getNr()] = getNodeType(node, true);
			}
		}

		rootTypeProbs = new double[parameterization.getNTypes()];
        storedRootTypeProbs = new double[parameterization.getNTypes()];

        // Determine which, if any, of the leaf ages correspond exactly to
        // rho sampling times.
        isRhoTip = new boolean[tree.getLeafNodeCount()];
        for (int nodeNr = 0; nodeNr < tree.getLeafNodeCount(); nodeNr++) {
            isRhoTip[nodeNr] = false;
            double nodeTime = parameterization.getOrigin()-tree.getNode(nodeNr).getHeight();
            for (double rhoSampTime : parameterization.getRhoSamplingTimes()) {
                if (rhoSampTime == nodeTime) {
                    isRhoTip[nodeNr] = true;
                    break;
                }
            }
        }
	}

	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {

		Node root = tree.getRoot();

		if (parameterization.getOrigin() < tree.getRoot().getHeight()) {
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		// update the threshold for parallelization
		//TODO only do it if tree shape changed
		updateParallelizationThreshold();

        pInitialConditions = getAllInitialConditionsForP(tree);

        double probNoSample = 0;
        if (conditionOnSurvival.get()) {

            double[] noSampleExistsProp = pInitialConditions[pInitialConditions.length-1];
            if (print) System.out.println("\nnoSampleExistsProp = "
                    + noSampleExistsProp[0] + ", " + noSampleExistsProp[1]);

            for (int rootType = 0; rootType<parameterization.getNTypes(); rootType++){
                probNoSample += frequenciesInput.get().getArrayValue(rootType) *  noSampleExistsProp[rootType] ;
            }

            if (probNoSample<0 || probNoSample>1)
                return Double.NEGATIVE_INFINITY;

        }

        P0GeState finalP0Ge = calculateSubtreeLikelihood(root,0,
                parameterization.getOrigin() - tree.getRoot().getHeight(), PG);

        if (print) System.out.print("final p per state = ");

        SmallNumber PrSN = new SmallNumber(0);
        for (int rootType = 0; rootType<parameterization.getNTypes(); rootType++){

            SmallNumber jointProb = finalP0Ge
                    .conditionsOnG[rootType]
                    .scalarMultiply(frequenciesInput.get().getArrayValue(rootType));
            if (jointProb.getMantissa()>0 ) {
                rootTypeProbs[rootType] = jointProb.log();
                PrSN = SmallNumber.add(PrSN, jointProb);
            } else {
                rootTypeProbs[rootType] = Double.NEGATIVE_INFINITY;
            }

            if (print) System.out.print(finalP0Ge.conditionsOnP[rootType] + "\t"
                    + finalP0Ge.conditionsOnG[rootType] + "\t");
        }

        // Normalize root type probs:
        for (int rootType = 0; rootType<parameterization.getNTypes(); rootType++) {
            rootTypeProbs[rootType] -= PrSN.log();
            rootTypeProbs[rootType] = Math.exp(rootTypeProbs[rootType]);
        }

        if (conditionOnSurvival.get()){
            PrSN = PrSN.scalarMultiply(1/(1-probNoSample));
        }

		logP = PrSN.log();

		if (print) System.out.println("\nlogP = " + logP);

		// TGV: Why is this necessary?
		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		return logP;
	}

	private int getNodeType(Node node, Boolean init){

		try {

			if (!storeNodeTypes.get() || init){

				int nodestate = -1;

				if (tiptypes.get() != null) {

                    nodestate = (int) tiptypes.get().getValue((node.getID()));

                } else if (typeLabel.get()!=null) {

                       Object d = node.getMetaData(typeLabel.get());

                       if (d instanceof Integer) nodestate = (Integer) d;
                       else if (d instanceof Double) nodestate = ((Double) d).intValue();
                       else if (d instanceof int[]) nodestate = ((int[]) d)[0];
                       else if (d instanceof String) nodestate = Integer.valueOf((String)d);
                       else
                           throw new RuntimeException("Error interpreting as type index: " + d);

                } else {

                    nodestate = (int) tipTypeArray.get().getArrayValue((node.getNr()));
                }

                if (nodestate < -1)
                    throw new ConstraintViolatedException("State assignment failed.");

				return nodestate;

			}
			else return nodeStates[node.getNr()];

		} catch(Exception e){
			throw new ConstraintViolatedException("Something went wrong with the assignment of types to the nodes (node ID="+node.getID()+"). Please check your XML file!");
		}
	}

	public int getLeafStateForLogging(Node node, long sampleNr) {
		if(!node.isLeaf()) {
			throw new IllegalArgumentException("Node should be a leaf");
		}
		return getNodeType(node, sampleNr==0);
	}


	/**
	 *
	 * @param node Node below edge.
	 * @param tStart Time of start (top) of edge.
	 * @param tEnd Time of end (bottom) of edge.
	 * @param PG Object describing ODEs to integrate.
     *
	 * @return State at top of edge.
	 */
    private P0GeState calculateSubtreeLikelihood(Node node, double tStart, double tEnd, P0GeSystem PG) {

		double[] pconditions = new double[parameterization.getNTypes()];
		SmallNumber[] gconditions = new SmallNumber[parameterization.getNTypes()];
		for (int i = 0; i<parameterization.getNTypes(); i++)
		    gconditions[i] = new SmallNumber();

		P0GeState init = new P0GeState(pconditions, gconditions);

		int intervalIdx = Utils.getIntervalIndex(tEnd, PG.intervalStartTimes);

		if (node.isLeaf()){ // sampling event

			int nodeType = getNodeType(node, false);

			if (nodeType==-1) { //unknown state

				//TODO test if SA model case is properly implemented (not tested!)
				for (int type = 0; type < parameterization.getNTypes(); type++) {

					if (!isRhoTip[node.getNr()]) {
						init.conditionsOnG[type] = new SmallNumber(
						        (PG.r[intervalIdx][type] + pInitialConditions[node.getNr()][type]*(1-PG.r[intervalIdx][type]))
										* PG.s[intervalIdx][type]); // with SA: ψ_i(r + (1 − r)p_i(τ))
					}
					else {
						init.conditionsOnG[type] = new SmallNumber(
						        (PG.r[intervalIdx][type] + pInitialConditions[node.getNr()][type]
                                        / (1 - PG.rho[intervalIdx][type]) * (1 - PG.r[intervalIdx][type]))
										* PG.rho[intervalIdx][type]);
					}
				}
			}
			else {

				if (!isRhoTip[node.getNr()]) {
					init.conditionsOnG[nodeType] = new SmallNumber(
					        (PG.r[intervalIdx][nodeType] + pInitialConditions[node.getNr()][nodeType]
                                    * (1-PG.r[intervalIdx][nodeType]))
                                    * PG.s[intervalIdx][nodeType]); // with SA: ψ_i(r + (1 − r)p_i(τ))
				} else {
					init.conditionsOnG[nodeType] = new SmallNumber(
					        (PG.r[intervalIdx][nodeType] + pInitialConditions[node.getNr()][nodeType]
                                    /(1-PG.rho[intervalIdx][nodeType])*(1-PG.r[intervalIdx][nodeType]))
									*PG.rho[intervalIdx][nodeType]);
				}

			}
			if (print) System.out.println("Sampling at time " + (parameterization.getOrigin()-tEnd));

			return getG(tStart, init, tEnd, PG, node);
		}


		else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

			if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

				int childIndex = 0;

				if (node.getChild(childIndex).isDirectAncestor()) childIndex = 1;

				P0GeState g = calculateSubtreeLikelihood(node.getChild(childIndex), tEnd, parameterization.getOrigin() - node.getChild(childIndex).getHeight(), PG);

				int saNodeType = getNodeType(node.getChild(childIndex ^ 1), false); // get state of direct ancestor, XOR operation gives 1 if childIndex is 0 and vice versa

				//TODO test if properly implemented (not tested!)
				if (saNodeType == -1) { // unknown state
					for (int type = 0; type < parameterization.getNTypes(); type++) {
						if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

							init.conditionsOnP[type] = g.conditionsOnP[type];
							init.conditionsOnG[type] = g.conditionsOnG[type].scalarMultiply(PG.s[intervalIdx][type]
									* (1 - PG.r[intervalIdx][type]));

						} else {
							// TODO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
							init.conditionsOnP[type] = g.conditionsOnP[type] * (1 - PG.rho[intervalIdx][type]);
							init.conditionsOnG[type] = g.conditionsOnG[type].scalarMultiply(PG.rho[intervalIdx][type]
									* (1 - PG.r[intervalIdx][type]));

						}
					}
				}
				else {
					if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

						init.conditionsOnP[saNodeType] = g.conditionsOnP[saNodeType];
						init.conditionsOnG[saNodeType] = g.conditionsOnG[saNodeType]
                                .scalarMultiply(PG.s[intervalIdx][saNodeType]
								* (1 - PG.r[intervalIdx][saNodeType]));

//					System.out.println("SA but not rho sampled");

					} else {
						// TODO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
						init.conditionsOnP[saNodeType] = g.conditionsOnP[saNodeType]
                                * (1 - PG.rho[intervalIdx][saNodeType]);
						init.conditionsOnG[saNodeType] = g.conditionsOnG[saNodeType]
                                .scalarMultiply(PG.rho[intervalIdx][saNodeType]
								* (1 - PG.r[intervalIdx][saNodeType]));

					}
				}
			}

			else {   // birth / infection event

				int indexFirstChild = 0;
				if (node.getChild(1).getNr() > node.getChild(0).getNr())
					indexFirstChild = 1; // always start with the same child to avoid numerical differences

				int indexSecondChild = Math.abs(indexFirstChild-1);

				//TODO refactor with more explicit names
				P0GeState g0 = new P0GeState();
				P0GeState g1 = new P0GeState();

				// evaluate if the next step in the traversal should be split between one new thread and the currrent thread and run in parallel.

				if(isParallelizedCalculation
						&& weightOfNodeSubTree[node.getChild(indexFirstChild).getNr()] >  parallelizationThreshold
                        && weightOfNodeSubTree[node.getChild(indexSecondChild).getNr()] > parallelizationThreshold){

				    try {
                        // start a new thread to take care of the second subtree
                        Future<P0GeState> secondChildTraversal = pool.submit(
                                new TraversalServiceUncoloured(node.getChild(indexSecondChild), tEnd,
                                        parameterization.getOrigin() - node.getChild(indexSecondChild).getHeight()));

                        g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), tEnd,
                                parameterization.getOrigin() - node.getChild(indexFirstChild).getHeight(), PG);
                        g1 = secondChildTraversal.get();
                    } catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();

                        System.exit(1);
                    }

                } else {
					g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), tEnd,
                            parameterization.getOrigin() - node.getChild(indexFirstChild).getHeight(), PG);
					g1 = calculateSubtreeLikelihood(node.getChild(indexSecondChild), tEnd,
                            parameterization.getOrigin() - node.getChild(indexSecondChild).getHeight(), PG);
				}


				if (print)
					System.out.println("Infection at time " + (parameterization.getOrigin() - tEnd));//+ " with p = " + p + "\tg0 = " + g0 + "\tg1 = " + g1);


				for (int childType = 0; childType < parameterization.getNTypes(); childType++) {

					if (print) {
						System.out.println("state " + childType + "\t p0 = " + g0.conditionsOnP[childType] + "\t p1 = " + g1.conditionsOnP[childType]);
						System.out.println("\t\t g0 = " + g0.conditionsOnG[childType] + "\t g1 = " + g1.conditionsOnG[childType]);
					}

					init.conditionsOnP[childType] = g0.conditionsOnP[childType];
					init.conditionsOnG[childType] = SmallNumber
                            .multiply(g0.conditionsOnG[childType], g1.conditionsOnG[childType])
                            .scalarMultiply(PG.b[intervalIdx][childType]);

                    for (int otherChildType = 0; otherChildType < parameterization.getNTypes(); otherChildType++) {
                        if (otherChildType == childType)
                            continue;

                        init.conditionsOnG[childType] = SmallNumber.add(
                                init.conditionsOnG[childType],
                                SmallNumber.add(
                                        SmallNumber.multiply(g0.conditionsOnG[childType], g1.conditionsOnG[otherChildType]),
                                        SmallNumber.multiply(g0.conditionsOnG[otherChildType], g1.conditionsOnG[childType]))
                                        .scalarMultiply(0.5 * PG.b_ij[intervalIdx][childType][otherChildType]));
                    }


					if (Double.isInfinite(init.conditionsOnP[childType])) {
						throw new RuntimeException("infinite likelihood");
					}
				}
			}
		}


		else {// found a single child node

			throw new RuntimeException("Error: Single child-nodes found (although not using sampled ancestors)");
		}

		if (print){
			System.out.print("p after subtree merge = ");
			for (int i = 0; i< parameterization.getNTypes(); i++) System.out.print(init.conditionsOnP[i] + "\t");
			for (int i = 0; i< parameterization.getNTypes(); i++) System.out.print(init.conditionsOnG[i] + "\t");
			System.out.println();
		}

		return getG(tStart, init, tEnd, PG, node);
	}


    /**
     * Used to indicate that the state assignment went wrong.
     */
	protected class ConstraintViolatedException extends RuntimeException {
		public ConstraintViolatedException(String s) {
			super(s);
		}

	}

	class TraversalServiceUncoloured extends TraversalService {

		public TraversalServiceUncoloured(Node root, double from, double to) {
			super(root, from, to);
		}

		@Override
		protected P0GeState calculateSubtreeLikelihoodInThread() {
			return calculateSubtreeLikelihood(rootSubtree, from, to, PG);
		}

	}

    /**
     * @return retrieve current set of root type probabilities.
     */
	public double[] getRootTypeProbs() {
	    return rootTypeProbs;
    }


    /**
	 *
	 * @param t
	 * @param PG0
	 * @param t0
	 * @param PG
	 * @param node
	 * @return
	 */
	public P0GeState getG(double t, P0GeState PG0, double t0, P0GeSystem PG, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		if (node.isLeaf()) {
			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, parameterization.getNTypes());
		}

		return getG(t,  PG0,  t0, PG);
	}

	/**
	 *
	 * @param t
	 * @param PG0 initial conditions for p0 (0..n-1) and for ge (n..2n-1)
	 * @param t0
	 * @param PG
	 * @return
	 */
    public P0GeState getG(double t, P0GeState PG0, double t0, P0GeSystem PG){

        if (Math.abs(PG.origin -t) < globalPrecisionThreshold|| Math.abs(t0-t) < globalPrecisionThreshold ||  PG.origin < t) {
            return PG0;
        }

        double from = t;
        double to = t0;
        double oneMinusRho;

        int indexFrom = Utils.getIntervalIndex(from, PG.intervalStartTimes);
        int index = Utils.getIntervalIndex(to, PG.intervalStartTimes);

        int steps = index - indexFrom;
        if (Math.abs(from-PG.intervalStartTimes[indexFrom]) < globalPrecisionThreshold ) steps--;
        if (index>0 && Math.abs(to-PG.intervalStartTimes[index-1]) < globalPrecisionThreshold ) {
            steps--;
            index--;
        }
        index--;

        // pgScaled contains the set of initial conditions scaled made to fit the requirements on the values 'double' can represent. It also contains the factor by which the numbers were multiplied
        ScaledNumbers pgScaled = SmallNumberScaler.scale(PG0);

        while (steps > 0){

            from = PG.intervalStartTimes[index];

            pgScaled = safeIntegrate(PG, to, pgScaled, from); // solve PG , store solution temporarily integrationResults

            // 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
            PG0 = SmallNumberScaler.unscale(pgScaled.getEquation(), pgScaled.getScalingFactor());


            for (int i = 0; i< parameterization.getNTypes(); i++){
                oneMinusRho = 1-PG.rho[index][i];
                PG0.conditionsOnP[i] *= oneMinusRho;
                PG0.conditionsOnG[i] = PG0.conditionsOnG[i].scalarMultiply(oneMinusRho);
            }

            to = PG.intervalStartTimes[index];

            steps--;
            index--;

            // 'rescale' the results of the last integration to prepare for the next integration step
            pgScaled = SmallNumberScaler.scale(PG0);
        }

        pgScaled = safeIntegrate(PG, to, pgScaled, t); // solve PG , store solution temporarily integrationResults

        // 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
        PG0 = SmallNumberScaler.unscale(pgScaled.getEquation(), pgScaled.getScalingFactor());


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


	private void updateParallelizationThreshold(){
		if(isParallelizedCalculation) {
			getAllSubTreesWeights(tree);
			// set 'parallelizationThreshold' to a fraction of the whole tree weight.
			// The size of this fraction is determined by a tuning parameter. This parameter should be adjusted (increased) if more computation cores are available
			parallelizationThreshold = weightOfNodeSubTree[tree.getRoot().getNr()] * minimalProportionForParallelization;
		}
	}


	private void setupIntegrators(){   // set up ODE's and integrators

		//TODO set these minstep and maxstep to be a class field
		if (minstep == null) minstep = parameterization.getOrigin()*1e-100;
		if (maxstep == null) maxstep = parameterization.getOrigin()/10;

		P = new P0System(parameterization);
		PG = new P0GeSystem(parameterization, P, maxEvaluations.get());

		P0GeSystem.globalPrecisionThreshold = globalPrecisionThreshold;

        FirstOrderIntegrator pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
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
	 * @return result of integration
	 */
	private static ScaledNumbers safeIntegrate(P0GeSystem PG, double to, ScaledNumbers pgScaled, double from){

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
			double[] absoluteToleranceVector = new double [2* PG.nTypes];
			double[] relativeToleranceVector = new double [2* PG.nTypes];

			for(int i = 0; i< PG.nTypes; i++) {
				absoluteToleranceVector[i] = absoluteToleranceConstant;
				if(pgScaled.getEquation()[i+ PG.nTypes] > 0) { // adapt absoluteTolerance to the values stored in pgScaled
					absoluteToleranceVector[i+ PG.nTypes] = Math.max(1e-310, pgScaled.getEquation()[i+PG.nTypes]*absoluteToleranceConstant);
				} else {
					absoluteToleranceVector[i+ PG.nTypes] = absoluteToleranceConstant;
				}
				relativeToleranceVector[i] = relativeToleranceConstant;
				relativeToleranceVector[i+ PG.nTypes] = relativeToleranceConstant;
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
			pgScaled = SmallNumberScaler.scale(new P0GeState(pConditions, geConditions));
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
		double[] leafTimes = new double[leafCount];
		int[] indicesSortedByLeafHeight  =new int[leafCount];

		for (int i=0; i<leafCount; i++){ // get all leaf heights
			leafTimes[i] = parameterization.getOrigin() - tree.getNode(i).getHeight();
			// System.out.println(nodeHeight[i]);
			indicesSortedByLeafHeight[i] = i;
		}

		HeapSort.sort(leafTimes, indicesSortedByLeafHeight); // sort leafs in order their time in the tree
		//"sort" sorts in ascending order, so we have to be careful since the integration starts from the leaves at time T and goes up to the root at time 0 (or >0)

		double[][] pInitialCondsAtLeaves = new double[leafCount + 1][PG.nTypes];

		double t = leafTimes[indicesSortedByLeafHeight[leafCount-1]];

		pInitialCondsAtLeaves[indicesSortedByLeafHeight[leafCount-1]] = PG.getP(t);
		double t0 = t;

		if (leafCount >1 ){
			for (int i = leafCount-2; i>-1; i--){
				t = leafTimes[indicesSortedByLeafHeight[i]];

				//If the next higher leaf is actually at the same height, store previous results and skip iteration
				if (Math.abs(t-t0) < globalPrecisionThreshold) {
					t0 = t;
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] = pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]];
					continue;
				} else {
					/* TODO the integration performed in getP is done before all
                       the other potentially-parallelized getG, so should not
                       matter that it has its own integrator, but if it does
                       (or to simplify the code), take care of passing an integrator
                       as a local variable. */
					pInitialCondsAtLeaves[indicesSortedByLeafHeight[i]] =
                            PG.getP(t, pInitialCondsAtLeaves[indicesSortedByLeafHeight[i+1]], t0);
					t0 = t;
				}

			}
		}

		pInitialCondsAtLeaves[leafCount] = PG.getP(0, pInitialCondsAtLeaves[indicesSortedByLeafHeight[0]], t0);

		return pInitialCondsAtLeaves;
	}

	static void executorBootUp(){
        ExecutorService executor = Executors.newCachedThreadPool();
		pool = (ThreadPoolExecutor) executor;
	}

	static void executorShutdown(){
		pool.shutdown();
	}

	abstract class TraversalService implements Callable<P0GeState> {

		protected Node rootSubtree;
		protected double from;
		protected double to;
		protected P0GeSystem PG;
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

			PG = new P0GeSystem(parameterization, P, maxEvaluations.get());

			P0GeSystem.globalPrecisionThreshold = globalPrecisionThreshold;

			pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
		}

		abstract protected P0GeState calculateSubtreeLikelihoodInThread();

		@Override
		public P0GeState call() throws Exception {
			// traverse the tree in a potentially-parallelized way
			return calculateSubtreeLikelihoodInThread();
		}
	}

    /* StateNode implementation */

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

    @Override
    public void store() {
        super.store();

        for (int i = 0; i< parameterization.getNTypes(); i++)
            storedRootTypeProbs[i] = rootTypeProbs[i];
    }

    @Override
    public void restore() {
        super.restore();

        double[] tmp = storedRootTypeProbs;
        rootTypeProbs = storedRootTypeProbs;
        storedRootTypeProbs = tmp;
    }


}
