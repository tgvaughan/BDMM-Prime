package bdmm.distributions;

import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.*;
import beast.core.Input;
import beast.core.Description;
import bdmm.util.Utils;

import java.util.concurrent.*;

/**
 * @author Denise Kuehnert
 * Date: Jul 2, 2013
 * Time: 10:28:16 AM
 *
 */

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
		"This implementation also works with sampled ancestor trees.")
public class BirthDeathMigrationDistribution extends PiecewiseBirthDeathMigrationDistribution implements Loggable {


	public Input<TraitSet> tiptypes = new Input<>("tiptypes", "trait information for initializing traits (like node types/locations) in the tree");
	public Input<String> typeLabel = new Input<>("typeLabel", "type label in tree for initializing traits (like node types/locations) in the tree");
	public Input<IntegerParameter> tipTypeArray = new Input<IntegerParameter>("tipTypeArray", "integer array of traits (like node types/locations) in the tree, index corresponds to node number in tree");

	public Input<Boolean> storeNodeTypes = new Input<>("storeNodeTypes", "store tip node types? this assumes that tip types cannot change (default false)", false);

	private int[] nodeStates;

	Boolean print = false;

	double[] rootTypeProbs, storedRootTypeProbs;

	boolean[] isRhoTip;

	@Override
	public void initAndValidate() {

		if ((tiptypes.get()==null?0:1) + (typeLabel.get()==null?0:1) + (tipTypeArray.get()==null?0:1) != 1 )
			throw new RuntimeException("Tip types need to be specified exactly once using either tiptypes OR typeLabel OR tipTypeArray.");


		super.initAndValidate();

		if (storeNodeTypes.get()) {

			nodeStates = new int[ntaxa];

			for (Node node : tree.getExternalNodes()){
				nodeStates[node.getNr()] = getNodeType(node, true);
			}
		}

		rootTypeProbs = new double[nStates];
        storedRootTypeProbs = new double[nStates];

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

	/**
	 *
	 * @param t
	 * @param PG0
	 * @param t0
	 * @param PG
	 * @param node
	 * @return
	 */
	public p0ge_InitialConditions getG(double t, p0ge_InitialConditions PG0, double t0, p0ge_ODE PG, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		if (node.isLeaf()) {
			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, nStates);
		}

		return getG(t,  PG0,  t0, PG);
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

		double[] noSampleExistsProp ;

		SmallNumber PrSN = new SmallNumber(0);
		double nosample = 0;

		try{  // start calculation

			pInitialConditions = getAllInitialConditionsForP(tree);

			if (conditionOnSurvival.get()) {

				noSampleExistsProp = pInitialConditions[pInitialConditions.length-1];

				if (print) System.out.println("\nnoSampleExistsProp = " + noSampleExistsProp[0] + ", " + noSampleExistsProp[1]);

				for (int root_state = 0; root_state< nStates; root_state++){
					nosample += freq[root_state] *  noSampleExistsProp[root_state] ;
				}

				if (nosample<0 || nosample>1)
					return Double.NEGATIVE_INFINITY;

			}

			p0ge_InitialConditions pSN;

			//if(isParallelizedCalculation) {executorBootUp();}

            pSN = calculateSubtreeLikelihood(root,0, parameterization.getOrigin() - tree.getRoot().getHeight(), PG);

			if (print) System.out.print("final p per state = ");

			for (int root_state = 0; root_state< nStates; root_state++){

                SmallNumber jointProb = pSN.conditionsOnG[root_state].scalarMultiply(freq[root_state]);
				if (jointProb.getMantissa()>0 ) {
				    rootTypeProbs[root_state] = jointProb.log();
                    PrSN = SmallNumber.add(PrSN, jointProb);
                } else {
                    rootTypeProbs[root_state] = Double.NEGATIVE_INFINITY;
                }

				if (print) System.out.print(pSN.conditionsOnP[root_state] + "\t" + pSN.conditionsOnG[root_state] + "\t");
			}

			// Normalize root type probs:
            for (int root_state = 0; root_state< nStates; root_state++) {
                rootTypeProbs[root_state] -= PrSN.log();
                rootTypeProbs[root_state] = Math.exp(rootTypeProbs[root_state]);
            }

			if (conditionOnSurvival.get()){
				PrSN = PrSN.scalarMultiply(1/(1-nosample));
			}

		}catch(Exception e){

			if (e instanceof ConstraintViolatedException){throw e;}

			logP =  Double.NEGATIVE_INFINITY;

			//if(isParallelizedCalculation) executorShutdown();
			return logP;
		}

		logP = PrSN.log();

		if (print) System.out.println("\nlogP = " + logP);

		// TGV: Why is this necessary?
		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		//if(isParallelizedCalculation) executorShutdown();
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
	p0ge_InitialConditions calculateSubtreeLikelihood(Node node, double tStart, double tEnd, p0ge_ODE PG) {

		double[] pconditions = new double[nStates];
		SmallNumber[] gconditions = new SmallNumber[nStates];
		for (int i = 0; i< nStates; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int intervalIdx = Utils.getIntervalIndex(tEnd, PG.intervalStartTimes);

		if (node.isLeaf()){ // sampling event

			int nodeType = getNodeType(node, false);

			//TODO potentially refactor to make few lines below more concise and clearer
			if (nodeType==-1) { //unknown state

				//TODO test if SA model case is properly implemented (not tested!)
				for (int type = 0; type< nStates; type++) {

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

				p0ge_InitialConditions g = calculateSubtreeLikelihood(node.getChild(childIndex), tEnd, parameterization.getOrigin() - node.getChild(childIndex).getHeight(), PG);

				int saNodeType = getNodeType(node.getChild(childIndex ^ 1), false); // get state of direct ancestor, XOR operation gives 1 if childIndex is 0 and vice versa

				//TODO test if properly implemented (not tested!)
				if (saNodeType == -1) { // unknown state
					for (int type = 0; type < nStates; type++) {
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
				p0ge_InitialConditions g0 = new p0ge_InitialConditions();
				p0ge_InitialConditions g1 = new p0ge_InitialConditions();

				// evaluate if the next step in the traversal should be split between one new thread and the currrent thread and run in parallel.

				if(isParallelizedCalculation
						&& weightOfNodeSubTree[node.getChild(indexFirstChild).getNr()] >  parallelizationThreshold
						&& weightOfNodeSubTree[node.getChild(indexSecondChild).getNr()] > parallelizationThreshold){

					try {
						// start a new thread to take care of the second subtree
						Future<p0ge_InitialConditions> secondChildTraversal = pool.submit(
								new TraversalServiceUncoloured(node.getChild(indexSecondChild), tEnd,
                                        parameterization.getOrigin() - node.getChild(indexSecondChild).getHeight()));

						g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), tEnd,
                                parameterization.getOrigin() - node.getChild(indexFirstChild).getHeight(), PG);
						g1 = secondChildTraversal.get();

					} catch (Exception e) {
						e.printStackTrace();
						//TODO deal with exceptions properly, maybe do the traversal serially if something failed.
					}
				} else {
					g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), tEnd,
                            parameterization.getOrigin() - node.getChild(indexFirstChild).getHeight(), PG);
					g1 = calculateSubtreeLikelihood(node.getChild(indexSecondChild), tEnd,
                            parameterization.getOrigin() - node.getChild(indexSecondChild).getHeight(), PG);
				}


				if (print)
					System.out.println("Infection at time " + (parameterization.getOrigin() - tEnd));//+ " with p = " + p + "\tg0 = " + g0 + "\tg1 = " + g1);


				for (int childType = 0; childType < nStates; childType++) {

					if (print) {
						System.out.println("state " + childType + "\t p0 = " + g0.conditionsOnP[childType] + "\t p1 = " + g1.conditionsOnP[childType]);
						System.out.println("\t\t g0 = " + g0.conditionsOnG[childType] + "\t g1 = " + g1.conditionsOnG[childType]);
					}

					init.conditionsOnP[childType] = g0.conditionsOnP[childType];
					init.conditionsOnG[childType] = SmallNumber
                            .multiply(g0.conditionsOnG[childType], g1.conditionsOnG[childType])
                            .scalarMultiply(PG.b[intervalIdx][childType]);

                    for (int otherChildType = 0; otherChildType < nStates; otherChildType++) {
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
			for (int i = 0; i< nStates; i++) System.out.print(init.conditionsOnP[i] + "\t");
			for (int i = 0; i< nStates; i++) System.out.print(init.conditionsOnG[i] + "\t");
			System.out.println();
		}

		return getG(tStart, init, tEnd, PG, node);
	}


	// used to indicate that the state assignment went wrong
	protected class ConstraintViolatedException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public ConstraintViolatedException(String s) {
			super(s);
		}

	}

	class TraversalServiceUncoloured extends TraversalService {

		public TraversalServiceUncoloured(Node root, double from, double to) {
			super(root, from, to);
		}

		@Override
		protected p0ge_InitialConditions calculateSubtreeLikelihoodInThread() {
			return calculateSubtreeLikelihood(rootSubtree, from, to, PG);
		}

	}

    /**
     * @return retrieve current set of root type probabilities.
     */
	public double[] getRootTypeProbs() {
	    return rootTypeProbs;
    }

	/* StateNode implementation */

    @Override
    public void store() {
        super.store();

        for (int i = 0; i< nStates; i++)
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

