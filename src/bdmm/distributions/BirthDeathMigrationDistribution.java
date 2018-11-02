package bdmm.distributions;

import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.*;
import beast.core.Input;
import beast.core.Description;
import bdmm.util.Utils;

import java.io.PrintStream;
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

	@Override
	public void initAndValidate() {

		if ((tiptypes.get()==null?0:1) + (typeLabel.get()==null?0:1) + (tipTypeArray.get()==null?0:1) != 1 )
			throw new RuntimeException("Tip types need to be specified exactly once using either tiptypes OR typeLabel OR tipTypeArray.");


		super.initAndValidate();

		if (storeNodeTypes.get()) {

			nodeStates = new int[ntaxa];

			for (Node node : tree.getExternalNodes()){
				nodeStates[node.getNr()] = getNodeState(node, true);
			}
		}

		rootTypeProbs = new double[n];
        storedRootTypeProbs = new double[n];
	}

	void computeRhoTips(){

		double tipTime;

		for (Node tip : treeInput.get().getExternalNodes()) {

			tipTime = T-tip.getHeight();
			isRhoTip[tip.getNr()] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TODO: make a warning that rho sampling precision is with 1e-10. Maybe do a threshold to the type of dating associated with the data?
				if (Math.abs(time-tipTime) <  globalPrecisionThreshold && rho[getNodeState(tip,false)*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

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
			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, n);
		}

		return getG(t,  PG0,  t0, PG);
	}

	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {

		Node root = tree.getRoot();

		if (origin.get()==null)
			T = root.getHeight();
		else
			updateOrigin(root);


		collectTimes(T);
		setRho();

		if ((orig < 0) || updateRates() < 0 ||  (times[totalIntervals-1] > T)) {
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

				for (int root_state=0; root_state<n; root_state++){
					nosample += freq[root_state] *  noSampleExistsProp[root_state] ;
				}

				if (nosample<0 || nosample>1)
					return Double.NEGATIVE_INFINITY;

			}

			p0ge_InitialConditions pSN;

			//if(isParallelizedCalculation) {executorBootUp();}

			if ( orig > 0 ) {
				pSN = calculateSubtreeLikelihood(root,0,orig, PG);}
			else {

				int childIndex = 0;
				if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

				pSN = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight(), PG);
				childIndex = Math.abs(childIndex-1);

				p0ge_InitialConditions p1SN;
				p1SN = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight(), PG);

				for (int i =0; i<pSN.conditionsOnG.length; i++) pSN.conditionsOnG[i] = SmallNumber.multiply(pSN.conditionsOnG[i], p1SN.conditionsOnG[i]);

			}

			if (print) System.out.print("final p per state = ");

			for (int root_state=0; root_state<n; root_state++){

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
            for (int root_state=0; root_state<n; root_state++) {
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

		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		if (SAModel && !(removalProbability.get().getDimension()==n && removalProbability.get().getValue()==1.)) {
			int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
			logP +=  Math.log(2)*internalNodeCount;
		}

		//if(isParallelizedCalculation) executorShutdown();
		return logP;
	}

	private int getNodeState(Node node, Boolean init){

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

		}catch(Exception e){
			throw new ConstraintViolatedException("Something went wrong with the assignment of types to the nodes (node ID="+node.getID()+"). Please check your XML file!");
		}
	}

	public int getLeafStateForLogging(Node node, long sampleNr) {
		if(!node.isLeaf()) {
			throw new IllegalArgumentException("Node should be a leaf");
		}
		return getNodeState(node, sampleNr==0);
	}


	/**
	 *
	 * @param node
	 * @param from
	 * @param to
	 * @param PG
	 * @return
	 */
	p0ge_InitialConditions calculateSubtreeLikelihood(Node node, double from, double to, p0ge_ODE PG) {

		double[] pconditions = new double[n];
		SmallNumber[] gconditions = new SmallNumber[n];
		for (int i=0; i<n; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int index = Utils.index(to,times, totalIntervals);

		if (node.isLeaf()){ // sampling event

			int nodestate = getNodeState(node, false);

			//TODO potentially refactor to make few lines below more concise and clearer
			if (nodestate==-1) { //unknown state

				//TODO remove entirely when tested
//				if (SAModel)
//					throw new ConstraintViolatedException("SA model not implemented with unknown states!");

				//TODO test if SA model case is properly implemented (not tested!)
				for (int i=0; i<n; i++) {

					if (!isRhoTip[node.getNr()]) {
						init.conditionsOnG[i] = SAModel?
								new SmallNumber((r[i * totalIntervals + index] + pInitialConditions[node.getNr()][i]*(1-r[i * totalIntervals + index]))
										*psi[i * totalIntervals + index]) // with SA: ψ_i(r + (1 − r)p_i(τ))
								: new SmallNumber(psi[i * totalIntervals + index]);
					}
					else {
						init.conditionsOnG[i] = SAModel ?
								new SmallNumber((r[i * totalIntervals + index] + pInitialConditions[node.getNr()][i] / (1 - rho[i * totalIntervals + index]) * (1 - r[i * totalIntervals + index]))
										* rho[i * totalIntervals + index]) :
								new SmallNumber(rho[i * totalIntervals + index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ+δ)) //the +δ is translated by dividing p_i with 1-ρ_i (otherwise there's one too many "*ρ_i" )
					}
				}
			}
			else {

				if (!isRhoTip[node.getNr()]) {

					init.conditionsOnG[nodestate] = SAModel?
							new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]*(1-r[nodestate * totalIntervals + index]))
									*psi[nodestate * totalIntervals + index]) // with SA: ψ_i(r + (1 − r)p_i(τ))
							: new SmallNumber(psi[nodestate * totalIntervals + index]);

				}	else {
					init.conditionsOnG[nodestate] = SAModel?
							new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]/(1-rho[nodestate*totalIntervals+index])*(1-r[nodestate * totalIntervals + index]))
									*rho[nodestate*totalIntervals+index])  :
							new SmallNumber(rho[nodestate*totalIntervals+index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ+δ)) //the +δ is translated by dividing p_i with 1-ρ_i (otherwise there's one too many "*ρ_i" )
				}

			}
			if (print) System.out.println("Sampling at time " + (T-to));

			return getG(from, init, to, PG, node);
		}


		else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

			if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

				if (r==null)
					throw new ConstraintViolatedException("Error: Sampled ancestor found, but removalprobability not specified!");

				int childIndex = 0;

				if (node.getChild(childIndex).isDirectAncestor()) childIndex = 1;

				p0ge_InitialConditions g = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight(), PG);

				int saNodeState = getNodeState(node.getChild(childIndex ^ 1), false); // get state of direct ancestor, XOR operation gives 1 if childIndex is 0 and vice versa

				//TODO test if properly implemented (not tested!)
				if (saNodeState == -1) { // unknown state
					for (int i = 0; i < n; i++) {
						if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

							init.conditionsOnP[i] = g.conditionsOnP[i];
							init.conditionsOnG[i] = g.conditionsOnG[i].scalarMultiply(psi[i * totalIntervals + index]
									* (1 - r[i * totalIntervals + index]));

						} else {
							// TODO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
							init.conditionsOnP[i] = g.conditionsOnP[i] * (1 - rho[i * totalIntervals + index]);
							init.conditionsOnG[i] = g.conditionsOnG[i].scalarMultiply(rho[i * totalIntervals + index]
									* (1 - r[i * totalIntervals + index]));

						}
					}
				}
				else {
					if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

						init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState];
						init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(psi[saNodeState * totalIntervals + index]
								* (1 - r[saNodeState * totalIntervals + index]));

//					System.out.println("SA but not rho sampled");

					} else {
						// TODO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
						init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState] * (1 - rho[saNodeState * totalIntervals + index]);
						init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(rho[saNodeState * totalIntervals + index]
								* (1 - r[saNodeState * totalIntervals + index]));

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
								new TraversalServiceUncoloured(node.getChild(indexSecondChild), to, T - node.getChild(indexSecondChild).getHeight()));

						g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), to, T - node.getChild(indexFirstChild).getHeight(), PG);
						g1 = secondChildTraversal.get();

					} catch (Exception e) {
						e.printStackTrace();
						//TODO deal with exceptions properly, maybe do the traversal serially if something failed.
					}
				} else {
					g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), to, T - node.getChild(indexFirstChild).getHeight(), PG);
					g1 = calculateSubtreeLikelihood(node.getChild(indexSecondChild), to, T - node.getChild(indexSecondChild).getHeight(), PG);
				}


				if (print)
					System.out.println("Infection at time " + (T - to));//+ " with p = " + p + "\tg0 = " + g0 + "\tg1 = " + g1);


				for (int childstate = 0; childstate < n; childstate++) {

					if (print) {
						System.out.println("state " + childstate + "\t p0 = " + g0.conditionsOnP[childstate] + "\t p1 = " + g1.conditionsOnP[childstate]);
						System.out.println("\t\t g0 = " + g0.conditionsOnG[childstate] + "\t g1 = " + g1.conditionsOnG[childstate]);
					}

					init.conditionsOnP[childstate] = g0.conditionsOnP[childstate];
					init.conditionsOnG[childstate] = SmallNumber.multiply(g0.conditionsOnG[childstate], g1.conditionsOnG[childstate]).scalarMultiply(birth[childstate * totalIntervals + index]);

					if (birthAmongDemes) {
						for (int j = 0; j < n; j++) {
							if (childstate != j) {
								init.conditionsOnG[childstate] = SmallNumber.add(init.conditionsOnG[childstate], SmallNumber.add(SmallNumber.multiply(g0.conditionsOnG[childstate], g1.conditionsOnG[j]) , SmallNumber.multiply(g0.conditionsOnG[j], g1.conditionsOnG[childstate]))
										.scalarMultiply(0.5 * b_ij[totalIntervals * (childstate * (n - 1) + (j < childstate ? j : j - 1)) + index]));
							}
						}

					}

					if (Double.isInfinite(init.conditionsOnP[childstate])) {
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
			for (int i=0;i<n;i++) System.out.print(init.conditionsOnP[i] + "\t");
			for (int i=0;i<n;i++) System.out.print(init.conditionsOnG[i] + "\t");
			System.out.println();
		}

		return getG(from, init, to, PG, node);
	}


	// used to indicate that the state assignment went wrong
	protected class ConstraintViolatedException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public ConstraintViolatedException(String s) {
			super(s);
		}

	}

	@Override
	public void init(PrintStream out){

		super.init(out);

		if (tipTypeArray.get()!=null) {
			IntegerParameter types = tipTypeArray.get();
			for (int i = 0; i < types.getDimension(); i++) {
				out.print(tipTypeArray.get().getID() + ("_node") + (i+1) + "\t");
			}

		}

	}

	@Override
	public void log(int sampleNr, PrintStream out) {

		super.log(sampleNr, out);

		if (tipTypeArray.get()!=null) {
			for (int i = 0; i < tipTypeArray.get().getDimension(); i++) {
				out.print(treeInput.get().getNode(i).getID() + "\t");

			}
		}
	}

	class TraversalServiceUncoloured extends TraversalService {

		public TraversalServiceUncoloured(Node root, double from, double to) {
			super(root, from, to, false);
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

        for (int i=0; i<n; i++)
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

