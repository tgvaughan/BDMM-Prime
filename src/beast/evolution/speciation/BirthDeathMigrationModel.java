package beast.evolution.speciation;

import beast.core.Description;
import beast.core.util.Utils;
import beast.evolution.tree.*;
import beast.core.Input;

import beast.math.*;

import java.util.concurrent.*;


/**
 * @author Denise Kuehnert
 *         Date: May 25, 2012
 *         Time: 11:38:27 AM
 */

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
		"This should only be used when the migration process along the phylogeny is important. Otherwise the computationally less intense BirthDeathMigrationModelUncoloured can be employed.")
public class BirthDeathMigrationModel extends PiecewiseBirthDeathMigrationDistribution {

	// !!! TODO: test birth among deme implementation!!!

	public Input<MultiTypeRootBranch> originBranchInput =
			new Input<>("originBranch", "MultiTypeRootBranch for origin coloring");

	MultiTypeRootBranch originBranch;

	Boolean print = false;

	@Override
	public void initAndValidate() {

		if (birthAmongDemes && migrationMatrix.get()!=null) throw new RuntimeException("Error in BDMM setup: When using MultiTypeTrees there can be migration OR transmission among types, but not  both.");
		super.initAndValidate();
	}

	@Override
	void computeRhoTips(){

		double tipTime;

		for (Node tip : treeInput.get().getExternalNodes()) {

			tipTime = T-tip.getHeight();
			isRhoTip[tip.getNr()] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TODO: DEAL WITH THE IMPLICIT THRESHOLD HERE, THAT SHOULD WORK WITH THE OTHERS (and probably same thing in uncoloured)
				// may need to change to 1e-10
				if (Math.abs(time-tipTime) < 1e-10 && rho[((MultiTypeNode)tip).getNodeType()*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

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
	 * @param isMigrationEvent
	 * @return
	 */
	public static p0ge_InitialConditions getG(double t, p0ge_InitialConditions PG0, double t0, p0ge_ODE PG, Node node, boolean isMigrationEvent){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		if (node.isLeaf() && !isMigrationEvent){ //TODO understand why the !isMigrationEvent here and document it (or remove it) //bc otherwise pb with getP ?
			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, n);
		}

		return getG(t,  PG0,  t0, PG);
	}

	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {

		if (SAModel && treeInput.isDirty()) throw new RuntimeException("Error: SA Model only implemented for fixed trees!");

		MultiTypeNode root = (MultiTypeNode) tree.getRoot();

		if (!((MultiTypeTree) tree).isValid() || (origin.get()!=null && !originBranchIsValid(root, birthAmongDemes))){
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		int node_state;
		if (origin.get()==null) {
			T = root.getHeight();
			node_state =  ((MultiTypeNode) tree.getRoot()).getNodeType();
		}
		else{
			updateOrigin(root);
			node_state = (originBranch.getChangeCount()>0) ?
					originBranch.getChangeType(originBranch.getChangeCount()-1) :
					((MultiTypeNode) tree.getRoot()).getNodeType();
			if (orig < 0)
				return Double.NEGATIVE_INFINITY;
		}

		collectTimes(T);
		setRho();

		if (updateRates() < 0 ||  (times[totalIntervals-1] > T)) {
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		double[] noSampleExistsProp =  new double[n];

		// update the threshold for parallelization
		//TODO only do it if tree shape changed
		updateParallelizationThreshold();

		try{  // start calculation

			pInitialConditions = getAllInitialConditionsForP(tree);

			if (conditionOnSurvival.get()) {

				noSampleExistsProp = pInitialConditions[pInitialConditions.length-1];

				if (print) System.out.println("\nnoSampleExistsProp = " + noSampleExistsProp[0]);// + ", " + noSampleExistsProp[1]);

				if ((noSampleExistsProp[node_state] < 0) || (noSampleExistsProp[node_state] > 1) || (Math.abs(1 - noSampleExistsProp[node_state]) < 1e-14)) {
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}
			}

			p0ge_InitialConditions pSN;

			//TODO remove these executorBootUp and shutdown if keeping the threadpool alive during the whole MCMC works
			//if(isParallelizedCalculation) {executorBootUp();}

			if (orig>0){
				if (originBranch.getChangeCount()>0) {
					pSN = calculateOriginLikelihood(originBranch.getChangeCount()-1, 0, T-originBranch.getChangeTime(originBranch.getChangeCount()-1) );
				} else {
					pSN = calculateSubtreeLikelihood(root, false, null, 0, orig , PG);
				}

			} else {
				int childIndex = 0;
				if (root.getChild(1).getNr() > root.getChild(0).getNr())
					childIndex = 1; // always start with the same child to avoid numerical differences

				double t0 = T - root.getChild(childIndex).getHeight();
				int childChangeCount = ((MultiTypeNode) root.getChild(childIndex)).getChangeCount();
				if (childChangeCount > 0)
					t0 = T - ((MultiTypeNode) root.getChild(childIndex)).getChangeTime(childChangeCount - 1);


				pSN = calculateSubtreeLikelihood(root.getChild(childIndex), false, null, 0., t0, PG);

				childIndex = Math.abs(childIndex - 1);

				t0 = T - root.getChild(childIndex).getHeight();
				childChangeCount = ((MultiTypeNode) root.getChild(childIndex)).getChangeCount(); // changeCounts[root.getChild(1).getNr()];
				if (childChangeCount > 0)
					t0 = T - ((MultiTypeNode) root.getChild(childIndex)).getChangeTime(childChangeCount - 1);

				p0ge_InitialConditions p1SN;

				p1SN = calculateSubtreeLikelihood(root.getChild(childIndex), false, null, 0., t0, PG);

				for (int i=0; i<pSN.conditionsOnG.length; i++) pSN.conditionsOnG[i] = SmallNumber.multiply(pSN.conditionsOnG[i], p1SN.conditionsOnG[i]);

			}
			if (conditionOnSurvival.get()) {
				pSN.conditionsOnG[node_state] = pSN.conditionsOnG[node_state].scalarMultiply(1/(1-noSampleExistsProp[node_state]));    // condition on survival
			}

			logP = Math.log(freq[node_state]) +  pSN.conditionsOnG[node_state].log();

		}catch(Exception e){
			logP =  Double.NEGATIVE_INFINITY;

			//if(isParallelizedCalculation) executorShutdown();

			return logP;
		}

		if (print) System.out.println("final logL = " + logP);

		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		if (SAModel && !(removalProbability.get().getDimension()==n && removalProbability.get().getValue()==1.)) {
			int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
			logP +=  Math.log(2)*internalNodeCount;
		}

		//if (isParallelizedCalculation) executorShutdown();

		return logP;
	}

	/**
	 *
	 * @param migIndex
	 * @param from
	 * @param to
	 * @return
	 */
	p0ge_InitialConditions calculateOriginLikelihood(Integer migIndex, double from, double to) {

		double[] pconditions = new double[n];
		SmallNumber[] gconditions = new SmallNumber[n];
		for (int i=0; i<n; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int index = Utils.index(to, times, totalIntervals);

		int prevcol = originBranch.getChangeType(migIndex);
		int col =  (migIndex > 0)?  originBranch.getChangeType(migIndex-1):  ((MultiTypeNode) tree.getRoot()).getNodeType();

		migIndex--;

		p0ge_InitialConditions g ;

		if (migIndex >= 0){

			g = calculateOriginLikelihood(migIndex, to, T - originBranch.getChangeTime(migIndex));

			System.arraycopy(g.conditionsOnP, 0, pconditions, 0, n);

			if (birthAmongDemes)
				init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(b_ij[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);
			else
				init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);


			return getG(from,  init,  to, PG);

		}
		else {
			g = calculateSubtreeLikelihood(tree.getRoot(), false, null, to, orig, PG);

			System.arraycopy(g.conditionsOnP, 0, pconditions, 0, n);
			if (birthAmongDemes)
				init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(b_ij[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);
			else
				init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);		// with ratechange in M

			return getG(from, init, to, PG, tree.getRoot(), false);
		}
	}

	p0ge_InitialConditions calculateSubtreeLikelihood(Node node, Boolean isMigrationEvent, Integer migrationIndex, double from, double to, p0ge_ODE PG) {

		double[] pconditions = new double[n];
		SmallNumber[] gconditions = new SmallNumber[n];
		for (int i=0; i<n; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int nodestate = ((MultiTypeNode)node).getNodeType();
		int index = Utils.index(to, times, totalIntervals);

		if (isMigrationEvent){ // migration event

			int prevcol = ((MultiTypeNode) node).getChangeType(migrationIndex);
			int col =  (migrationIndex > 0)?  ((MultiTypeNode) node).getChangeType(migrationIndex-1):  ((MultiTypeNode) node).getNodeType();
			double time ;

			migrationIndex--;

			time = (migrationIndex >= 0)? ((MultiTypeNode) node).getChangeTime(migrationIndex) :node.getHeight();
			p0ge_InitialConditions g = calculateSubtreeLikelihood(node, (migrationIndex >= 0), migrationIndex, to, T-time, PG);

			System.arraycopy(g.conditionsOnP, 0, init.conditionsOnP, 0, n);
			if (birthAmongDemes) // this might be a birth among demes where only the child with the different type got sampled
				init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(b_ij[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);
			if (M[0]!=null)     // or it really is a migration event
				init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);

			return getG(from, init, to, PG, node, true);
		}

		else {

			if (migrationIndex==null &&  ((MultiTypeNode)node).getChangeCount()>0){ // node has migration event(psi)
				return calculateSubtreeLikelihood(node, true, ((MultiTypeNode)node).getChangeCount()-1, from, to, PG) ;
			}

			else{

				if (node.isLeaf()){ // sampling event

					if (!isRhoTip[node.getNr()]){

						init.conditionsOnG[nodestate] = SAModel
								? new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]*(1-r[nodestate * totalIntervals + index]))
								*psi[nodestate * totalIntervals + index])

								: new SmallNumber(psi[nodestate * totalIntervals + index]);

					} else {
						init.conditionsOnG[nodestate] = SAModel?
								new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]/(1-rho[nodestate*totalIntervals+index])*(1-r[nodestate * totalIntervals + index]))
										*rho[nodestate*totalIntervals+index])  :
								new SmallNumber(rho[nodestate*totalIntervals+index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ+δ)) //the +δ is translated by dividing p_i with 1-ρ_i (otherwise there's one too many "*ρ_i" )
					}

					if (print) System.out.println("Sampling at time " + to);

					return getG(from, init, to, PG, node, false);
				}

				else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

					if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

						if (r==null)
							throw new RuntimeException("Error: Sampled ancestor found, but removalprobability not specified!");

						int childIndex = 0;

						if (node.getChild(childIndex).isDirectAncestor()) childIndex = 1;

						p0ge_InitialConditions g = calculateSubtreeLikelihood(node.getChild(childIndex), false, null, to, T - node.getChild(childIndex).getHeight(), PG);

						int saNodeState = ((MultiTypeNode) node.getChild(childIndex ^ 1)).getNodeType(); // get state of direct ancestor, XOR operation gives 1 if childIndex is 0 and vice versa

						if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

							init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState];
							init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(psi[saNodeState * totalIntervals + index]
									* (1-r[saNodeState * totalIntervals + index]));

							//							System.out.println("SA but not rho sampled");

						} else {
							// TODO Change: can be dealt with with getAllPInitialConds
							init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState]*(1-rho[saNodeState*totalIntervals+index]) ;
							init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(rho[saNodeState*totalIntervals+index]
									* (1-r[saNodeState * totalIntervals + index]));

						}

					}

					else {   // birth / infection event

						int indexFirstChild = 0;
						if (node.getChild(1).getNr() > node.getChild(0).getNr()) indexFirstChild = 1; // always start with the same child to avoid numerical differences

						int indexSecondChild = Math.abs(indexFirstChild-1);

						double t0 = T - node.getChild(indexFirstChild).getHeight();
						int childChangeCount = ((MultiTypeNode)node.getChild(indexFirstChild)).getChangeCount();
						if (childChangeCount > 0)
							t0 = T - ((MultiTypeNode)node.getChild(indexFirstChild)).getChangeTime(childChangeCount-1);


						double t1 = T - node.getChild(indexSecondChild).getHeight();
						childChangeCount = ((MultiTypeNode)node.getChild(indexSecondChild)).getChangeCount();
						if (childChangeCount > 0)
							t1 = T - ((MultiTypeNode)node.getChild(indexSecondChild)).getChangeTime(childChangeCount-1);

						p0ge_InitialConditions g0 = new p0ge_InitialConditions();
						p0ge_InitialConditions g1 = new p0ge_InitialConditions();

						// if the calculations are parallelized,
						// evaluate if the next step in the traversal should be split between one new thread and the currrent thread and run in parallel,
						// the split is made if the two subtrees of the current node are bigger than a set threshold.
						if(isParallelizedCalculation
								&& weightOfNodeSubTree[node.getChild(indexFirstChild).getNr()] >  parallelizationThreshold
								&& weightOfNodeSubTree[node.getChild(indexSecondChild).getNr()] > parallelizationThreshold){

							try {
								// start a new thread to take care of the second subtree
								Future<p0ge_InitialConditions> secondChildTraversal = pool.submit(
										new TraversalServiceColoured(node.getChild(indexSecondChild), false, null, to, t1));

								g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), false, null, to, t0, PG);
								g1 = secondChildTraversal.get();

							} catch (Exception e) {
								e.printStackTrace();
								//TODO deal with exceptions properly, maybe do the traversal serially if something failed.
							}
						} else {
							g0 = calculateSubtreeLikelihood(node.getChild(indexFirstChild), false, null, to, t0, PG);
							g1 = calculateSubtreeLikelihood(node.getChild(indexSecondChild), false, null, to, t1, PG);
						}

						System.arraycopy(g0.conditionsOnP, 0, init.conditionsOnP, 0, n);

						if (((MultiTypeNode) node.getChild(0)).getFinalType() == nodestate && nodestate == ((MultiTypeNode) node.getChild(1)).getFinalType()) { // within type transmission event

							init.conditionsOnG[nodestate] = SmallNumber.multiply(g0.conditionsOnG[nodestate], g1.conditionsOnG[nodestate]).scalarMultiply(birth[nodestate * totalIntervals + index]);

						} else { // among type transmission event

							if 	(((MultiTypeNode) node.getChild(0)).getFinalType() != nodestate && nodestate != ((MultiTypeNode) node.getChild(1)).getFinalType())
								throw new RuntimeException("Error: Invalid tree (both children have typeChange event at parent node!");

							int child = (((MultiTypeNode) node.getChild(0)).getFinalType() != nodestate) ? 0 : 1;
							int childstate = ((MultiTypeNode)node.getChild(child)).getFinalType();

							init.conditionsOnG[nodestate] =
									SmallNumber.multiply(g0.conditionsOnG[child==0? childstate : nodestate], g1.conditionsOnG[child==1? childstate : nodestate]).scalarMultiply(b_ij[totalIntervals * (childstate * (n - 1) + (nodestate < childstate ? nodestate : nodestate - 1)) + index]);

						}
					}
				}
			}
		}

		return getG(from, init, to, PG, node, false);
	}


	@Override
	void checkOrigin(TreeInterface tree){
		if (origin.get()==null){
			T = tree.getRoot().getHeight();
		}
		else {
			originBranch = originBranchInput.get();

			if (originBranch==null)  throw new RuntimeException("Error: Origin specified but originBranch missing!");

			updateOrigin(tree.getRoot());

			if (!Boolean.valueOf(System.getProperty("beast.resume")) && orig < 0)
				throw new RuntimeException("Error: origin("+T+") must be larger than tree height("+tree.getRoot().getHeight()+")!");

		}
	}

	public Boolean originBranchIsValid(MultiTypeNode root, Boolean allowChangeAtNode){

		int count = originBranch.getChangeCount();

		if (count>0){

			if (originBranch.getChangeTime(0) < root.getHeight() || originBranch.getChangeTime(count-1) > origin.get().getValue() )
				return false;

			if (!allowChangeAtNode && originBranch.getChangeType(0) == root.getFinalType())
				return false;

			for (int i=1; i<count; i++){
				if (originBranch.getChangeType(i-1) == originBranch.getChangeType(i))
					return false;
			}
		}
		return true;
	}

	class TraversalServiceColoured extends TraversalService{

		private Boolean isMigrationEvent;
		private Integer migrationIndex;

		public TraversalServiceColoured(Node root, Boolean isMigrationEvent, Integer migrationIndex, double from, double to) {

			super(root, from, to, true);
			this.isMigrationEvent = isMigrationEvent;
			this.migrationIndex = migrationIndex;
		}

		@Override
		protected p0ge_InitialConditions calculateSubtreeLikelihoodInThread() {

			return calculateSubtreeLikelihood(rootSubtree, isMigrationEvent, migrationIndex, from, to, PG);
		}

	}

}

