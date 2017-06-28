package beast.evolution.speciation;

import beast.evolution.tree.*;
import beast.core.Input;
import beast.core.Description;
import beast.core.util.Utils;

import beast.math.SmallNumber;
import beast.math.p0ge_InitialConditions;
import beast.util.HeapSort;

// currently cleaning
// removing the overflowing way (done)
// refactor to rename methods appropriately (done)


/**
 * @author Denise Kuehnert
 * Date: Jul 2, 2013
 * Time: 10:28:16 AM
 *
 */

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
		"This should be used when the migration process along the phylogeny is irrelevant. Otherwise the BirthDeathMigrationModel can be employed." +
		"This implementation also works with sampled ancestor trees." +
		"Two implementations are available. The first is the fast classic one; the second one prevents underflowing, using so-called 'SmallNumbers', with the cost of additional computational complexity")
public class BirthDeathMigrationModelUncoloured extends PiecewiseBirthDeathMigrationDistribution {

	public Input<TraitSet> tiptypes = new Input<>("tiptypes", "trait information for initializing traits (like node types/locations) in the tree",  Input.Validate.REQUIRED);
	public Input<String> typeLabel = new Input<>("typeLabel", "type label in tree for initializing traits (like node types/locations) in the tree",  Input.Validate.XOR, tiptypes);

	public Input<Boolean> storeNodeTypes = new Input<>("storeNodeTypes", "store tip node types? this assumes that tip types cannot change (default false)", false);

	private int[] nodeStates;

	Boolean print = false;

	@Override
	public void initAndValidate() {

		super.initAndValidate();

		TreeInterface tree = treeInput.get();

		checkOrigin(tree);

		ntaxa = tree.getLeafNodeCount();

		birthAmongDemes = (birthRateAmongDemes.get() !=null || R0AmongDemes.get()!=null);

		if (storeNodeTypes.get()) {

			nodeStates = new int[ntaxa];

			for (Node node : tree.getExternalNodes()){
				nodeStates[node.getNr()] = getNodeState(node, true);
			}
		}

		int contempCount = 0;
		for (Node node : tree.getExternalNodes())
			if (node.getHeight()==0.)
				contempCount++;


		if (checkRho.get() && contempCount>1 && rho==null)
			throw new RuntimeException("Error: multiple tips given at present, but sampling probability \'rho\' is not specified.");

		collectTimes(T);
		setRho();
	}

	protected Double updateRates(TreeInterface tree) {

		birth = new Double[n*totalIntervals];
		death = new Double[n*totalIntervals];
		psi = new Double[n*totalIntervals];
		b_ij = new Double[totalIntervals*(n*(n-1))];
		M = new Double[totalIntervals*(n*(n-1))];
		if (SAModel) r =  new Double[n * totalIntervals];

		if (transform) {
			transformParameters();
		}
		else {

			Double[] birthAmongDemesRates = new Double[1];

			if (birthAmongDemes) birthAmongDemesRates = birthRateAmongDemes.get().getValues();

			updateBirthDeathPsiParams();

			if (birthAmongDemes) {

				updateAmongParameter(b_ij, birthAmongDemesRates, b_ij_Changes, b_ijChangeTimes);
			}
		}

		Double[] migRates = migrationMatrix.get().getValues();

		updateAmongParameter(M, migRates, migChanges, migChangeTimes);

		updateRho();

		freq = frequencies.get().getValues();

		setupIntegrators();

		return 0.;
	}

	void computeRhoTips(){

		double tipTime;

		for (Node tip : treeInput.get().getExternalNodes()) {

			tipTime = T-tip.getHeight();
			isRhoTip[tip.getNr()] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TO DO: make a warning that rho sampling precision is with 1e-10. Maybe do a threshold to the type of dating associated with the data?
				if (Math.abs(time-tipTime) <  globalPrecisionThreshold && rho[getNodeState(tip,false)*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

			}
		}
	}

	/**
	 * Implementation of getG with Small Number structure for ge equations. Avoids underflowing of integration results.
	 * WARNING: getG and getGSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param t
	 * @param PG0
	 * @param t0
	 * @param node
	 * @return
	 */
	public p0ge_InitialConditions getG(double t, p0ge_InitialConditions PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)


		if (node.isLeaf()) {

			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, n);
			//TO DO clean up
			//System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0.conditionsOnP, 0, n);
		}

		return getG(t,  PG0,  t0, pg_integrator, PG, T, maxEvalsUsed);

	}

	/**
	 * WARNING: calculateTreeLogLikelihood allows use of both classic and non-underflowing methods. Some chunks of code are therefore present in two similar versions in this method.
	 * When modifying one of the versions, one should check if the other version also needs the corresponding changes.
	 */
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {

		Node root = tree.getRoot();

		if (origin.get()==null)
			T = root.getHeight();
		else
			updateOrigin(root);


		collectTimes(T);
		setRho();

		if ((orig < 0) || updateRates(tree) < 0 ||  (times[totalIntervals-1] > T)) {
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

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

			p0ge_InitialConditions pSN = new p0ge_InitialConditions();

			if ( orig > 0 ) {
					pSN = calculateSubtreeLikelihood(root,0,orig);
			} else {
				int childIndex = 0;
				if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

					pSN = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					childIndex = Math.abs(childIndex-1);

					p0ge_InitialConditions p1SN = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					for (int i =0; i<pSN.conditionsOnG.length; i++) pSN.conditionsOnG[i] = SmallNumber.multiply(pSN.conditionsOnG[i], p1SN.conditionsOnG[i]);
			}

			if (print) System.out.print("final p per state = ");

				for (int root_state=0; root_state<n; root_state++){

					if (pSN.conditionsOnG[root_state].getMantissa()>0 ) 
						PrSN = SmallNumber.add(PrSN, pSN.conditionsOnG[root_state].scalarMultiply(freq[root_state]));

					if (print) System.out.print(pSN.conditionsOnP[root_state] + "\t" + pSN.conditionsOnG[root_state] + "\t");
				}

				if (conditionOnSurvival.get()){
					PrSN = PrSN.scalarMultiply(1/(1-nosample));
				}




		}catch(Exception e){

			if (e instanceof ConstraintViolatedException){throw e;}

			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		maxEvalsUsed = Math.max(maxEvalsUsed, PG.maxEvalsUsed);

		logP = PrSN.log();

		if (print) System.out.println("\nlogP = " + logP);

		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		if (SAModel && !(removalProbability.get().getDimension()==n && removalProbability.get().getValue()==1.)) {
			int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
			logP +=  Math.log(2)*internalNodeCount;
		}
		return logP;
	}

	private int getNodeState(Node node, Boolean init){

		try {

			if (!storeNodeTypes.get() || init){

				int nodestate = tiptypes.get() != null ?
						(int) tiptypes.get().getValue((node.getID())) :
							((node instanceof MultiTypeNode) ? ((MultiTypeNode) node).getNodeType() : -2);

						if (nodestate == -2) {
							Object d = node.getMetaData(typeLabel.get());

							if (d instanceof Integer) nodestate = (Integer) node.getMetaData(typeLabel.get());
							else if
							(d instanceof Double) nodestate = (((Double) node.getMetaData(typeLabel.get())).intValue());
							else if
							(d instanceof int[]) nodestate = (((int[]) node.getMetaData(typeLabel.get()))[0]);
						}

						return nodestate;

			}
			else return nodeStates[node.getNr()];

		}catch(Exception e){
			throw new ConstraintViolatedException("Something went wrong with the assignment of types to the nodes (node ID="+node.getID()+"). Please check your XML file!");
		}
	}

	/**
	 * Implementation of calculateSubtreeLikelihood with Small Number structure. Avoids underflowing of integration results.
	 * WARNING: calculateSubTreeLikelihood and calculateSubTreeLikelihoodSmalNumber are very similar. A modification made in one of the two would likely be needed in the other one also
	 * @param node
	 * @param from
	 * @param to
	 * @return
	 */
	p0ge_InitialConditions calculateSubtreeLikelihood(Node node, double from, double to) {

		double[] pconditions = new double[n];
		SmallNumber[] gconditions = new SmallNumber[n];
		for (int i=0; i<n; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int index = Utils.index(to,times, totalIntervals);

		if (node.isLeaf()){ // sampling event

			int nodestate = getNodeState(node, false);

			if (nodestate==-1) { //unknown state

				if (SAModel)
					throw new ConstraintViolatedException("SA model not implemented with unknown states!");

				for (int i=0; i<n; i++) {

					if (!isRhoTip[node.getNr()]) {
						init.conditionsOnG[i] = new SmallNumber(psi[i * totalIntervals + index]);
					}
					else
						init.conditionsOnG[i] = new SmallNumber(rho[i*totalIntervals+index]);
				}
			}
			else {

				if (!isRhoTip[node.getNr()]) {

					init.conditionsOnG[nodestate] = SAModel?
							new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]*(1-r[nodestate * totalIntervals + index]))
									*psi[nodestate * totalIntervals + index]) // with SA: ψ_i(r + (1 − r)p_i(τ))
							: new SmallNumber(psi[nodestate * totalIntervals + index]);


				}	else {

					//TO DO make the modif in the manuscript (for the "/(1-rho)" thing)
					init.conditionsOnG[nodestate] = SAModel? 
							new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]/(1-rho[nodestate*totalIntervals+index])*(1-r[nodestate * totalIntervals + index]))
									*rho[nodestate*totalIntervals+index])  :
							new SmallNumber(rho[nodestate*totalIntervals+index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ+δ)) //the +δ is translated by dividing p_i with 1-ρ_i (otherwise there's one too many "*ρ_i" )
					
				}

			}
			if (print) System.out.println("Sampling at time " + (T-to));

			return getG(from, init, to, node);
		}


		else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

			if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

				if (r==null)
					throw new ConstraintViolatedException("Error: Sampled ancestor found, but removalprobability not specified!");

				int childIndex = 0;

				if (node.getChild(childIndex).isDirectAncestor()) childIndex = 1;

				p0ge_InitialConditions g = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				//				int saNodeState = getNodeState(node.getChild(Math.abs(childIndex - 1)), false); // get state of direct ancestor // TO DO REMOVE if below works
				int saNodeState = getNodeState(node.getChild(childIndex ^ 1), false); // get state of direct ancestor, XOR operation gives 1 if childIndex is 0 and vice versa

				if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

					init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState];
					init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(psi[saNodeState * totalIntervals + index]
							* (1-r[saNodeState * totalIntervals + index]));

//					System.out.println("SA but not rho sampled");

				} else {
					// TO DO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
					init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState]*(1-rho[saNodeState*totalIntervals+index]) ;
					init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(rho[saNodeState*totalIntervals+index] 
							* (1-r[saNodeState * totalIntervals + index]));
					
					//TO DO working on below, probably doesn't work
//					init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState];
//					init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(rho[saNodeState*totalIntervals+index]/(1-rho[saNodeState*totalIntervals+index])
//							* (1-r[saNodeState * totalIntervals + index]));
					
//					System.out.println("SA and rho sampled and rho is: " + rho[saNodeState*totalIntervals+index] );
				}
			}

			else {   // birth / infection event

				int childIndex = 0;
				if (node.getChild(1).getNr() > node.getChild(0).getNr())
					childIndex = 1; // always start with the same child to avoid numerical differences

				p0ge_InitialConditions g0 = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				childIndex = Math.abs(childIndex - 1);

				p0ge_InitialConditions g1 = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

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

					// TO DO actually test this works with a tree with rho sampling at a branching event
					//TO DO i should be able to remove this
					//TO DO TAKE INTO ACCOUNT THE MODIFS IN THE BDMM MANUSCRIPT


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

		return getG(from, init, to, node);
	}

	public void transformParameters(){

		transformWithinParameters();
		transformAmongParameters();
	}

	// used to indicate that the state assignment went wrong
	protected class ConstraintViolatedException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public ConstraintViolatedException(String s) {
			super(s);
		}

	}

}

