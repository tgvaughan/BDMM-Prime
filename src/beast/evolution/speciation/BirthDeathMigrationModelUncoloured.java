package beast.evolution.speciation;

import beast.evolution.tree.*;
import beast.core.Input;
import beast.core.Description;
import beast.core.util.Utils;

import beast.math.SmallNumber;
import beast.math.p0ge_InitialConditions;
import beast.util.HeapSort;


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

	//TO DO REMOVE IF USELESS
	/*
	void computeRhoInternalNodes(){
		double nodeTime;
		int tipCount = treeInput.get().getLeafNodeCount();

		for (Node internalNode : treeInput.get().getInternalNodes()) {

			nodeTime = T-internalNode.getHeight();
			isRhoInternalNode[internalNode.getNr()-tipCount] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TO DO: make a warning that rho sampling precision is with 1e-10. Maybe do a relative threshold?
				if (Math.abs(time-nodeTime) < globalPrecisionThreshold && rho[getNodeState(internalNode,false)*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoInternalNode[internalNode.getNr()-tipCount] = true;

			}
		}
	}
	*/

	/**
	 * WARNING: getG and getGSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param t
	 * @param PG0
	 * @param t0
	 * @param node
	 * @return
	 */
	public double[] getG(double t, double[] PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		if (node.isLeaf()) {
			System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0, 0, n);
		}

		return getG(t,  PG0,  t0, pg_integrator, PG, T, maxEvalsUsed);

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
	public p0ge_InitialConditions getGSmallNumber(double t, p0ge_InitialConditions PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)


		if (node.isLeaf()) {

			//			// TO DO CLEAN UP
			//			//System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0.conditionsOnP, 0, n);
			//			double h = T - node.getHeight();
			//			double[] temp = PG.getP(t0, m_rho.get()!=null, rho);
			//			double[] temp2 = pInitialConditions[node.getNr()];
			//			
			//			if (h!=t0) {
			//				throw new RuntimeException("t0 est pas comme height");
			//			}
			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, n);
			//System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0.conditionsOnP, 0, n);
		}

		return getGSmallNumber(t,  PG0,  t0, pg_integrator, PG, T, maxEvalsUsed);

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
		double Pr = 0;
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

			// alternatively p or pSN will be used, depending if the user wants to use the classic implementation or the one with SmallNumbers
			p0ge_InitialConditions pSN = new p0ge_InitialConditions();
			double [] p = new double[] {0};

			if ( orig > 0 ) {
				if (useSmallNumbers.get())
					pSN = calculateSubtreeLikelihoodSmallNumber(root,0,orig);
				else
					p = calculateSubtreeLikelihood(root,0,orig);
			} else {
				int childIndex = 0;
				if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences


				// depending on the method chosen by the user, SmallNumbers or doubles will be used to store the results between calculations step
				if (useSmallNumbers.get()) {
					pSN = calculateSubtreeLikelihoodSmallNumber(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					childIndex = Math.abs(childIndex-1);

					p0ge_InitialConditions p1SN = calculateSubtreeLikelihoodSmallNumber(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					for (int i =0; i<pSN.conditionsOnG.length; i++) pSN.conditionsOnG[i] = SmallNumber.multiply(pSN.conditionsOnG[i], p1SN.conditionsOnG[i]);

				} else { 
					p = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					childIndex = Math.abs(childIndex-1);

					double [] p1 = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

					for (int i =0; i<p.length; i++) p[i]*=p1[i];
				}
			}

			if (print) System.out.print("final p per state = ");

			if (useSmallNumbers.get()) {
				for (int root_state=0; root_state<n; root_state++){

					if (pSN.conditionsOnG[root_state].getMantissa()>0 ) 
						PrSN = SmallNumber.add(PrSN, pSN.conditionsOnG[root_state].scalarMultiply(freq[root_state]));

					if (print) System.out.print(pSN.conditionsOnP[root_state] + "\t" + pSN.conditionsOnG[root_state] + "\t");
				}

				if (conditionOnSurvival.get()){
					PrSN = PrSN.scalarMultiply(1/(1-nosample));
				}

			} else {
				for (int root_state=0; root_state<n; root_state++){

					if (p[n+root_state]>0 ) Pr += freq[root_state]* p[n+root_state];

					if (print) System.out.print(p[root_state] + "\t" + p[root_state+n] + "\t");
				}

				if (conditionOnSurvival.get()){
					Pr /= (1-nosample);
				}
			}



		}catch(Exception e){

			if (e instanceof ConstraintViolatedException){throw e;}

			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		maxEvalsUsed = Math.max(maxEvalsUsed, PG.maxEvalsUsed);

		if (useSmallNumbers.get()) logP = PrSN.log();

		else logP = Math.log(Pr);

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
	 * WARNING: calculateSubTreeLikelihood and calculateSubTreeLikelihoodSmalNumber are very similar. A modification made in one of the two would likely be needed in the other one also
	 * @param node
	 * @param from
	 * @param to
	 * @return
	 */
	double[] calculateSubtreeLikelihood(Node node, double from, double to) {


		double[] init = new double[2*n];

		int index = Utils.index(to,times, totalIntervals);

		if (node.isLeaf()){ // sampling event

			int nodestate = getNodeState(node, false);

			if (nodestate==-1) { //unknown state

				if (SAModel)
					throw new ConstraintViolatedException("SA model not implemented with unknown states!");

				for (int i=0; i<n; i++) {

					if (!isRhoTip[node.getNr()]) {
						init[n + i] = psi[i * totalIntervals + index];
					}
					else
						init[n + i] = rho[i*totalIntervals+index];
				}
			}
			else {

				if (!isRhoTip[node.getNr()])
					init[n + nodestate] = SAModel
					? psi[nodestate * totalIntervals + index]* (r[nodestate * totalIntervals + index] + (1-r[nodestate * totalIntervals + index])*PG.getP(to, m_rho.get()!=null, rho)[nodestate]) // with SA: ψ_i(r + (1 − r)p_i(τ))
							: psi[nodestate * totalIntervals + index];

					else
						init[n+nodestate] = rho[nodestate*totalIntervals+index];

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

				double[] g = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				int saNodeState = getNodeState(node.getChild(Math.abs(childIndex - 1)), false); // get state of direct ancestor

				init[saNodeState] = g[saNodeState];
				//initial condition for SA: ψ_i(1 − r_i) g :
				init[n + saNodeState] = psi[saNodeState * totalIntervals + index] * (1-r[saNodeState * totalIntervals + index]) * g[n + saNodeState];
			}

			else {   // birth / infection event

				int childIndex = 0;
				if (node.getChild(1).getNr() > node.getChild(0).getNr())
					childIndex = 1; // always start with the same child to avoid numerical differences

				double[] g0 = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				childIndex = Math.abs(childIndex - 1);

				double[] g1 = calculateSubtreeLikelihood(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				if (print)
					System.out.println("Infection at time " + (T - to));//+ " with p = " + p + "\tg0 = " + g0 + "\tg1 = " + g1);


				for (int childstate = 0; childstate < n; childstate++) {

					if (print) {
						System.out.println("state " + childstate + "\t p0 = " + g0[childstate] + "\t p1 = " + g1[childstate]);
						System.out.println("\t\t g0 = " + g0[n + childstate] + "\t g1 = " + g1[n + childstate]);
					}

					init[childstate] = g0[childstate]; //Math.floor(((g0[childstate]+g1[childstate])/2.)/tolerance.get())*tolerance.get();  // p0 is the same for both sides of the tree, but there might be tiny numerical differences
					init[n + childstate] = birth[childstate * totalIntervals + index] * g0[n + childstate] * g1[n + childstate];

					if (birthAmongDemes) {
						for (int j = 0; j < n; j++) {
							if (childstate != j) {

								init[n + childstate] += 0.5 * b_ij[totalIntervals * (childstate * (n - 1) + (j < childstate ? j : j - 1)) + index] * (g0[n + childstate] * g1[n + j] + g0[n + j] * g1[n + childstate]);
							}
						}
					}

					if (Double.isInfinite(init[childstate])) {
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
			for (int i=0;i<2*n;i++) System.out.print(init[i] + "\t");
			System.out.println();
		}

		return getG(from, init, to, node);
	}

	/**
	 * Implementation of calculateSubtreeLikelihood with Small Number structure. Avoids underflowing of integration results.
	 * WARNING: calculateSubTreeLikelihood and calculateSubTreeLikelihoodSmalNumber are very similar. A modification made in one of the two would likely be needed in the other one also
	 * @param node
	 * @param from
	 * @param to
	 * @return
	 */
	p0ge_InitialConditions calculateSubtreeLikelihoodSmallNumber(Node node, double from, double to) {

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

					// TO DO change Threshold of 1e-10 when threshold in computeRhoTips is changed
					if((node.getHeight())< 1e-10) 
						init.conditionsOnG[nodestate] = new SmallNumber(rho[nodestate*totalIntervals+index]);
					else {
						init.conditionsOnG[nodestate] = SAModel? 
								new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]*(1-r[nodestate * totalIntervals + index]))
										*rho[nodestate*totalIntervals+index]):
											new SmallNumber(rho[nodestate*totalIntervals+index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ))

					}
				}


			}
			if (print) System.out.println("Sampling at time " + (T-to));

			return getGSmallNumber(from, init, to, node);
		}


		else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

			if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

				if (r==null)
					throw new ConstraintViolatedException("Error: Sampled ancestor found, but removalprobability not specified!");

				int childIndex = 0;

				if (node.getChild(childIndex).isDirectAncestor()) childIndex = 1;

				p0ge_InitialConditions g = calculateSubtreeLikelihoodSmallNumber(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				int saNodeState = getNodeState(node.getChild(Math.abs(childIndex - 1)), false); // get state of direct ancestor

				init.conditionsOnP[saNodeState] = g.conditionsOnP[saNodeState];
				init.conditionsOnG[saNodeState] = g.conditionsOnG[saNodeState].scalarMultiply(psi[saNodeState * totalIntervals + index] * (1-r[saNodeState * totalIntervals + index]));
			}

			else {   // birth / infection event

				int childIndex = 0;
				if (node.getChild(1).getNr() > node.getChild(0).getNr())
					childIndex = 1; // always start with the same child to avoid numerical differences

				p0ge_InitialConditions g0 = calculateSubtreeLikelihoodSmallNumber(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

				childIndex = Math.abs(childIndex - 1);

				p0ge_InitialConditions g1 = calculateSubtreeLikelihoodSmallNumber(node.getChild(childIndex), to, T - node.getChild(childIndex).getHeight());

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

		return getGSmallNumber(from, init, to, node);
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

