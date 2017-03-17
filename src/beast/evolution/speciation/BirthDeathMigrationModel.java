package beast.evolution.speciation;

import beast.core.Description;
import beast.core.util.Utils;
import beast.evolution.speciation.BirthDeathMigrationModelUncoloured.ConstraintViolatedException;
import beast.evolution.tree.*;
import beast.core.Input;

import beast.math.SmallNumber;
import beast.math.p0ge_InitialConditions;
import beast.util.HeapSort;



/**
 * @author Denise Kuehnert
 *         Date: May 25, 2012
 *         Time: 11:38:27 AM
 */

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
		"This should only be used when the migration process along the phylogeny is important. Otherwise the computationally less intense BirthDeathMigrationModelUncoloured can be employed." +
		"Two implementations are available. The first is the fast classic one; the second one prevents underflowing, using so-called 'SmallNumbers', with the cost of additional computational complexity")
public class BirthDeathMigrationModel extends PiecewiseBirthDeathMigrationDistribution {

	public Input<MultiTypeRootBranch> originBranchInput =
			new Input<>("originBranch", "MultiTypeRootBranch for origin coloring");



	MultiTypeTree coltree;
	MultiTypeRootBranch originBranch;

	double[][] pInitialConditionsOnMigrationEvents;

	Boolean print = false;

	@Override
	public void initAndValidate() {

		super.initAndValidate();

		if (birthRateAmongDemes.get() !=null || R0AmongDemes.get()!=null)
			throw new RuntimeException("Error: You've specified birthRateAmongDemes or R0AmongDemes, but transmission among demes is currently not possible in MultiTypeTrees. " +
					"Please use BirthDeathMigrationModelUncoloured instead.");

		coltree = (MultiTypeTree) treeInput.get();

		if (origin.get()==null){

			T = coltree.getRoot().getHeight();
		}
		else {

			originBranch = originBranchInput.get();

			if (originBranch==null)  throw new RuntimeException("Error: Origin specified but originBranch missing!");

			checkOrigin(coltree);
		}

		ntaxa = coltree.getLeafNodeCount();

		int contempCount = 0;
		for (Node node : coltree.getExternalNodes())
			if (node.getHeight()==0.)
				contempCount++;
		if (checkRho.get() && contempCount>1 && rho==null)
			throw new RuntimeException("Error: multiple tips given at present, but sampling probability \'rho\' is not specified.");

		collectTimes(T);
		setRho();

	}

	double updateRates(){

		birth = new Double[n*totalIntervals];
		death = new Double[n*totalIntervals];
		psi = new Double[n*totalIntervals];
		M = new Double[totalIntervals*(n*(n-1))];
		if (SAModel) r =  new Double[n * totalIntervals];

		if (transform)
			transformParameters();

		else
			updateBirthDeathPsiParams();

		Double[] migRates = migrationMatrix.get().getValues();

		updateAmongParameter(M, migRates, migChanges, migChangeTimes);

		updateRho();

		freq = frequencies.get().getValues();

		setupIntegrators();

		return 0.;
	}

	@Override
	void computeRhoTips(){

		double tipTime;

		for (Node tip : treeInput.get().getExternalNodes()) {

			tipTime = T-tip.getHeight();
			isRhoTip[tip.getNr()] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TO DO: DEAL WITH THE IMPLICIT THRESHOLD HERE, THAT SHOULD WORK WITH THE OTHERS (and probably same thing in coloured)
				// may need to change to 1e-10
				if (Math.abs(time-tipTime) < 1e-10 && rho[((MultiTypeNode)tip).getNodeType()*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

			}
		}
	}

	void computeRhoInternalNodes(){
		double nodeTime;
		int tipCount = treeInput.get().getLeafNodeCount();

		for (Node internalNode : treeInput.get().getInternalNodes()) {

			nodeTime = T-internalNode.getHeight();
			isRhoInternalNode[internalNode.getNr()-tipCount] = false;

			for (Double time:rhoSamplingChangeTimes){

				// TO DO: DEAL WITH THE IMPLICIT THRESHOLD HERE, THAT SHOULD WORK WITH THE OTHERS (and probably same thing in coloured)
				if (Math.abs(time-nodeTime) < 1e-10 && rho[((MultiTypeNode)internalNode).getNodeType()*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoInternalNode[internalNode.getNr()-tipCount] = true;

			}
		}
	}

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
	public p0ge_InitialConditions getGSmallNumber(double t, p0ge_InitialConditions PG0, double t0, Node node, boolean isMigrationEvent){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		if (node.isLeaf() && !isMigrationEvent){
			//			// TO DO CLEAN UP
						//System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0.conditionsOnP, 0, n);
//						double h = T - node.getHeight();
//						double[] temp = PG.getP(t0, m_rho.get()!=null, rho);
//						double[] temp2 = pInitialConditions[node.getNr()];
//						
//						if (h!=t0) {
//							throw new RuntimeException("t0 est pas comme height");
//						}
//TO DO REMOVE IF IT WORKS
			//System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0.conditionsOnP, 0, n);
			System.arraycopy(pInitialConditions[node.getNr()], 0, PG0.conditionsOnP, 0, n);
		}

		return getGSmallNumber(t,  PG0,  t0, pg_integrator, PG, T, maxEvalsUsed);
	}


	/**
	 * WARNING: calculateTreeLogLikelihood allows use of both classic and non-underflowing methods. Some chunks of code are therefore present in two similar versions in this method.
	 * When modifying one of the two versions, one should check if the other version also need the corresponding changes.
	 */
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {

		if (SAModel && treeInput.isDirty()) throw new RuntimeException("Error: SA Model only implemented for fixed trees!");

		coltree = (MultiTypeTree) tree;

		MultiTypeNode root = (MultiTypeNode) coltree.getRoot();


		if (!coltree.isValid() || (origin.get()!=null && !originBranchIsValid(root))){
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		int node_state;
		if (origin.get()==null) {
			T = root.getHeight();
			node_state =  ((MultiTypeNode) coltree.getRoot()).getNodeType(); 
		}
		else{
			updateOrigin(root);
			node_state = (originBranch.getChangeCount()>0) ? originBranch.getChangeType(originBranch.getChangeCount()-1) : ((MultiTypeNode) coltree.getRoot()).getNodeType();

			if (orig < 0){
				return Double.NEGATIVE_INFINITY;
			}
		}

		collectTimes(T);
		setRho();

		if (updateRates() < 0 ||  (times[totalIntervals-1] > T)) { 
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		double[] noSampleExistsProp =  new double[n];

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

			// alternatively p or pSN will be used, depending on whether the user wants to use the classic implementation or the one with SmallNumbers
			p0ge_InitialConditions pSN = new p0ge_InitialConditions();
			double [] p = new double[] {0};

			if (orig>0){
				if (originBranch.getChangeCount()>0) {
					if (useSmallNumbers.get())
						pSN = calculateOriginLikelihoodSmallNumber(originBranch.getChangeCount()-1, 0, T-originBranch.getChangeTime(originBranch.getChangeCount()-1) );
					else
						p = calculateOriginLikelihood(originBranch.getChangeCount()-1, 0, T-originBranch.getChangeTime(originBranch.getChangeCount()-1) );
				} else {
					if (useSmallNumbers.get())
						pSN = calculateSubtreeLikelihoodSmallNumber(root, false, null, 0, orig);
					else
						p = calculateSubtreeLikelihood(root, false, null, 0, orig);
				}
			} else {
				int childIndex = 0;
				if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

				double t0 = T - root.getChild(childIndex).getHeight();
				int childChangeCount = ((MultiTypeNode)root.getChild(childIndex)).getChangeCount();
				if (childChangeCount > 0)
					t0 = T - ((MultiTypeNode)root.getChild(childIndex)).getChangeTime(childChangeCount-1);

				// SmallNumbers or doubles will be used to store the results between calculations step, depending on the method chosen by the user
				if (useSmallNumbers.get())
					pSN = calculateSubtreeLikelihoodSmallNumber(root.getChild(childIndex), false, null, 0., t0);
				else
					p = calculateSubtreeLikelihood(root.getChild(childIndex), false, null, 0., t0);

				childIndex = Math.abs(childIndex-1);

				t0 = T - root.getChild(childIndex).getHeight();
				childChangeCount = ((MultiTypeNode)root.getChild(childIndex)).getChangeCount(); // changeCounts[root.getChild(1).getNr()];
				if (childChangeCount > 0)
					t0 = T - ((MultiTypeNode)root.getChild(childIndex)).getChangeTime(childChangeCount-1);

				if (useSmallNumbers.get()) {
					p0ge_InitialConditions p1SN = calculateSubtreeLikelihoodSmallNumber(root.getChild(childIndex), false, null, 0., t0);

					for (int i=0; i<pSN.conditionsOnG.length; i++) pSN.conditionsOnG[i] = SmallNumber.multiply(pSN.conditionsOnG[i], p1SN.conditionsOnG[i]);

				} else {
					double[] p1 = calculateSubtreeLikelihood(root.getChild(childIndex), false, null, 0., t0);

					for (int i=0; i<p.length; i++) p[i]*=p1[i];
				}
			}
			if (conditionOnSurvival.get()) {
				if (useSmallNumbers.get())
					pSN.conditionsOnG[node_state] = pSN.conditionsOnG[node_state].scalarMultiply(1/(1-noSampleExistsProp[node_state]));    // condition on survival
				else
					p[n+node_state] /= (1-noSampleExistsProp[node_state]);
			}

			if (useSmallNumbers.get())
				logP = Math.log(freq[node_state]) +  pSN.conditionsOnG[node_state].log();
			else
				logP = Math.log(freq[node_state]) +  Math.log(p[n+node_state]);


			maxEvalsUsed = Math.max(maxEvalsUsed, PG.maxEvalsUsed);

		}catch(Exception e){
			logP =  Double.NEGATIVE_INFINITY;
			return logP;
		}

		if (print) System.out.println("final logL = " + logP);

		if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

		if (SAModel && !(removalProbability.get().getDimension()==n && removalProbability.get().getValue()==1.)) {
			int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
			logP +=  Math.log(2)*internalNodeCount;
		}

		return logP;
	}

	/**
	 * WARNING: calculateOriginLikelihood and calculateOriginLikelihoodSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param migIndex
	 * @param from
	 * @param to
	 * @return
	 */
	double[] calculateOriginLikelihood(Integer migIndex, double from, double to) {

		double[] init = new double[2*n];
		int index = Utils.index(to, times, totalIntervals);

		int prevcol = originBranch.getChangeType(migIndex);
		int col =  (migIndex > 0)?  originBranch.getChangeType(migIndex-1):  ((MultiTypeNode) coltree.getRoot()).getNodeType();

		migIndex--;

		double[] g ;

		if (migIndex >= 0){

			g = calculateOriginLikelihood(migIndex, to, T - originBranch.getChangeTime(migIndex));

			System.arraycopy(g, 0, init, 0, n);
			init[n+prevcol] = M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index] * g[n + col];       // with ratechange in M

			return getG(from,  init,  to, pg_integrator, PG, T, maxEvalsUsed);

		}
		else {

			g = calculateSubtreeLikelihood(coltree.getRoot(), false, null, to, orig);

			System.arraycopy(g, 0, init, 0, n);
			init[n+prevcol] = M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index] * g[n + col];       // with ratechange in M

			return getG(from, init, to, coltree.getRoot());
		}
	}

	/**
	 * Implementation of calculateOriginLikelihood with Small Number structure. Avoids underflowing of integration results.
	 * WARNING: calculateOriginLikelihood and calculateOriginLikelihoodSmallNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param migIndex
	 * @param from
	 * @param to
	 * @return
	 */
	p0ge_InitialConditions calculateOriginLikelihoodSmallNumber(Integer migIndex, double from, double to) {

		double[] pconditions = new double[n];
		SmallNumber[] gconditions = new SmallNumber[n];
		for (int i=0; i<n; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int index = Utils.index(to, times, totalIntervals);

		int prevcol = originBranch.getChangeType(migIndex);
		int col =  (migIndex > 0)?  originBranch.getChangeType(migIndex-1):  ((MultiTypeNode) coltree.getRoot()).getNodeType();

		migIndex--;

		p0ge_InitialConditions g ;

		if (migIndex >= 0){

			g = calculateOriginLikelihoodSmallNumber(migIndex, to, T - originBranch.getChangeTime(migIndex));

			System.arraycopy(g.conditionsOnP, 0, pconditions, 0, n);
			init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);		// with ratechange in M


			return getGSmallNumber(from,  init,  to, pg_integrator, PG, T, maxEvalsUsed);

		}
		else {

			g = calculateSubtreeLikelihoodSmallNumber(coltree.getRoot(), false, null, to, orig);

			System.arraycopy(g.conditionsOnP, 0, pconditions, 0, n);
			init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]);		// with ratechange in M

			// TO DO CHECK THAT isMigrationEvent should really be set to false here (otherwise pb with getP calc in getG
			// but should be ok bc coltree.getRoot should not be a leaf (except if tree with one tip maybe, which is not very interesting) 
			return getGSmallNumber(from, init, to, coltree.getRoot(), false);
		}
	}

	/**
	 * WARNING: calculateSubTreeLikelihood and calculateSubTreeLikelihoodSmalNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param node
	 * @param migration
	 * @param migIndex
	 * @param from
	 * @param to
	 * @return
	 */
	double[] calculateSubtreeLikelihood(Node node, Boolean migration, Integer migIndex, double from, double to) {

		double[] init = new double[2*n];
		int nodestate = ((MultiTypeNode)node).getNodeType();
		int index = Utils.index(to, times, totalIntervals);

		if (migration){ // migration event

			int prevcol = ((MultiTypeNode) node).getChangeType(migIndex);
			int col =  (migIndex > 0)?  ((MultiTypeNode) node).getChangeType(migIndex-1):  ((MultiTypeNode) node).getNodeType();
			double time ;

			migIndex--;

			time = (migIndex >= 0)? ((MultiTypeNode) node).getChangeTime(migIndex) :node.getHeight();
			double[] g = calculateSubtreeLikelihood(node, (migIndex >= 0), migIndex, to, T-time);

			System.arraycopy(g, 0, init, 0, n);
			init[n+prevcol] = M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index] * g[n + col];       // with ratechange in M

			return getG(from, init, to, node);
		}

		else {

			if (migIndex==null &&  ((MultiTypeNode)node).getChangeCount()>0){ // node has migration event(psi)

				return calculateSubtreeLikelihood(node, true, ((MultiTypeNode)node).getChangeCount()-1, from, to) ;
			}

			else{

				if (node.isLeaf()){ // sampling event

					if (!isRhoTip[node.getNr()])

						init[n + nodestate] = SAModel
						? psi[nodestate * totalIntervals + index]* (r[nodestate * totalIntervals + index] + (1-r[nodestate * totalIntervals + index])*PG.getP(to, m_rho.get()!=null, rho)[nodestate]) // with SA: Ïˆ_i(r + (1 âˆ’ r)p_i(Ï„))
								: psi[nodestate * totalIntervals + index];

						else
							init[n+nodestate] = rho[nodestate*totalIntervals+index];


					if (print) System.out.println("Sampling at time " + to);

					return getG(from, init, to, node);
				}

				else if (node.getChildCount()==2){  // birth / infection event

					int childIndex = 0;
					if (node.getChild(1).getNr() > node.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

					double t0 = T - node.getChild(childIndex).getHeight();
					int childChangeCount = ((MultiTypeNode)node.getChild(childIndex)).getChangeCount();
					if (childChangeCount > 0)
						t0 = T - ((MultiTypeNode)node.getChild(childIndex)).getChangeTime(childChangeCount-1);

					double[] g0 = calculateSubtreeLikelihood(node.getChild(childIndex), false, null, to, t0);

					childIndex = Math.abs(childIndex-1);

					double t1 = T - node.getChild(childIndex).getHeight();
					childChangeCount = ((MultiTypeNode)node.getChild(childIndex)).getChangeCount();
					if (childChangeCount > 0)
						t1 = T - ((MultiTypeNode)node.getChild(childIndex)).getChangeTime(childChangeCount-1);

					double[] g1 = calculateSubtreeLikelihood(node.getChild(childIndex), false, null, to, t1);

					System.arraycopy(g0, 0, init, 0, n);
					init[n+nodestate] =  birth[nodestate*totalIntervals+index] * g0[n+nodestate] * g1[n+nodestate];
				}
			}
		}

		return getG(from, init, to, node);
	}

	/**
	 * Implementation of calculateSubtreeLikelihood with Small Number structure. Avoids underflowing of integration results.
	 * WARNING: calculateSubTreeLikelihood and calculateSubTreeLikelihoodSmalNumber are very similar. A modification made in one of the two would likely be needed in the other one also.
	 * @param node
	 * @param migration
	 * @param migIndex
	 * @param from
	 * @param to
	 * @return
	 */
	p0ge_InitialConditions calculateSubtreeLikelihoodSmallNumber(Node node, Boolean migration, Integer migIndex, double from, double to) {

		double[] pconditions = new double[n];
		SmallNumber[] gconditions = new SmallNumber[n];
		for (int i=0; i<n; i++) gconditions[i] = new SmallNumber();

		p0ge_InitialConditions init = new p0ge_InitialConditions(pconditions, gconditions);

		int nodestate = ((MultiTypeNode)node).getNodeType();
		int index = Utils.index(to, times, totalIntervals);

		if (migration){ // migration event

			int prevcol = ((MultiTypeNode) node).getChangeType(migIndex);
			int col =  (migIndex > 0)?  ((MultiTypeNode) node).getChangeType(migIndex-1):  ((MultiTypeNode) node).getNodeType();
			double time ;

			migIndex--;

			time = (migIndex >= 0)? ((MultiTypeNode) node).getChangeTime(migIndex) :node.getHeight();
			p0ge_InitialConditions g = calculateSubtreeLikelihoodSmallNumber(node, (migIndex >= 0), migIndex, to, T-time);

			System.arraycopy(g.conditionsOnP, 0, init.conditionsOnP, 0, n);
			init.conditionsOnG[prevcol] = g.conditionsOnG[col].scalarMultiply(M[totalIntervals * (prevcol * (n - 1) + (col < prevcol ? col : col - 1)) + index]); // with ratechange in M

			return getGSmallNumber(from, init, to, node, true);
		}

		else {

			if (migIndex==null &&  ((MultiTypeNode)node).getChangeCount()>0){ // node has migration event(psi)

				return calculateSubtreeLikelihoodSmallNumber(node, true, ((MultiTypeNode)node).getChangeCount()-1, from, to) ;
			}

			else{

				if (node.isLeaf()){ // sampling event

					if (!isRhoTip[node.getNr()]){
						
						init.conditionsOnG[nodestate] = SAModel
								? new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]*(1-r[nodestate * totalIntervals + index]))
										*psi[nodestate * totalIntervals + index])

										: new SmallNumber(psi[nodestate * totalIntervals + index]);
								
								//TO DO REMOVE IF ABOVE WORKS
//						init.conditionsOnG[nodestate] = SAModel
//								? new SmallNumber((r[nodestate * totalIntervals + index] + PG.getP(to, m_rho.get()!=null, rho)[nodestate]*(1-r[nodestate * totalIntervals + index]))
//										*psi[nodestate * totalIntervals + index])
//
//										: new SmallNumber(psi[nodestate * totalIntervals + index]);

					} else {
						// TO DO change Threshold of 1e-10 when threshold in computeRhoTips is changed
						if((node.getHeight())< 1e-10) 
							init.conditionsOnG[nodestate] = new SmallNumber(rho[nodestate*totalIntervals+index]);
						else
							
							init.conditionsOnG[nodestate] = SAModel? 
							new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions[node.getNr()][nodestate]*(1-r[nodestate * totalIntervals + index]))
									*rho[nodestate*totalIntervals+index]):
										new SmallNumber(rho[nodestate*totalIntervals+index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ))

							
							//TO DO REMOVE IF ABOVE WORK
//							init.conditionsOnG[nodestate] = SAModel? 
//									new SmallNumber((r[nodestate * totalIntervals + index] + PG.getP(to, m_rho.get()!=null, rho)[nodestate]*(1-r[nodestate * totalIntervals + index]))
//											*rho[nodestate*totalIntervals+index]):
//												new SmallNumber(rho[nodestate*totalIntervals+index]); // rho-sampled leaf in the past: ρ_i(τ)(r + (1 − r)p_i(τ))

					}

					if (print) System.out.println("Sampling at time " + to);

					return getGSmallNumber(from, init, to, node, false);
				}

				else if (node.getChildCount()==2){  // birth / infection event or sampled ancestor

					// TO DO make test to check that the sampled ancestor thing here works
					//NO IDEA if this part is actually reached by the code, check that
					if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

						if (r==null)
							throw new RuntimeException("Error: Sampled ancestor found, but removalprobability not specified!");

						int childIndex = 0;

						if (node.getChild(childIndex).isDirectAncestor()) childIndex = 1;

						p0ge_InitialConditions g = calculateSubtreeLikelihoodSmallNumber(node.getChild(childIndex), false, null, to, T - node.getChild(childIndex).getHeight());

						init.conditionsOnP[nodestate] = g.conditionsOnP[nodestate];
						init.conditionsOnG[nodestate] = g.conditionsOnG[nodestate].scalarMultiply(psi[nodestate * totalIntervals + index] * (1-r[nodestate * totalIntervals + index]));
					}

					else {   // birth / infection event

						int childIndex = 0;
						if (node.getChild(1).getNr() > node.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

						double t0 = T - node.getChild(childIndex).getHeight();
						int childChangeCount = ((MultiTypeNode)node.getChild(childIndex)).getChangeCount();
						if (childChangeCount > 0)
							t0 = T - ((MultiTypeNode)node.getChild(childIndex)).getChangeTime(childChangeCount-1);


						p0ge_InitialConditions g0 = calculateSubtreeLikelihoodSmallNumber(node.getChild(childIndex), false, null, to, t0);

						childIndex = Math.abs(childIndex-1);

						double t1 = T - node.getChild(childIndex).getHeight();
						childChangeCount = ((MultiTypeNode)node.getChild(childIndex)).getChangeCount(); 
						if (childChangeCount > 0)
							t1 = T - ((MultiTypeNode)node.getChild(childIndex)).getChangeTime(childChangeCount-1);

						p0ge_InitialConditions g1 = calculateSubtreeLikelihoodSmallNumber(node.getChild(childIndex), false, null, to, t1);

						System.arraycopy(g0.conditionsOnP, 0, init.conditionsOnP, 0, n);
						init.conditionsOnG[nodestate] = SmallNumber.multiply(g0.conditionsOnG[nodestate], g1.conditionsOnG[nodestate]).scalarMultiply(birth[nodestate*totalIntervals+index]);

						// TO DO actually test this works with a tree with rho sampling at a branching event
						// TO DO check that this part of the code is actually reached
						if (m_rho.get()!=null && isRhoInternalNode[node.getNr()-treeInput.get().getLeafNodeCount()]) {

							init.conditionsOnG[nodestate]= init.conditionsOnG[nodestate].scalarMultiply(1 - rho[nodestate*totalIntervals+index]);
							// TO DO REMOVE PRINT
							System.out.println("state " + nodestate + "\t 1-rho[childstate] " + (1 - rho[nodestate*totalIntervals+index]));

						}
					}
				}
			}
		}

		//TO DO: again, check that this can never be starting from a migration event, but it shouldn't
		return getGSmallNumber(from, init, to, node, false);
	}

	public void transformParameters(){

		transformWithinParameters();
	}

	public Boolean originBranchIsValid(MultiTypeNode root){

		int count = originBranch.getChangeCount();

		if (count>0){

			if (originBranch.getChangeTime(0) < root.getHeight() || originBranch.getChangeTime(count-1) > origin.get().getValue() )
				return false;

			if (originBranch.getChangeType(0) == root.getFinalType())
				return false;

			for (int i=1; i<count; i++){
				if (originBranch.getChangeType(i-1) == originBranch.getChangeType(i))
					return false;
			}
		}
		return true;
	}

}
