package beast.evolution.speciation;

import gnu.*;
import gnu.trove.map.hash.TDoubleObjectHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.lang3.ArrayUtils;

import beast.core.Description;
import beast.core.util.Utils;
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

	TDoubleObjectHashMap<double[]> pInitialConditions;

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

				if (Math.abs(time-tipTime) < 1e-20 && rho[((MultiTypeNode)tip).getNodeType()*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

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
	public p0ge_InitialConditions getGSmallNumber(double t, p0ge_InitialConditions PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

		if (node.isLeaf()){

			// TO DO remove this part
			//			if (Math.abs(node.getHeight() - (T-t0))>1e-8) {
			//				System.out.println("Yes it is different indeed: " + "t\t" + t + " t0\t" + t0 + " node heigth\t" + (T -  node.getHeight()));
			//				if (node.getHeight() < t0)
			//					System.out.println("Yes it is indeed above");
			//				else
			//					System.out.println("No it is not above");
			//			}
			
			
			// TO DO CLEAN UP THIS SECTION
			
			// Temporary fix to pass the unit tests, it is probably not enough to actually make the code run properly
			// (because this bug should not be there in the first place, if it there probably is sth wrong someplace else) 
			if(m_rho.get()!=null && !isRhoTip[node.getNr()]){
				//if (!isRhoTip[node.getNr()]){
				//System.out.println("Yes it happened\t" + node.getNr());
				double start = System.nanoTime();
				System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0.conditionsOnP, 0, n);
				double timeElapsed = System.nanoTime() - start;
				//System.out.println("Time needed to recalculate: " + timeElapsed/1000000 + " milliseconds");

//				// TO DO remove this part
//				double[] getParray = pInitialConditions.get(t0);
//				double[] getParrayLong = PG.getP(t0, m_rho.get()!=null, rho);
//
//				// TO DO USE THIS PART TO KNOW WHERE THE NEXT BUG IS
//				for(int i=0; i<getParray.length ;i++){
//					if(Math.abs(getParrayLong[i] - getParray[i]) > 1e-5) {
//						System.out.println("Yeah the arrays are not the same dude, node number: " + node.getNr());
//						System.out.println("Pause");
//					}
//				}

			} else {
				double start = System.nanoTime();
				System.arraycopy(pInitialConditions.get(t0), 0, PG0.conditionsOnP, 0, n);
				double timeElapsed = System.nanoTime() - start;
				//System.out.println("Time needed to fetch result: " + timeElapsed/1000000 + " milliseconds");
			}


			// TO DO remove this part
			//			double[] getParray = PG.getP(t0, m_rho.get()!=null, rho);
			//
			//			// TO DO USE THIS PART TO KNOW WHERE THE NEXT BUG IS
			//			for(int i=0; i<getParray.length ;i++){
			//				if(Math.abs(pInitialConditions.get(t0)[i] - getParray[i]) > 1e-5) {
			//					System.out.println("Yeah the arrays are not the same dude, node number: " + node.getNr());
			//					System.out.println("Pause");
			//				}
			//			}
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

		//TO DO rajouter le try/catch
		//try{  // start calculation

		pInitialConditions = getAllInitialConditionsForP(tree);

		if (conditionOnSurvival.get()) {

			// TO DO remove comment
			//noSampleExistsProp = PG.getP(0,m_rho.get()!= null,rho);
			// COULD TRIGGER AN EXCEPTION
			//				for (int i=0; i<n; i++) {
			//					noSampleExistsProp[i] = pInitialConditions.get(0.)[i];
			//				}
			//				
			System.arraycopy(pInitialConditions.get(0.), 0, noSampleExistsProp, 0, n);


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

		//		}catch(Exception e){
		//			logP =  Double.NEGATIVE_INFINITY;
		//			return logP;
		//		}

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

			return getGSmallNumber(from, init, to, coltree.getRoot());
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

					// TO DO remove this part
					if (node.getHeight() != to) {

					}

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

			return getGSmallNumber(from, init, to, node);
		}

		else {

			if (migIndex==null &&  ((MultiTypeNode)node).getChangeCount()>0){ // node has migration event(psi)

				return calculateSubtreeLikelihoodSmallNumber(node, true, ((MultiTypeNode)node).getChangeCount()-1, from, to) ;
			}

			else{

				if (node.isLeaf()){ // sampling event

					if (!isRhoTip[node.getNr()]){
						// TO DO remove comments
						//												init.conditionsOnG[nodestate] = SAModel
						//														? new SmallNumber((r[nodestate * totalIntervals + index] + PG.getP(to, m_rho.get()!=null, rho)[nodestate]*(1-r[nodestate * totalIntervals + index]))
						//																*psi[nodestate * totalIntervals + index])
						//						
						//																: new SmallNumber(psi[nodestate * totalIntervals + index]);

						// TO DO CEST CA QUI BUG !!!!!!!!
						// TO DO CHECKER CA
						double[] test = pInitialConditions.get(T-node.getHeight());
						double[] testOrig = PG.getP(to, m_rho.get()!=null, rho);
						double tobis = T-node.getHeight();

						// TO DO remove this part
						//						if (Math.abs(test[nodestate]-testOrig[nodestate])> 1e-7) {
						//							System.out.println(to);
						//							throw new RuntimeException("Wouh there was an issue here");
						//						}
						init.conditionsOnG[nodestate] = SAModel
								? new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions.get(T - node.getHeight())[nodestate]*(1-r[nodestate * totalIntervals + index]))
										*psi[nodestate * totalIntervals + index])

										: new SmallNumber(psi[nodestate * totalIntervals + index]);
								//								 TENTATIVE TO DO REMOVE CE QUI SERT A RIEN
								//								init.conditionsOnG[nodestate] = SAModel
								//										? new SmallNumber((r[nodestate * totalIntervals + index] + pInitialConditions.get(node.getHeight())[nodestate]*(1-r[nodestate * totalIntervals + index]))
								//												*psi[nodestate * totalIntervals + index])
								//
								//												: new SmallNumber(psi[nodestate * totalIntervals + index]);

					} else {
						init.conditionsOnG[nodestate] = new SmallNumber(rho[nodestate*totalIntervals+index]);
					}

					if (print) System.out.println("Sampling at time " + to);

					return getGSmallNumber(from, init, to, node);
				}

				else if (node.getChildCount()==2){  // birth / infection event

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
				}
			}
		}

		return getGSmallNumber(from, init, to, node);
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


	// TO DO CLEAN UP METHOD ONCE IT WORKS
	//  TO DO REMOVE timeCount if it is really useless ( it should be)
	public TDoubleObjectHashMap<double[]> getAllInitialConditionsForP(TreeInterface tree){
		int leafCount = tree.getLeafNodeCount();
		// TO DO find a better alternative than SortedSet
		// should be able to go back to heap and heapsort as in uncoloured trees
		SortedSet<Double> breakingPointsTimes = new TreeSet<Double>();

		// TO DO remove timeCount, was only for debugging
		int timeCount = 0;
		for (int i=0; i<leafCount; i++){
			breakingPointsTimes.add(T - tree.getNode(i).getHeight());
			timeCount ++;
			int migrationCount = ((MultiTypeNode) tree.getNode(i)).getChangeCount();

			// TO DO: comment/uncomment section to see if migration events are called by getP
			if(migrationCount>0){ //if there are some migration events on the branch above this leaf, get the migration times
				// TO DO REMOVE COMMENT
				//				 VERIFIER COMMENT EST CALCULE LA HAUTEUR DU CHANGE TIME
				for(int j = 0; j<migrationCount; j++) {
					// TO DO REMOVE THIS PART
					Node n = tree.getNode(i);

					breakingPointsTimes.add(T - ((MultiTypeNode) tree.getNode(i)).getChangeTime(j));
					timeCount ++;
				}
			}

			// TO DO REMOVE THIS PART
			//			// System.out.println(nodeHeight[i]);
			//			indicesSortedByLeafHeight[i] = i;
		}

		Double[] breakingPoints = breakingPointsTimes.toArray(new Double[breakingPointsTimes.size()]);

		int nTimes = breakingPoints.length;
		double t = breakingPoints[nTimes-1];

		boolean rhoSampling =  (m_rho.get()!=null);

		TDoubleObjectHashMap<double[]>	 pInitialCondsAtLeaves = new TDoubleObjectHashMap<double[]>();

		double[] tempInitialConditions = PG.getP(t, rhoSampling, rho);

		pInitialCondsAtLeaves.put(t, tempInitialConditions); //ArrayUtils.toObject() makes the conversion from double[] to Double[]. As HAshMap cannot be used with primitive type double.
		double t0 = t;

		if (nTimes> 1){
			// TO DO faire un do while pour rendre plus joli le fait que on commence et ensuite on fait la boucle for
			for (int i = nTimes-2; i>-1; i--){
				t = breakingPoints[i];




				//If the next higher breakpoint is actually at the same height, skip iteration
				// WARNING I JUST CHANGED THE THRESHOLD FROM -20 to -50
				if (Math.abs(t-t0) < globalPrecisionThreshold) {
					t0 = t;
					continue;
				} else {
					// TO DO CLEAN UP THIS SECTION
					// double[] benchmarkValue2 = PG.getP(t0, rhoSampling, rho);
					tempInitialConditions = PG.getP(t, tempInitialConditions, t0, rhoSampling, rho);
					
					// TO DO remove if, only for debugging
					//				if((t - 0.6260071807800001)<1e-10){
					//					System.out.println(rhoSampling);
					//					System.out.println(rho.toString());
//					double[] benchmarkValue = PG.getP(t, rhoSampling, rho);
//					System.out.println(benchmarkValue[0]);
//					
//					if (benchmarkValue[0] != tempInitialConditions[0]) {
//						System.out.println("Issue here! ");
//					}
					//	
					pInitialCondsAtLeaves.put(t, tempInitialConditions);
					t0 = t;
					
					
				}
			}
		}

		if(!pInitialCondsAtLeaves.containsKey(0.)) {
			pInitialCondsAtLeaves.put(0., PG.getP(0, tempInitialConditions, t0, rhoSampling, rho));
		}


		//        String sortedIndices = new String();
		//        for (int d: indicesSortedByNodeHeight) {
		//        	sortedIndices += (d + "\t");
		//        }
		//        System.out.println("Node indices sorted: " + sortedIndices);
		//        for (int i = 0; i<nodeCount; i++) {
		//        	System.out.println("Value of node " + i + ":\t" + allPInitials[i][0]);
		//        }

		return pInitialCondsAtLeaves;
	}



	// working copy of the function for uncoloured trees
	//	/**
	//	 * Find all initial conditions for all future integrations on p0 equations 
	//	 * @param tree
	//	 * @return an array of arrays storing the initial conditions values
	//	 */
	//	public double[][] getAllInitialConditionsForP(TreeInterface tree){
	//        int nodeCount = tree.getNodeCount();
	//        double[] nodeHeights = new double[nodeCount];
	//        int[] indicesSortedByNodeHeight  =new int[nodeCount];
	//        for (int i=0; i<nodeCount; i++){
	//        	nodeHeights[i] = T - tree.getNode(i).getHeight();
	//        	// System.out.println(nodeHeight[i]);
	//        	indicesSortedByNodeHeight[i] = i;
	//        }
	//        
	//        HeapSort.sort(nodeHeights, indicesSortedByNodeHeight);
	//        //"sort" sorts in ascending order, so we have to be careful since the integration starts from the leaves at height T and goes up to the root at height 0 (or >0)
	//        
	//        int allPlength = (nodeHeights[indicesSortedByNodeHeight[0]] == 0)? nodeCount: nodeCount + 1; // in case the origin is not at zero, an extra space is left for the value of integration till 0. 
	//        double[][] allPInitials = new double[allPlength][n]; 
	//        
	//        double t = nodeHeights[indicesSortedByNodeHeight[nodeCount-1]];
	//        
	//        boolean rhoSampling =  (m_rho.get()!=null);
	//        
	//        allPInitials[indicesSortedByNodeHeight[nodeCount-1]] = PG.getP(t, rhoSampling, rho);
	//        double t0 = t;
	//        
	//        if (nodeCount >1 ){
	//            for (int i = nodeCount-2; i>-1; i--){
	//            	t = nodeHeights[indicesSortedByNodeHeight[i]];
	//            	
	//            	//If the next higher node is actually at the same height, store previous results and skip iteration
	//            	if (Math.abs(t-t0) < 1e-10) {
	//            		t0 = t;
	//            		allPInitials[indicesSortedByNodeHeight[i]] = allPInitials[indicesSortedByNodeHeight[i+1]];
	//            		continue;
	//            	} else {
	//            		allPInitials[indicesSortedByNodeHeight[i]] = PG.getP(t, allPInitials[indicesSortedByNodeHeight[i+1]], t0, rhoSampling, rho);
	//            		t0 = t;
	//            	}
	//                
	//            }
	//        }
	//        
	//        if (allPlength > nodeCount) {
	//        	allPInitials[nodeCount] = PG.getP(0, allPInitials[indicesSortedByNodeHeight[0]], t0, rhoSampling, rho);
	//        }
	//
	////        String sortedIndices = new String();
	////        for (int d: indicesSortedByNodeHeight) {
	////        	sortedIndices += (d + "\t");
	////        }
	////        System.out.println("Node indices sorted: " + sortedIndices);
	////        for (int i = 0; i<nodeCount; i++) {
	////        	System.out.println("Value of node " + i + ":\t" + allPInitials[i][0]);
	////        }
	//
	//		return allPInitials;
	//	}
}
