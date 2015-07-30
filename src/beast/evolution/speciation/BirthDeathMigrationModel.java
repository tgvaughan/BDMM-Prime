package beast.evolution.speciation;

import beast.core.Description;
import beast.core.util.Utils;
import beast.evolution.tree.*;
import beast.core.parameter.RealParameter;
import beast.core.Input;

import math.p0_ODE;
import math.p0ge_ODE;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;


/**
 * @author Denise Kuehnert
 *         Date: May 25, 2012
 *         Time: 11:38:27 AM
 */

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
        "This should only be used when the migration process along the phylogeny is important. Otherwise the computationally less intense BirthDeathMigrationModelUncoloured can be employed.")
public class BirthDeathMigrationModel extends PiecewiseBirthDeathSamplingDistribution {

    public Input<RealParameter> migrationMatrix =
            new Input<>("migrationMatrix", "Flattened migration matrix, can be asymmetric, diagnonal entries omitted",  Input.Validate.REQUIRED);

    public Input<RealParameter> frequencies =
            new Input<>("frequencies", "state frequencies",  Input.Validate.REQUIRED);

    public Input<RealParameter> origin =
            new Input<>("origin", "The origin of infection x1");

    public Input<MultiTypeRootBranch> originBranchInput =
            new Input<>("originBranch", "MultiTypeRootBranch for origin coloring");

    public Input<Boolean> originIsRootEdge =
            new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);

    public Input<Integer> maxEvaluations =
            new Input<>("maxEvaluations", "The maximum number of evaluations for ODE solver", 20000);

    public Input<Boolean> conditionOnSurvival =
            new Input<>("conditionOnSurvival", "condition on at least one survival? Default true.", true);

    public Input<Double> tolerance =
            new Input<>("tolerance", "tolerance for numerical integration", 1e-14);


    MultiTypeTree coltree;

    Double[] freq;
    Double[] M;
    double T;
    double orig;
    MultiTypeRootBranch originBranch;
    int ntaxa;

    p0_ODE P;
    p0ge_ODE PG;

    FirstOrderIntegrator pg_integrator;
    public int maxEvalsUsed;
    public Double minstep;
    public Double maxstep;

    Boolean print = false;

    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();

        coltree = (MultiTypeTree) treeInput.get();

        if (origin.get()==null){
            T = coltree.getRoot().getHeight();
        }

        else{

            originBranch = originBranchInput.get();

            if (originBranch==null)  throw new RuntimeException("Error: Origin specified but originBranch missing!");

                updateOrigin(coltree.getRoot());

            if (!Boolean.valueOf(System.getProperty("beast.resume")) && treeInput.get().getRoot().getHeight() >= origin.get().getValue())
                throw new RuntimeException("Error: origin("+T+") must be larger than tree height("+coltree.getRoot().getHeight()+")!");
        }

        ntaxa = coltree.getLeafNodeCount();

        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null){

            transform = false;
            death = deathRate.get().getValues();
            psi = samplingRate.get().getValues();
            birth = birthRate.get().getValues();
        }
        else if (R0.get() != null && becomeUninfectiousRate.get() != null && samplingProportion.get() != null){

            transform = true;
        }

        else{
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }

        M = migrationMatrix.get().getValues();

        if (n>1 && M.length != n*(n-1)) throw new RuntimeException("Migration matrix must have dimension stateNumber x (stateNumber-1) ("+n*(n-1)+")!");

        freq = frequencies.get().getValues();

        collectTimes(T);
        setRho();


    }


    void setupIntegrators(){   // set up ODE's and integrators

        if (minstep == null) minstep = tolerance.get();
        if (maxstep == null) maxstep = 1000.;

        P = new p0_ODE(birth,null, death,psi,M, n, totalIntervals, times);
        PG = new p0ge_ODE(birth, null, death,psi,M, n, totalIntervals, T, times, P, maxEvaluations.get(), true);

        pg_integrator = new DormandPrince853Integrator(minstep, maxstep, tolerance.get(), tolerance.get()); //
        pg_integrator.setMaxEvaluations(maxEvaluations.get());

        PG.p_integrator = new  DormandPrince853Integrator(minstep, maxstep, tolerance.get(), tolerance.get()); //
        PG.p_integrator.setMaxEvaluations(maxEvaluations.get());
    }


    double updateRates(){

        birth = new Double[n*totalIntervals];
        death = new Double[n*totalIntervals];
        psi = new Double[n*totalIntervals];
        if (SAModel) r =  new Double[n * totalIntervals];

        if (transform) {
            transformParameters();
        }
        else {

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

        M = migrationMatrix.get().getValues();

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

                if (Math.abs(time-tipTime) < 1e-10 && rho[((MultiTypeNode)tip).getNodeType()*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

            }
        }
    }

    public double[] getG(double t, double[] PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

        if (node.isLeaf()) {

            System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0, 0, n);
        }
//        else if (rhoSamplingChangeTimes.contains(t)){
//
//            int nodestate = ((MultiTypeNode)node).getNodeType();
//            PG0[nodestate] *= (1-rho[nodestate*totalIntervals+ Utils.index(t,times,totalIntervals)]);
//        }

        return getG(t,  PG0,  t0);
    }


    public double[] getG(double t, double[] PG0, double t0){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

        try {

            if (Math.abs(T-t)<1e-10 || Math.abs(t0-t)<1e-10 ||  T < t) {
                return PG0;
            }

            double from = t;
            double to = t0;
            double oneMinusRho;

            int indexFrom = Utils.index(from, times, times.length);
            int index = Utils.index(to, times, times.length);

            int steps = index - indexFrom;
            if (Math.abs(from-times[indexFrom])<1e-10) steps--;
            if (index>0 && Math.abs(to-times[index-1])<1e-10) {
                steps--;
                index--;
            }
            index--;

            while (steps > 0){

                from = times[index];// + 1e-14;

                pg_integrator.integrate(PG, to, PG0, from, PG0); // solve PG , store solution in PG0

                if (rhoChanges>0){
                    for (int i=0; i<n; i++){
                        oneMinusRho = (1-rho[i*totalIntervals + Utils.index(times[index], times, totalIntervals)]);
                        PG0[i] *= oneMinusRho;
                        PG0[i+n] *= oneMinusRho;
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


    void updateOrigin(Node root){

        T = origin.get().getValue();
        orig = T - root.getHeight();

        if (originIsRootEdge.get()) {

            orig = origin.get().getValue();
            T = orig + coltree.getRoot().getHeight();
        }

    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        if (SAModel && treeInput.isDirty()) throw new RuntimeException("Error: SA Model only implemented for fixed trees!");

        coltree = (MultiTypeTree) tree;

        MultiTypeNode root = (MultiTypeNode) coltree.getRoot();
//        int node_state =  ((MultiTypeNode) coltree.getRoot()).getNodeType(); //.getNodeColour(root);


        if (!coltree.isValid() || (origin.get()!=null && !originBranchIsValid(root))){
            logP =  Double.NEGATIVE_INFINITY;
            return logP;
        }

        int node_state;
        if (origin.get()==null) {
            T = root.getHeight();
            node_state =  ((MultiTypeNode) coltree.getRoot()).getNodeType(); //.getNodeColour(root);

        }
        else{
            updateOrigin(root);
            node_state = (originBranch.getChangeCount()>0) ? originBranch.getChangeType(originBranch.getChangeCount()-1) : ((MultiTypeNode) coltree.getRoot()).getNodeType();
            //originBranch.getFinalType();

            if (orig < 0){
                logP =  Double.NEGATIVE_INFINITY;
                return logP;
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

            if (conditionOnSurvival.get()) {
                noSampleExistsProp = PG.getP(0,m_rho.get()!= null,rho);
                if (print) System.out.println("\nnoSampleExistsProp = " + noSampleExistsProp[0]);// + ", " + noSampleExistsProp[1]);

                if ((noSampleExistsProp[node_state] < 0) || (noSampleExistsProp[node_state] > 1) || (Math.abs(1 - noSampleExistsProp[node_state]) < 1e-14)) {
                    logP = Double.NEGATIVE_INFINITY;
                    return logP;
                }
            }

            double[] p;

            if (orig>0){
                if (originBranch.getChangeCount()>0) {
                    p = calculateOriginLikelihood(originBranch.getChangeCount()-1, 0, T-originBranch.getChangeTime(originBranch.getChangeCount()-1) );
                }
                else
                    p = calculateSubtreeLikelihood(root, false, null, 0, orig);
            }
            else {
                int childIndex = 0;
                if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

                double t0 = T - root.getChild(childIndex).getHeight();
                int childChangeCount = ((MultiTypeNode)root.getChild(childIndex)).getChangeCount();
                if (childChangeCount > 0)
                    t0 = T - ((MultiTypeNode)root.getChild(childIndex)).getChangeTime(childChangeCount-1);

                p = calculateSubtreeLikelihood(root.getChild(childIndex), false, null, 0., t0);

                childIndex = Math.abs(childIndex-1);

                t0 = T - root.getChild(childIndex).getHeight();
                childChangeCount = ((MultiTypeNode)root.getChild(childIndex)).getChangeCount(); // changeCounts[root.getChild(1).getNr()];
                if (childChangeCount > 0)
                    t0 = T - ((MultiTypeNode)root.getChild(childIndex)).getChangeTime(childChangeCount-1);

                double[] p1 = calculateSubtreeLikelihood(root.getChild(childIndex), false, null, 0., t0);

                for (int i=0; i<p.length; i++) p[i]*=p1[i];

            }
            if (conditionOnSurvival.get()) p[n+node_state] /= (1-noSampleExistsProp[node_state]);    // condition on survival

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

    double[] calculateOriginLikelihood(Integer migIndex, double from, double to) {

        double[] init = new double[2*n];

        int prevcol = originBranch.getChangeType(migIndex);
        int col =  (migIndex > 0)?  originBranch.getChangeType(migIndex-1):  ((MultiTypeNode) coltree.getRoot()).getNodeType();

        migIndex--;

        double[] g ;

        if (migIndex >= 0){

            g = calculateOriginLikelihood(migIndex, to, T - originBranch.getChangeTime(migIndex));

            System.arraycopy(g, 0, init, 0, n);
            init[n+prevcol] = M[prevcol*(n-1)+(col<prevcol?col:col-1)] * g[n+col];

            return getG(from, init, to);

        }
        else {

            g = calculateSubtreeLikelihood(coltree.getRoot(), false, null, to, orig);

            System.arraycopy(g, 0, init, 0, n);
            init[n+prevcol] = M[prevcol*(n-1)+(col<prevcol?col:col-1)] * g[n+col];

            return getG(from, init, to, coltree.getRoot());
        }
    }

    double[] calculateSubtreeLikelihood(Node node, Boolean migration, Integer migIndex, double from, double to) {

        double[] init = new double[2*n];
        int nodestate = ((MultiTypeNode)node).getNodeType();
        int index = Utils.index(to, times, totalIntervals);

        if (migration){ // migration event

            int prevcol = ((MultiTypeNode) node).getChangeType(migIndex); //coltree.getChangeColour(node, migIndex);
            int col =  (migIndex > 0)?  ((MultiTypeNode) node).getChangeType(migIndex-1):  ((MultiTypeNode) node).getNodeType(); //  (migIndex > 0)? coltree.getChangeColour(node, migIndex-1): coltree.getNodeColour(node);
            double time ;

            migIndex--;

            time = (migIndex >= 0)? ((MultiTypeNode) node).getChangeTime(migIndex) :node.getHeight();// (migIndex >= 0)?coltree.getChangeTime(node, migIndex):node.getHeight();
            double[] g = calculateSubtreeLikelihood(node, (migIndex >= 0), migIndex, to, T-time);

            System.arraycopy(g, 0, init, 0, n);
            init[n+prevcol] = M[prevcol*(n-1)+(col<prevcol?col:col-1)] * g[n+col];

            return getG(from, init, to, node);
        }

        else {

            if (migIndex==null &&  ((MultiTypeNode)node).getChangeCount()>0){ // node has migration event(psi)

                return calculateSubtreeLikelihood(node, true, ((MultiTypeNode)node).getChangeCount()-1, from, to) ;
            }

            else{

                if (node.isLeaf()){ // sampling event

                    if (!isRhoTip[node.getNr()])
//                        init[n+nodestate] = psi[nodestate*totalIntervals+index];
                        init[n + nodestate] = SAModel
                                ? psi[nodestate * totalIntervals + index]* (r[nodestate * totalIntervals + index] + (1-r[nodestate * totalIntervals + index])*PG.getP(to, m_rho.get()!=null, rho)[nodestate]) // with SA: ψ_i(r + (1 − r)p_i(τ))
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
                    childChangeCount = ((MultiTypeNode)node.getChild(childIndex)).getChangeCount(); // changeCounts[node.getChild(1).getNr()];
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




    public void transformParameters(){

        Double[] p = samplingProportion.get().getValues();
        Double[] ds = becomeUninfectiousRate.get().getValues();
        Double[] R = R0.get().getValues();
        Double[] removalProbabilities = new Double[1];
        if (SAModel) removalProbabilities = removalProbability.get().getValues();

        int state;

        for (int i = 0; i < totalIntervals*n; i++){

            state =  i/totalIntervals;

            birth[i] = R[R.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state]
                    * ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state];

//            psi[i] = p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]
//                    * ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] ;
//
//            death[i] = ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] - psi[i];
//
            if (!SAModel || removalAffectsSamplingProportion.get())
                psi[i] = p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]
                        * ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] ;

            if (!SAModel)
                death[i] = ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] - psi[i];

            else {
                r[i] = removalProbabilities[removalProbabilities.length > n ? (rChanges+1)*state+index(times[i%totalIntervals], rChangeTimes) : state];

                if (!removalAffectsSamplingProportion.get())
                    psi[i] = p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]
                            * ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state]
                            / (1+(r[i]-1)*p[p.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state]);


                death[i] = ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] - psi[i]*r[i];
            }
        }

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
