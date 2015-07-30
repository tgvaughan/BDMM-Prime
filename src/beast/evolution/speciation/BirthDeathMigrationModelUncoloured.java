package beast.evolution.speciation;

import beast.evolution.tree.*;
import beast.core.parameter.RealParameter;
import beast.core.Input;
import beast.core.Description;
import beast.core.util.Utils;

import math.p0_ODE;
import math.p0ge_ODE;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;


/**
 * @author Denise Kuehnert
 * Date: Jul 2, 2013
 * Time: 10:28:16 AM
 *
 */

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
        "This should be used when the migration process along the phylogeny is irrelevant. Otherwise the BirthDeathMigrationModel can be employed." +
        "This implementation also works with sampled ancestor trees.")
public class BirthDeathMigrationModelUncoloured extends PiecewiseBirthDeathSamplingDistribution {


    public Input<RealParameter> migrationMatrix =
            new Input<>("migrationMatrix", "Flattened migration matrix, can be asymmetric, diagnonal entries omitted",  Input.Validate.REQUIRED);

    public Input<RealParameter> frequencies =
            new Input<>("frequencies", "state frequencies",  Input.Validate.REQUIRED);

    public Input<RealParameter> origin =
            new Input<>("origin", "The origin of infection x1");

    public Input<Boolean> originIsRootEdge =
            new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);

    public Input<Integer> maxEvaluations =
            new Input<>("maxEvaluations", "The maximum number of evaluations for ODE solver", 20000);

    public Input<Boolean> conditionOnSurvival =
            new Input<>("conditionOnSurvival", "condition on at least one survival? Default true.", true);

    public Input<Double> tolerance =
            new Input<>("tolerance", "tolerance for numerical integration", 1e-14);

    public Input<TraitSet> tiptypes = new Input<>("tiptypes", "trait information for initializing traits (like node types/locations) in the tree",  Input.Validate.REQUIRED);
    public Input<String> typeLabel = new Input<>("typeLabel", "type label in tree for initializing traits (like node types/locations) in the tree",  Input.Validate.XOR, tiptypes);

    public Input<Boolean> storeNodeTypes = new Input<>("storeNodeTypes", "store tip node types? this assumes that tip types cannot change (default false)", false);
    public Input<Boolean> checkRho = new Input<>("checkRho", "check if rho is set if multiple tips are given at present (default true)", true);

    Double[] freq;
    Double[] M;
    double T;
    double orig;
    int ntaxa;

    p0_ODE P;
    p0ge_ODE PG;

    FirstOrderIntegrator pg_integrator;
    public int maxEvalsUsed;
    public Double minstep;
    public Double maxstep;

    private int[] nodeStates;

    Boolean print = false;

    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();

        TreeInterface tree = treeInput.get();

        if (origin.get()==null){
            T = tree.getRoot().getHeight();
        }

        else{

            T = origin.get().getValue();
            orig = T -  tree.getRoot().getHeight();


            if (originIsRootEdge.get()) {

                orig = origin.get().getValue();
                T = orig + tree.getRoot().getHeight();
            }

            if (!Boolean.valueOf(System.getProperty("beast.resume")) && orig < 0)
                throw new RuntimeException("Error: origin("+T+") must be larger than tree height("+tree.getRoot().getHeight()+")!");
        }

        ntaxa = tree.getLeafNodeCount();

        birthAmongDemes = (birthRateAmongDemes.get() !=null || R0AmongDemes.get()!=null);

        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null){

            if (birthAmongDemes) b_ij = birthRateAmongDemes.get().getValues();

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

        int contempCount = 0;
        for (Node node : tree.getExternalNodes())
            if (node.getHeight()==0.)
                contempCount++;
        if (checkRho.get() && contempCount>1 && rho==null)
            throw new RuntimeException("Error: multiple tips given at present, but sampling probability \'rho\' is not specified.");

        M = migrationMatrix.get().getValues();

        if (n>1 && M.length != n*(n-1)) throw new RuntimeException("Migration matrix must have dimension stateNumber x (stateNumber-1) ("+n*(n-1)+")!");

        freq = frequencies.get().getValues();

        double freqSum = 0;
        for (double f : freq) freqSum+= f;
        if (freqSum!=1.) throw new RuntimeException("Error: frequencies must add up to 1.");

//        // calculate equilibrium frequencies for 2 types:
//        double LambMu = -b[0]-b[1]-(d[0]+s[0])+(d[1]+s[1]);
//        double c = Math.sqrt(Math.pow(LambMu,2) +4*b_ij[0]*b_ij[1]);
//        freq[0] = (c+LambMu)/(c+LambMu+2*b_ij[0]) ;
//        freq[1] = 1 - freq[0];


        collectTimes(T);
        setRho();

        maxEvalsUsed = 0;

        if (storeNodeTypes.get()) {

            nodeStates = new int[ntaxa];

            for (Node node : tree.getExternalNodes()){
                nodeStates[node.getNr()] = getNodeState(node, true);
            }
        }

    }


    void setupIntegrators(){   // set up ODE's and integrators

        if (minstep == null) minstep = tolerance.get();
        if (maxstep == null) maxstep = 1000.;

        P = new p0_ODE(birth, (birthAmongDemes ? b_ij : null), death,psi,M, n, totalIntervals, times);
        PG = new p0ge_ODE(birth, (birthAmongDemes ? b_ij : null), death,psi,M, n, totalIntervals, T, times, P, maxEvaluations.get(), false);

        pg_integrator = new DormandPrince853Integrator(minstep, maxstep, tolerance.get(), tolerance.get()); //
        pg_integrator.setMaxEvaluations(maxEvaluations.get());

        PG.p_integrator = new  DormandPrince853Integrator(minstep, maxstep, tolerance.get(), tolerance.get()); //
        PG.p_integrator.setMaxEvaluations(maxEvaluations.get());
    }

    protected Double updateRates(TreeInterface tree) {

        if (origin.get()==null){
            T = tree.getRoot().getHeight();
        }

        else{

            T = origin.get().getValue();
            orig = T -  tree.getRoot().getHeight();


            if (originIsRootEdge.get()) {

                orig = origin.get().getValue();
                T = orig + tree.getRoot().getHeight();
            }

            if (!Boolean.valueOf(System.getProperty("beast.resume")) && orig < 0)
                throw new RuntimeException("Error: origin("+T+") must be larger than tree height("+tree.getRoot().getHeight()+")!");
        }

        birth = new Double[n*totalIntervals];
        death = new Double[n*totalIntervals];
        psi = new Double[n*totalIntervals];
        b_ij = new Double[totalIntervals*(n*(n-1))];
        if (SAModel) r =  new Double[n * totalIntervals];

        if (transform) {
            transformParameters();
        }
        else {

            Double[] birthRates = birthRate.get().getValues();
            Double[] deathRates = deathRate.get().getValues();
            Double[] samplingRates = samplingRate.get().getValues();
            Double[] birthAmongDemesRates = new Double[1];
            if (birthAmongDemes) birthAmongDemesRates = birthRateAmongDemes.get().getValues();
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

            if (birthAmongDemes)    {

                for (int i = 0; i < n; i++){
                    for (int j=0; j<n ; j++){
                        for (int dt=0; dt<totalIntervals; dt++){
                            if (i!=j){
                                b_ij[(i*(n-1)+(j<i?j:j-1))*totalIntervals+dt]
                                        = birthAmongDemesRates[(birthAmongDemesRates.length>(n*(n-1)))
                                        ?  (b_ij_Changes+1)*(n-1)*i + index(times[dt], b_ijChangeTimes)
                                        : (i*(n-1)+(j<i?j:j-1))];
                            }
                        }
                    }
                }
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

                if (Math.abs(time-tipTime) < 1e-10 && rho[getNodeState(tip,false)*totalIntervals + Utils.index(time, times, totalIntervals)]>0) isRhoTip[tip.getNr()] = true;

            }
        }
    }

    public double[] getG(double t, double[] PG0, double t0, Node node){ // PG0 contains initial condition for p0 (0..n-1) and for ge (n..2n-1)

        try {

            if (node.isLeaf()) {

                System.arraycopy(PG.getP(t0, m_rho.get()!=null, rho), 0, PG0, 0, n);
            }

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


    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        Node root = tree.getRoot();

        if (origin.get()==null){
            T =tree.getRoot().getHeight();
        }
        else{

            T = origin.get().getValue();
            orig = T - root.getHeight();


            if (originIsRootEdge.get()) {

                orig = origin.get().getValue();
                T = orig + tree.getRoot().getHeight();
            }

            if (orig < 0)
                return Double.NEGATIVE_INFINITY;
        }

        collectTimes(T);
        setRho();

        if (updateRates(tree) < 0 ||  (times[totalIntervals-1] > T)) {
            logP =  Double.NEGATIVE_INFINITY;
            return logP;
        }

        double[] noSampleExistsProp ;

        double Pr = 0;
        double nosample = 0;

        try{  // start calculation

            if (conditionOnSurvival.get()) {
                noSampleExistsProp = PG.getP(0,m_rho.get()!=null,rho);
                if (print) System.out.println("\nnoSampleExistsProp = " + noSampleExistsProp[0] + ", " + noSampleExistsProp[1]);


                for (int root_state=0; root_state<n; root_state++){
                    nosample += freq[root_state] *  noSampleExistsProp[root_state] ;
                }

                if (nosample<0 || nosample>1)
                    return Double.NEGATIVE_INFINITY;
            }

            double[] p;
            if ( orig > 0 ) {
                p = calculateSubtreeLikelihood(root,0,orig);
            }
            else {
                int childIndex = 0;
                if (root.getChild(1).getNr() > root.getChild(0).getNr()) childIndex = 1; // always start with the same child to avoid numerical differences

                p = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

                childIndex = Math.abs(childIndex-1);

                double [] p1 = calculateSubtreeLikelihood(root.getChild(childIndex),0., T - root.getChild(childIndex).getHeight());

                for (int i =0; i<p.length; i++) p[i]*=p1[i];
            }

            if (print) System.out.print("final p per state = ");

            for (int root_state=0; root_state<n; root_state++){

                if (p[n+root_state]>0 ) Pr += freq[root_state]* p[n+root_state];

                if (print) System.out.print(p[root_state] + "\t" + p[root_state+n] + "\t");
            }

            if (conditionOnSurvival.get()){
                Pr /= (1-nosample);
            }

        }catch(Exception e){

            if (e instanceof RuntimeException){throw e;}

            logP =  Double.NEGATIVE_INFINITY;
            return logP;
        }

        maxEvalsUsed = Math.max(maxEvalsUsed, PG.maxEvalsUsed);

        logP = Math.log(Pr);
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
            throw new RuntimeException("Something went wrong with the assignment of types to the nodes (node ID="+node.getID()+"). Please check your XML file!");
        }
    }


    double[] calculateSubtreeLikelihood(Node node, double from, double to) {

        double[] init = new double[2*n];

        int index = Utils.index(to,times, totalIntervals);

        if (node.isLeaf()){ // sampling event

            int nodestate = getNodeState(node, false);

            if (nodestate==-1) { //unknown state

                if (SAModel) throw new RuntimeException("SA model not implemented with unknown states!");

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
//                    init[n+nodestate] = psi[nodestate*totalIntervals+index];
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

//                            if (b_ij.length>(n*(n-1))){ // b_ij can change over time
//                                throw new RuntimeException("ratechanges in b_ij not implemented!");
                                init[n + childstate] += 0.5 * b_ij[totalIntervals * (childstate * (n - 1) + (j < childstate ? j : j - 1)) + index] * (g0[n + childstate] * g1[n + j] + g0[n + j] * g1[n + childstate]);
//                            } else { // b_ij cannot change over time
//                                init[n+childstate] += 0.5 * b_ij[childstate*(n-1)+(j<childstate?j:j-1)] * (g0[n+childstate] * g1[n+j] + g0[n+j] * g1[n+childstate]);
//                            }
                            }
                        }
                    }


                    if (Double.isInfinite(init[childstate])) {
                        throw new RuntimeException("infinite likelihood");
                    }
                }
            }
//            init[0] = g0[0];
//            init[1] = g0[1];
//            init[2] = b[0]*g0[2]*g1[2] + .5*b_ij[0]*(g0[2]*g1[3] + g0[3]*g1[2]);
//            init[3] = b[1]*g0[3]*g1[3] + .5*b_ij[1]*(g0[2]*g1[3] + g0[3]*g1[2]);
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

    public void transformParameters(){

        Double[] p = samplingProportion.get().getValues();
        Double[] ds = becomeUninfectiousRate.get().getValues();
        Double[] R = R0.get().getValues();
        Double[] RaD = (birthAmongDemes) ? R0AmongDemes.get().getValues() : new Double[1];
        Double[] removalProbabilities = new Double[1];
        if (SAModel) removalProbabilities = removalProbability.get().getValues();

//        birth = new Double[n * totalIntervals];
//        death = new Double[n * totalIntervals];
//        psi = new Double[n * totalIntervals];
//        b_ij = new Double[totalIntervals * (n * (n - 1))];

        if (coupledR0Changes.get()!=null){

            Double[] changes = coupledR0Changes.get().getValues();
            int c = changes.length / n;
            R = new Double[(c+1)*n];
            RaD = new Double[(c+1)*n*(n-1)];



            for (int i=0; i<n; i++){

                R[i*(c+1)] =  R0.get().getValue(i);

                for (int k=0; k<n-1; k++){

                    RaD[(i*(n-1)+k)*(c+1)] =  R0AmongDemes.get().getValue((i*(n-1)+k));
                }

                for (int j=1; j<=c; j++){

                    R[i*(c+1)+j] = R[i*(c+1)] * changes[i*c+j-1];
//                    RaD[j] = RaD[0] * changes[j-1];

                    for (int k=0; k<n-1; k++){

                        RaD[(i*(n-1)+k)*(c+1)+j] =  RaD[(i*(n-1)+k)*(c+1)] * changes[i*(n-1)*c+j-1];
                    }

                }
            }
        }

        int state;

        for (int i = 0; i < totalIntervals*n; i++){

            state =  i/totalIntervals;

//            if (birthAmongDemes)    {
//                 for (int j=1; j<n && state!=j ; j++){
//                        b_ij[((state*(n-1)+(j<state?j:j-1))-1)*totalIntervals+i%totalIntervals]
//                                = RaD[(RaD.length>(n*(n-1))) ? (b_ij_Changes+1)*n*(n-1)*state + index(times[i%totalIntervals], b_ijChangeTimes) : (state*(n-1)+(j<state?j:j-1))]
//                                    * ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] ;
//                 }
//            }

            birth[i] = R[R.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state]
                    * ds[ds.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state] ;

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

        if (birthAmongDemes)    {

            for (int i = 0; i < n; i++){

                for (int j=0; j<n ; j++){

                    for (int dt=0; dt<totalIntervals; dt++){

                        if (i!=j){
                            b_ij[(i*(n-1)+(j<i?j:j-1))*totalIntervals+dt]
                                    = RaD[(RaD.length>(n*(n-1)))
                                    ?  (b_ij_Changes+1)*(n-1)*i + index(times[dt], b_ijChangeTimes)
                                    : (i*(n-1)+(j<i?j:j-1))]
                                    * ds[ds.length > n ? (deathChanges+1)*i+index(times[dt], deathRateChangeTimes) : i];
                        }
                    }
                }

            }
        }
    }

}

