package bdmm.util;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.XMLProducer;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author: denise.kuehnert@gmail.com adapted from SimulatedAlignment by remco@cs.waikato.ac.nz
 * @date:   07.03.17
 */
@Description("A NUCLEOTIDE alignment containing sequences randomly generated using a"
        + "given site model down a given tree with automatic taxon creation.")
public class SimulatedAlignmentAutoTaxa extends Alignment {
    final public Input<Alignment> m_data = new Input<>("data", "alignment data which specifies datatype and taxa of the beast.tree", Input.Validate.REQUIRED);
    final public Input<Tree> m_treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Input.Validate.REQUIRED);
    final public Input<SiteModel.Base> m_pSiteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Input.Validate.REQUIRED);
    final public Input<BranchRateModel.Base> m_pBranchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<Integer> m_sequenceLengthInput = new Input<>("sequencelength", "nr of samples to generate (default 1000).", 1000);
    final public Input<String> m_outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");

//    DataType dataType;
    /**
     * nr of samples to generate *
     */
    protected int m_sequenceLength;
    /**
     * tree used for generating samples *
     */
    protected Tree m_tree;
    /**
     * site model used for generating samples *
     */
    protected SiteModel.Base m_siteModel;
    /**
     * branch rate model used for generating samples *
     */
    protected BranchRateModel m_branchRateModel;
    /**
     * nr of categories in site model *
     */
    int m_categoryCount;
    /**
     * nr of states in site model *
     */
    int m_stateCount;

    /**
     * name of output file *
     */
    String m_outputFileName;

    /**
     * an array used to transfer transition probabilities
     */
    protected double[][] m_probabilities;

    public SimulatedAlignmentAutoTaxa() {

        // Override the sequence input requirement.
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        m_tree = m_treeInput.get();
        m_siteModel = m_pSiteModelInput.get();
        m_branchRateModel = m_pBranchRateModelInput.get();
        m_sequenceLength = m_sequenceLengthInput.get();
        m_stateCount = m_data.get().getMaxStateCount();
        m_categoryCount = m_siteModel.getCategoryCount();
        m_probabilities = new double[m_categoryCount][m_stateCount * m_stateCount];
        m_outputFileName = m_outputFileNameInput.get();

        sequenceInput.get().clear();

        m_dataType = m_data.get().getDataType();

        int taxonCount = m_treeInput.get().getLeafNodeCount();
        List<Sequence> sequenceList = new ArrayList<Sequence>(taxonCount);

        Sequence s;
        for (int i=0; i<taxonCount; i++) {
            s = new Sequence();
            s.initByName("totalcount","4",
                            "taxon", Integer.toString(i+1),
                            "value","?");
            sequenceList.add(s);
        }

        initializeWithSequenceList(sequenceList,true);

//        for (int i=2; i<taxonCount; i++) {
//            m_data.get().setInputValue("sequence", sequenceList.get(i));
//        }
//        m_data.get().sequenceInput.get().clear();
        m_data.get().initByName("sequenceInput", sequenceList);

        simulate();

        // Write simulated alignment to disk if requested:
        if (m_outputFileName != null) {
            PrintStream pstream;
            try {
                pstream = new PrintStream(m_outputFileName);
                pstream.println(new XMLProducer().toRawXML(this));
                pstream.close();
            } catch (FileNotFoundException e) {
                throw new IllegalArgumentException(e.getMessage());
            }
        }
        if (m_outputFileName != null) {
            PrintStream pstream;
            try {
                pstream = new PrintStream(m_outputFileName+".fasta");
                for (Sequence sequence : sequenceInput.get()){
                        pstream.println(">" + sequence.toString().replace(":","\n"));
                }
                pstream.close();
            } catch (FileNotFoundException e) {
                throw new IllegalArgumentException(e.getMessage());
            }
        }

        super.initAndValidate();
    }

    /**
     * Initializes the alignment given the provided list of sequences and no other information.
     * It site weights and/or data type have been previously set up with initAndValidate then they
     * remain in place. This method is used mainly to re-order the sequences to a new taxon order
     * when an analysis of multiple alignments on the same taxa are undertaken.
     *
     * @param sequences
     */
    private void initializeWithSequenceList(List<Sequence> sequences, boolean log) {
        this.sequences = sequences;
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        try {
            for (Sequence seq : sequences) {

                counts.add(seq.getSequence(m_dataType));
                if (taxaNames.contains(seq.getTaxon())) {
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon());
                }
                taxaNames.add(seq.getTaxon());
                tipLikelihoods.add(seq.getLikelihoods());
                // if seq.isUncertain() == false then the above line adds 'null'
                // to the list, indicating that this particular sequence has no tip likelihood information

                if (seq.totalCountInput.get() != null) {
                    stateCounts.add(seq.totalCountInput.get());
                } else {
                    stateCounts.add(m_dataType.getStateCount());
                }
            }
            if (counts.size() == 0) {
                // no sequence data
                throw new RuntimeException("Sequence data expected, but none found");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    /**
     * Convert integer representation of sequence into a Sequence
     *
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     */
    Sequence intArray2Sequence(int[] seq, Node node) {
        String seqString = m_dataType.state2string(seq);

        String taxon = getTaxaNames().get(node.getNr());
        return new Sequence(taxon, seqString);
    } // intArray2Sequence

    /**
     * perform the actual sequence generation
     *
     * @return alignment containing randomly generated sequences for the nodes in the
     *         leaves of the tree
     */
    public void simulate() {
        Node root = m_tree.getRoot();


        double[] categoryProbs = m_siteModel.getCategoryProportions(root);
        int[] category = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            category[i] = Randomizer.randomChoicePDF(categoryProbs);
        }

        double[] frequencies = m_siteModel.getSubstitutionModel().getFrequencies();
        int[] seq = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            seq[i] = Randomizer.randomChoicePDF(frequencies);
        }

        traverse(root, seq, category);

    } // simulate

    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     */
    void traverse(Node node, int[] parentSequence, int[] category) {
        for (int childIndex = 0; childIndex < 2; childIndex++) {
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());
            for (int i = 0; i < m_categoryCount; i++) {
                getTransitionProbabilities(m_tree, child, i, m_probabilities[i]);
            }

            int[] seq = new int[m_sequenceLength];
            double[] cProb = new double[m_stateCount];
            for (int i = 0; i < m_sequenceLength; i++) {
                System.arraycopy(m_probabilities[category[i]], parentSequence[i] * m_stateCount, cProb, 0, m_stateCount);
                seq[i] = Randomizer.randomChoicePDF(cProb);
            }

            if (child.isLeaf()) {
                sequenceInput.setValue(intArray2Sequence(seq, child), this);
            } else {
                traverse(child, seq, category);
            }
        }
    } // traverse

    /**
     * get transition probability matrix for particular rate category *
     */
    void getTransitionProbabilities(Tree tree, Node node, int rateCategory, double[] probs) {

        Node parent = node.getParent();
        double branchRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(node));
        branchRate *= m_siteModel.getRateForCategory(rateCategory, node);

        //m_siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probs);
        m_siteModel.getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);

    } // getTransitionProbabilities


} // class SequenceAlignment

