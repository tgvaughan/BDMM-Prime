package bdmm.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.MultiTypeNode;
import bdmm.tree.MultiTypeRootBranch;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import multitypetree.operators.RandomRetypeOperator;

import java.util.ArrayList;
import java.util.List;

/**
 * User: Denise
 * Date: 08.07.14
 * Time: 17:57
 */
@Description("Retypes a randomly chosen node and its attached branches. "
        + "This variant uses the uniformization branch retyping procedure.")
public class RootBranchRetype extends RandomRetypeOperator {

    public Input<MultiTypeRootBranch> multiTypeRootBranchInput =
            new Input<MultiTypeRootBranch>("multiTypeRootBranch", "MultiTypeRootBranch for origin coloring", Input.Validate.REQUIRED);

     public Input<RealParameter> originInput =
             new Input<RealParameter>("origin", "The origin of infection x1", Input.Validate.REQUIRED);

    MultiTypeRootBranch rootBranch;

    @Override
    public void initAndValidate() {

        super.initAndValidate();
    }

    @Override
    public double proposal() {

        rootBranch = multiTypeRootBranchInput.get();
        mtTree = multiTypeTreeInput.get();

        Node rootNode = mtTree.getRoot();

        double logHR = 0.0;

        // Select new node type:
        int newType = Randomizer.nextInt(migModel.getNTypes());

        logHR += getBranchTypeProb(rootBranch);
        logHR += getBranchTypeProb(rootNode.getLeft()) + getBranchTypeProb(rootNode.getRight());

        // Set new root type
        ((MultiTypeNode) rootNode).setNodeType(newType);

        //  Retype attached branches:
        logHR -= retypeBranch(rootBranch);

        logHR -= retypeBranch(rootNode.getLeft()) + retypeBranch(rootNode.getRight());

        if ((((MultiTypeNode)rootNode.getLeft()).getFinalType() != ((MultiTypeNode) rootNode).getNodeType())
                || (((MultiTypeNode)rootNode.getRight()).getFinalType() != ((MultiTypeNode) rootNode).getNodeType()))
            return Double.NEGATIVE_INFINITY;

        return logHR;
    }



    /**
     * Retype branch between srcNode and its parent with rate fixed by the
     * tuning parameter mu. (adapted from Tim's RandomRetypeOperator)
     *
     * @param srcNode
     * @return Probability of branch typing
     */
    protected double retypeBranch(MultiTypeRootBranch srcNode) {

        double mu = muInput.get();

        double t_srcNode = mtTree.getRoot().getHeight();
        double t_srcNodeParent = originInput.get().getValue();

        int srcNodeType = ((MultiTypeNode) mtTree.getRoot()).getNodeType();

        // Clear existing changes in preparation for adding replacements:
        srcNode.clearChanges();

        double t = 0.;
        double length =  t_srcNodeParent - t_srcNode;
        List<Double> times = new ArrayList<>();

        // generate forward times of type changes:
        while (t < length) {

            // Determine time to next type change event:
            t += Randomizer.nextExponential(mu);

            if (t < length) times.add(t);

        }

        // select change types backwards:
        int lastType = srcNodeType;
        for (int i=times.size()-1; i>=0; i--){

                int newType = Randomizer.nextInt(migModel.getNTypes() - 1);
                if (newType >= lastType)
                    newType += 1;
                srcNode.addChange(newType, (t_srcNodeParent-times.get(i)));

                lastType = newType;

        }
        // Return log of branch type probability:
        return -mu*(t_srcNodeParent - t_srcNode)
                + srcNode.getChangeCount()*Math.log(mu/(migModel.getNTypes() - 1));

    }


    /**
     * Get probability of the colouring along the branch between srcNode
     * and its parent. (adapted from Tim's RandomRetypeOperator)
     *
     * @param srcNode
     * @return Probability of the colouring.
     */
    protected double getBranchTypeProb(MultiTypeRootBranch srcNode) {

        double mu = muInput.get();
        int n = srcNode.getChangeCount();
        int N = migModel.getNTypes();

        double T = originInput.get().getValue() - mtTree.getRoot().getHeight();

        if (N == 0)
            return 0.0;
        else
            return -mu*T + n*Math.log(mu/(N-1));
    }

}
