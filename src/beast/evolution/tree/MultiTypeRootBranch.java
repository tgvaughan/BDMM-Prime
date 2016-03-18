package beast.evolution.tree;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * User: Denise
 * Date: 08.07.14
 * Time: 14:44
 */
@Description("The branch that would be above the root, if BEAST trees allowed that. " +
        "It is needed for BirthDeathMigrationModel, to allow type change events along the root branch.")
public class MultiTypeRootBranch extends StateNode {

    // Type metadata:
    int nTypeChanges = 0;
    int stored_nTypeChanges = 0;
    List<Integer> changeTypes = new ArrayList<Integer>();
    List<Double> changeTimes = new ArrayList<Double>();

    List<Integer> stored_changeTypes = new ArrayList<Integer>();
    List<Double> stored_changeTimes = new ArrayList<Double>();

    @Override
    public void initAndValidate() {}

    /**
     * Retrieve the total number of changes on the branch above this node.
     *
     * @return type change count
     */
    public int getChangeCount() {
        return nTypeChanges;
    }

    /**
     * Obtain destination type (reverse time) of the change specified by idx on
     * the branch between this node and its parent.
     *
     * @param idx
     * @return change type
     */
    public int getChangeType(int idx) {
        return changeTypes.get(idx);
    }

    /**
     * Obtain time of the change specified by idx on the branch between this
     * node and its parent.
     *
     * @param idx
     * @return time of change
     */
    public double getChangeTime(int idx) {
        return changeTimes.get(idx);
    }

    /**
     * Add a new type change to branch above node.
     *
     * @param newType Destination (reverse time) type of change
     * @param time Time at which change occurs.
     */
    public void addChange(int newType, double time) {
        startEditing(null);
        changeTypes.add(newType);
        changeTimes.add(time);
        nTypeChanges += 1;
    }

    /**
     * Remove all type changes from branch above node.
     */
    public void clearChanges() {
        startEditing(null);
        changeTypes.clear();
        changeTimes.clear();
        nTypeChanges = 0;
    }

    @Override
    public MultiTypeRootBranch copy() {

        MultiTypeRootBranch node = new MultiTypeRootBranch();

        node.ID = ID;
        node.nTypeChanges = nTypeChanges;

        node.changeTimes.addAll(changeTimes);
        node.changeTypes.addAll(changeTypes);

        return node;
    }


    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
     protected boolean requiresRecalculation() {
         return true;
     }

    @Override
    public void assignTo(StateNode other) {

        MultiTypeRootBranch otherBranch = (MultiTypeRootBranch) other;

        otherBranch.ID = getID();
        otherBranch.nTypeChanges = nTypeChanges;

        otherBranch.changeTimes.clear();
        otherBranch.changeTimes.addAll(changeTimes);

        otherBranch.changeTypes.clear();
        otherBranch.changeTypes.addAll(changeTypes);
    }

    @Override
    public void assignFrom(StateNode other) {

        MultiTypeRootBranch otherBranch = (MultiTypeRootBranch) other;

        ID = otherBranch.getID();
        nTypeChanges = otherBranch.nTypeChanges;

        changeTimes.clear();
        changeTimes.addAll(otherBranch.changeTimes);

        changeTypes.clear();
        changeTypes.addAll(otherBranch.changeTypes);
    }

    @Override
    public void assignFromFragile(StateNode other) {
        assignFrom(other);
    }

    @Override
     public String toString() { // format: nTypeChanges:(type1,time1):(type2,time2) etc.

        String string = nTypeChanges + "";

        for (int iValue = 0; iValue < nTypeChanges; iValue++) {
            string += ":(" + getChangeType(iValue) + "," + getChangeTime(iValue) + ")";
        }
        return string;
     }


    @Override
    public void fromXML(Node node) {

        String[] s = node.getTextContent().split(":");
        String[] change;

        nTypeChanges = Integer.parseInt(s[0]);

        for (int i=0; i<nTypeChanges; i++){

            change = s[i+1].substring(1,s[i+1].length()-1).split(",");

            changeTypes.add(Integer.parseInt(change[0]));
            changeTimes.add(Double.parseDouble(change[1]));
        }
    }

    @Override
    public int scale(double fScale) {
        throw new RuntimeException("MultiTypeRootBranch cannot be scaled.");
    }

    @Override
    protected void store() {

        stored_nTypeChanges = nTypeChanges;

        stored_changeTimes.clear();
        stored_changeTimes.addAll(changeTimes);

        stored_changeTypes.clear();
        stored_changeTypes.addAll(changeTypes);

        hasStartedEditing = true;
    }

    @Override
    public void restore() {

        nTypeChanges = stored_nTypeChanges;

        changeTimes.clear();
        changeTimes.addAll(stored_changeTimes);

        changeTypes.clear();
        changeTypes.addAll(stored_changeTypes);

        hasStartedEditing = false;
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        return nTypeChanges;
    }

    @Override
    public double getArrayValue(int iDim) {
        return nTypeChanges;
    }

    @Override
    public void init(PrintStream out) {
        out.print(getID() + "\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(getChangeCount()  + "\t"); // just log number of type changes
//        out.print(toString()  + "\t");    // this would be ideal to store in an extra log file to report details of type changes
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

}
