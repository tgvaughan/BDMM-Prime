package bdmmprime.trajectories.obsevents;

public class CoalescenceEvent extends ObservedEvent {

    public int parentType, childType1, childType2;

    public CoalescenceEvent(double time, int parentType, int childType1, int childType2, int multiplicity) {
        this.time = time;
        this.parentType = parentType;
        this.childType1 = childType1;
        this.childType2 = childType2;
        this.multiplicity = multiplicity;
    }
}
