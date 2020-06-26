package bdmmprime.trajectories.obsevents;

public class ObservedTypeChangeEvent extends ObservedEvent {

    public int parentType, childType;

    public ObservedTypeChangeEvent(double time, int parentType, int childType, int multiplicity) {
        this.time = time;
        this.parentType = parentType;
        this.childType = childType;
        this.multiplicity = multiplicity;
    }
}
