package bdmmprime.util;

import beast.base.core.*;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("A Function representing a number of elements of another Function.")
public class Slice extends CalculationNode implements RealVector<Real>, Loggable {

    public Input<RealVector<?>> realVectorInput = new Input<>("arg",
            "Argument to extract element from.", Input.Validate.REQUIRED);

    public Input<Integer> startIndexInput = new Input<>("index",
            "Index of first element to extract.", Input.Validate.REQUIRED);

    public Input<Integer> countInput = new Input<>("count",
            "Number of elements to extract.", 1);

    protected int indexStart, indexEnd, count;

    @Override
    public void initAndValidate() {
        indexStart = startIndexInput.get();
        count = countInput.get();
        indexEnd = indexStart + count - 1;

        if (indexEnd >= realVectorInput.get().size())
            throw new IllegalArgumentException("Index and count arguments to" +
                    " Slice are out of bounds.");
    }

    @Override
    public List<Double> getElements() {
        List<Double> elements = new ArrayList<>();
        for (int i=0; i<count; i++)
            elements.add(get(i));

        return elements;
    }

    @Override
    public Real getDomain() {
        return realVectorInput.get().getDomain();
    }

    @Override
    public int size() {
        return count;
    }

    @Override
    public double get(int iDim) {
        if (iDim < count)
            return realVectorInput.get().get(indexStart + iDim);
        else
            return 0;
    }

    @Override
    public void init(PrintStream out) {
        for (int i=0; i<count; i++)
            out.print(((BEASTObject) realVectorInput.get()).getID()
                    + "[" + (indexStart + i) + "]\t");
    }

    @Override
    public void log(long nSample, PrintStream out) {
        for (int i=0; i<count; i++)
            out.print(get(i) + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }
}
