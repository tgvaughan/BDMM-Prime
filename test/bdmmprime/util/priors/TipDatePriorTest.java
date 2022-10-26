package bdmmprime.util.priors;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import junit.framework.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class TipDatePriorTest {


    @Test
    public void test() {

        List<Taxon> taxonList = new ArrayList<>();
        taxonList.add(new Taxon("t1"));
        taxonList.add(new Taxon("t2"));
        taxonList.add(new Taxon("t3"));
        taxonList.add(new Taxon("t4"));
        taxonList.add(new Taxon("t5"));
        TaxonSet taxonSet = new TaxonSet(taxonList);


        TraitSet tipDates = new TraitSet();
        tipDates.initByName(
                "traitname", "date-forward",
                "taxa", taxonSet,
                "dateFormat", "yyyy-M-dd",
                "value", "t1=2020-05-15," +
                        "t2=2020-02-15," +
                        "t3=2021-05-15," +
                        "t4=2020-05-15," +
                        "t5=2020-05-15");

        TraitSet earlierBounds = new TraitSet();
        earlierBounds.initByName(
                "traitname", "date-forward",
                "taxa", taxonSet,
                "dateFormat", "yyyy-M-dd",
                "value", "t1=2020-05-01," +
                        "t2=2020-02-01," +
                        "t3=2021-05-01," +
                        "t4=2020-05-14," +
                        "t5=2020-05-14");

        TraitSet laterBounds = new TraitSet();
        laterBounds.initByName(
                "traitname", "date-forward",
                "taxa", taxonSet,
                "dateFormat", "yyyy-M-dd",
                "value", "t1=2020-06-01," +
                        "t2=2020-03-01," +
                        "t3=2021-06-01," +
                        "t4=2020-05-15," +
                        "t5=2020-05-15");

        Tree tree = new Tree();
        tree.initByName("taxonset", taxonSet,
                "trait", tipDates);

        System.out.println(tree);

        double fso = laterBounds.getDate(0) - tipDates.getDate(0);
        RealParameter fsoParam = new RealParameter(String.valueOf(fso));

        TipDatePrior prior = new TipDatePrior();
        prior.initByName(
                "tree", tree,
                "finalSampleOffset", fsoParam,
                "earlierBound", earlierBounds,
                "laterBound", laterBounds);

        Assert.assertEquals(0.0, prior.calculateLogP(), 1e-10);

        fsoParam.setValue(fso+0.99/365.25);
        prior.initAndValidate();
        Assert.assertEquals(0.0, prior.calculateLogP(), 1e-10);

        fsoParam.setValue(fso-0.99/365.25);
        prior.initAndValidate();
        Assert.assertEquals(Double.NEGATIVE_INFINITY, prior.calculateLogP(), 1e-10);

    }
}
