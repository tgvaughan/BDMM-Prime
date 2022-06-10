package bdmmprime.util.priors;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import feast.fileio.TraitSetFromTaxonSet;
import junit.framework.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class TipDatePriorTest {


    @Test
    public void test() {

        List<Taxon> taxonList = new ArrayList<>();
        taxonList.add(new Taxon("t1|2020-05-15|2020-05-01|2020-06-01"));
        taxonList.add(new Taxon("t2|2020-02-15|2020-02-01|2020-03-01"));
        taxonList.add(new Taxon("t3|2021-05-15|2021-05-01|2021-06-01"));
        taxonList.add(new Taxon("t4|2020-05-15|2020-05-14|2020-05-15"));
        taxonList.add(new Taxon("t5|2020-05-15|2020-05-14|2020-05-15"));
        TaxonSet taxonSet = new TaxonSet(taxonList);

        TraitSet tipDates = new TraitSetFromTaxonSet();
        tipDates.initByName(
                "traitname", "date-forward",
                "taxa", taxonSet,
                "delimiter", "|",
                "takeGroup", 1,
                "dateFormat", "yyyy-M-dd");

        TraitSet earlierBounds = new TraitSetFromTaxonSet();
        earlierBounds.initByName(
                "traitname", "date-forward",
                "taxa", taxonSet,
                "delimiter", "|",
                "takeGroup", 2,
                "dateFormat", "yyyy-M-dd");

        TraitSet laterBounds = new TraitSetFromTaxonSet();
        laterBounds.initByName(
                "traitname", "date-forward",
                "taxa", taxonSet,
                "delimiter", "|",
                "takeGroup", 3,
                "dateFormat", "yyyy-M-dd");


        Tree tree = new Tree();
        tree.initByName("taxonset", taxonSet,
                "trait", tipDates);

        System.out.println(tree);

        RealParameter fso = new RealParameter("0.0");

        TipDatePrior prior = new TipDatePrior();
        prior.initByName(
                "tree", tree,
                "finalSampleOffset", fso,
                "initialTipDates", tipDates,
                "earlierBound", earlierBounds,
                "laterBound", laterBounds);

        Assert.assertEquals(0.0, prior.calculateLogP(), 1e-10);

        fso.setValue(0.99/365.25);
        prior.initAndValidate();
        Assert.assertEquals(0.0, prior.calculateLogP(), 1e-10);

        fso.setValue(-0.99/365.25);
        prior.initAndValidate();
        Assert.assertEquals(Double.NEGATIVE_INFINITY, prior.calculateLogP(), 1e-10);

    }
}
