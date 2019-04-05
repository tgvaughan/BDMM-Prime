package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineParameter;
import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import org.knowm.xchart.XChartPanel;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.internal.series.Series;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.Marker;
import org.knowm.xchart.style.markers.SeriesMarkers;

import javax.swing.*;
import java.awt.*;
import java.util.Collections;

public class EpochVisualizerPanel extends JPanel {

    Tree tree;
    TraitSet typeTraitSet;
    SkylineParameter param;

    XYChart chart;
    XChartPanel chartPanel;

    JCheckBox showVisualizerCheckBox;

    public EpochVisualizerPanel() {

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        showVisualizerCheckBox = new JCheckBox("Display visualizer");
        showVisualizerCheckBox.setSelected(false);
        showVisualizerCheckBox.addItemListener( e -> {
            chartPanel.setVisible(showVisualizerCheckBox.isSelected());
        });
        add(showVisualizerCheckBox, Component.LEFT_ALIGNMENT);

        chart = new XYChartBuilder()
                .height(getFontMetrics(getFont()).getHeight()*10)
                .width(getFontMetrics(getFont()).getHeight()*30)
                .build();

        chart.getStyler().setLegendVisible(false);

        chartPanel = new XChartPanel<>(chart);
        chartPanel.setVisible(false);
        add(chartPanel, Component.LEFT_ALIGNMENT);
    }



    public void drawChart(Tree tree, TraitSet typeTraitSet, SkylineParameter param) {

        boolean useAges = param.timesAreAgesInput.get();

        double[] leafAges = new double[tree.getLeafNodeCount()];
        if (tree.hasDateTrait()) {
            tree.getDateTrait().initAndValidate();
            for (int nodeNr = 0; nodeNr < tree.getLeafNodeCount(); nodeNr++)
                leafAges[nodeNr] = tree.getDateTrait().getValue(tree.getNode(nodeNr).getID());
        } else {
            for (int nodeNr = 0; nodeNr < tree.getLeafNodeCount(); nodeNr++)
                leafAges[nodeNr] = 0.0;

        }

        RealParameter origin = param.originInput.get();
        TypeSet typeSet = param.typeSetInput.get();

        int epoch = 1;
        for (Double epocBoundaryTime : param.getChangeTimes()) {
            String seriesName = "Epoch " + epoch + "->" + (epoch+1) + " Boundary";
            chart.addSeries(seriesName,
                    new double[] {epocBoundaryTime, epocBoundaryTime},
                    new double[] {-1, param.getNTypes()});

            XYSeries series = chart.getSeriesMap().get(seriesName);
            series.setLineColor(Color.BLACK);
            series.setMarker(SeriesMarkers.NONE);

            epoch += 1;
        }

        double originXvalue = useAges ? origin.getValue() : 0.0;

        chart.addSeries("Origin",
                new double[] {originXvalue, originXvalue},
                new double[] {-1, param.getNTypes()});
        XYSeries series = chart.getSeriesMap().get("Origin");
        series.setLineColor(Color.BLACK);
        series.setMarker(SeriesMarkers.NONE);

        int nLeaves = tree.getLeafNodeCount();

        double[] x = new double[nLeaves];
        double[] y = new double[nLeaves];


        if (useAges)
            chart.setXAxisTitle("Age before most recent sample");
        else
            chart.setXAxisTitle("Time after the origin");

        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {

            Node node = tree.getNode(nodeNr);
            String typeName = typeTraitSet.getStringValue(node.getID());
            int typeIdx = typeSet.getTypeIndex(typeName);

            x[nodeNr] = useAges ? leafAges[nodeNr] : origin.getValue() - leafAges[nodeNr];
            y[nodeNr] = typeSet.getNTypes() - (typeIdx + 1)
                    + 0.1*(Randomizer.nextDouble() - 0.5);
        }

        chart.addSeries("Samples", x, y);
        chart.getSeriesMap().get("Samples").setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        series = chart.getSeriesMap().get("Samples");
        series.setMarker(SeriesMarkers.DIAMOND);
        series.setMarkerColor(Color.ORANGE);

    }
}
