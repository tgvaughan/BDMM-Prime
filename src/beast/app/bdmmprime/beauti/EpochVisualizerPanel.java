package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.SkylineParameter;
import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import javax.swing.*;
import java.awt.*;

public class EpochVisualizerPanel extends JPanel {

    Tree tree;
    TraitSet typeTraitSet;
    SkylineParameter param;

    boolean isScalar = false;

    final int HEADER_HEIGHT = 1;
    final int FOOTER_HEIGHT = 2;
    final int HEIGHT_PER_TYPE = 1;
    final int MARGIN_WIDTH = 2;


    public EpochVisualizerPanel(Tree tree, TraitSet typeTraitSet, SkylineParameter param) {
        this.tree = tree;
        this.typeTraitSet = typeTraitSet;
        this.param = param;
    }

    @Override
    protected void paintComponent(Graphics g) {

        Graphics2D g2 = (Graphics2D)g;
        FontMetrics fm = getFontMetrics(getFont());

        boolean useAges = param.timesAreAgesInput.get();
        double origin = param.originInput.get().getValue();
        TypeSet typeSet = param.typeSetInput.get();

        // Draw rectangles highlighting type bands

        if (param.getNTypes()>1 && !isScalar) {
            for (int typeIdx=0; typeIdx<param.getNTypes(); typeIdx++) {

                // Draw bar
                if (typeIdx % 2 == 0) {
                    g2.setColor(Color.LIGHT_GRAY);
                    g2.fillRect(fm.getHeight()*MARGIN_WIDTH,
                            getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE*(typeIdx+1))*fm.getHeight(),
                            getWidth()-2*MARGIN_WIDTH*fm.getHeight(),
                            HEIGHT_PER_TYPE*fm.getHeight());
                }

                // Draw label
                g2.setColor(Color.DARK_GRAY);
                g2.drawString(param.typeSetInput.get().getTypeName(typeIdx),
                        fm.getHeight()*MARGIN_WIDTH,
                        getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE*typeIdx)*fm.getHeight()
                                - (HEIGHT_PER_TYPE-1)*fm.getHeight()/2);
            }
        }

        // Draw axis

        int axisBaseY = getHeight() - fm.getHeight()*FOOTER_HEIGHT;
        int boundaryLabelY = fm.getHeight();

        g2.setColor(Color.BLACK);
        g2.drawLine(getHorizontalPixel(0.0), axisBaseY,
                getHorizontalPixel(origin), axisBaseY);

        String axisLabel = useAges
                ? "Age before most recent sample"
                : "Time relative to start of BD process";

        g2.drawString(axisLabel,
                (int) ((getWidth()-fm.stringWidth(axisLabel))/2.0),
                (int) (getHeight()-0.5*fm.getHeight()));

        // Mark Origin

        String originLabel = "Origin";
        int originPosition = getHorizontalPixel(0.0);
        g2.drawLine(originPosition, axisBaseY, originPosition, fm.getHeight()*HEADER_HEIGHT);
        g2.drawString(originLabel,
                originPosition-fm.stringWidth(originLabel)/2,
                boundaryLabelY);

        // Mark Epochs

        for (int epoch=0; epoch<param.getChangeCount(); epoch++) {
            String boundaryLabel = (epoch+1) + " -> " + (epoch+2);
            g2.setColor(Color.BLUE);

            double changeTime = useAges
                    ? param.getChangeTimes()[param.getChangeCount()-epoch-1]
                    : param.getChangeTimes()[epoch];

            int boundaryPosition = getHorizontalPixel(changeTime);
            g2.drawLine(boundaryPosition, axisBaseY, boundaryPosition, fm.getHeight()*HEADER_HEIGHT);

            g2.drawString(boundaryLabel,
                    boundaryPosition - fm.stringWidth(boundaryLabel)/2,
                    fm.getHeight()*HEADER_HEIGHT);
        }

        // Compute sample ages

        int nLeaves = tree.getLeafNodeCount();
        double[] leafTimes = new double[nLeaves];
        if (tree.hasDateTrait()) {
            tree.getDateTrait().initAndValidate();
            for (int nodeNr = 0; nodeNr < nLeaves; nodeNr++)
                leafTimes[nodeNr] = origin - tree.getDateTrait().getValue(tree.getNode(nodeNr).getID());
        } else {
            for (int nodeNr = 0; nodeNr < nLeaves; nodeNr++)
                leafTimes[nodeNr] = origin;
        }

        // Draw samples

        g2.setColor(Color.RED);
        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {


            int typeIdx;
            if (isScalar) {
                typeIdx = 0;
            } else {
                Node node = tree.getNode(nodeNr);
                String typeName = typeTraitSet.getStringValue(node.getID());
                typeIdx = typeSet.getTypeIndex(typeName);
            }

            int circleRad = fm.getHeight()/4;
            g2.fillOval(getHorizontalPixel(leafTimes[nodeNr]) - circleRad,
                    axisBaseY - fm.getHeight()*HEIGHT_PER_TYPE/2 - fm.getHeight()*HEIGHT_PER_TYPE*typeIdx - circleRad,
                    circleRad*2, circleRad*2);
        }
    }

    int getHorizontalPixel(double time) {
        boolean useAges = param.timesAreAgesInput.get();

        int charHeight = getFontMetrics(getFont()).getHeight();
        int axisXStart = charHeight*2;
        int axisXEnd = getWidth() - charHeight*2;

        int scaledTime = (int)Math.round((axisXEnd-axisXStart)
                *time/param.originInput.get().getValue());

        return useAges
                ? axisXEnd - scaledTime
                : axisXStart + scaledTime;
    }

    public void setScalar(boolean isScalar) {
        this.isScalar = isScalar;
    }

    @Override
    public Dimension getPreferredSize() {
        int units = HEADER_HEIGHT + FOOTER_HEIGHT
                + HEIGHT_PER_TYPE*(isScalar ? 1 : param.getNTypes());

        return new Dimension(super.getPreferredSize().width,
                units*getFontMetrics(getFont()).getHeight());
    }
}
