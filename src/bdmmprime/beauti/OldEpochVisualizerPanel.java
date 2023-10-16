package bdmmprime.beauti;

import bdmmprime.parameterization.SkylineParameter;
import bdmmprime.parameterization.TypeSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;

import javax.swing.*;
import java.awt.*;
import java.math.BigDecimal;
import java.math.MathContext;

public class OldEpochVisualizerPanel extends JPanel {

    Tree tree;
    TraitSet typeTraitSet;
    SkylineParameter param;

    boolean isScalar = false;

    final int HEADER_HEIGHT = 1;
    final int FOOTER_HEIGHT = 3;
    final int HEIGHT_PER_TYPE = 1;
    final int MARGIN_WIDTH = 2;


    public OldEpochVisualizerPanel(Tree tree, TraitSet typeTraitSet, SkylineParameter param) {
        this.tree = tree;
        this.typeTraitSet = typeTraitSet;
        this.param = param;
    }

    @Override
    protected void paintComponent(Graphics g) {

        Graphics2D g2 = (Graphics2D)g;
        FontMetrics fm = getFontMetrics(getFont());

        boolean useAges = param.timesAreAgesInput.get();
        double origin = param.processLengthInput.get().getArrayValue();
        TypeSet typeSet = param.typeSetInput.get();

        if (origin <= 0.0) {
            g2.drawString("Can't visualize: invalid origin value.", 0, getHeight()/2);
            return;
        }

        // Draw rectangles highlighting type bands

        if (param.getNTypes()>1 && !isScalar) {
            for (int rowNum=0; rowNum<param.getNTypes(); rowNum++) {

                // Draw bar

                if (rowNum % 2 == 1) {
                    g2.setColor(Color.LIGHT_GRAY);
                    g2.fillRect(fm.getHeight()*MARGIN_WIDTH,
                            getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE*(rowNum+1))*fm.getHeight(),
                            getWidth()-2*MARGIN_WIDTH*fm.getHeight(),
                            HEIGHT_PER_TYPE*fm.getHeight());
                }

                int typeIdx = param.getNTypes()-1-rowNum;

                // Draw label

                String typeName = param.typeSetInput.get().getTypeName(typeIdx);

                g2.setColor(Color.DARK_GRAY);

                g2.drawString(typeName,
                        fm.getHeight()*MARGIN_WIDTH,
                        getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE * rowNum) * fm.getHeight());

                g2.drawString(typeName,
                        getWidth() - fm.getHeight()*MARGIN_WIDTH - fm.stringWidth(typeName),
                        getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE * rowNum) * fm.getHeight());
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

        // Draw ticks

        double delta = Math.pow(10, Math.ceil(Math.log10(origin/10)));

        Font oldFont = getFont();
        g2.setColor(Color.black);
        for (double t=0; t<=origin; t += delta) {
            g2.drawLine(getHorizontalPixel(t), axisBaseY, getHorizontalPixel(t), axisBaseY+fm.getHeight()/4);

            double val = useAges
                    ? origin - t
                    : t;

            String tickLabel;
            if (Math.ceil(val) == val)
                tickLabel = String.valueOf(Math.round(val));
            else if (Math.abs(val)<1e-10)
                tickLabel = "0";
            else
                tickLabel = new BigDecimal(val, new MathContext(4))
                        .stripTrailingZeros().toEngineeringString();

            g2.drawString(tickLabel,
                    getHorizontalPixel(t) - fm.stringWidth(tickLabel)/2,
                    axisBaseY+fm.getHeight()/4+fm.getHeight());
        }
        setFont(oldFont);

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
            if (tree.getDateTrait().getStringValue(tree.getNode(0).getID()) == null)
                tree.getDateTrait().initAndValidate();
            for (int nodeNr = 0; nodeNr < nLeaves; nodeNr++)
                leafTimes[nodeNr] = origin - tree.getDateTrait().getValue(tree.getNode(nodeNr).getID());
        } else {
            for (int nodeNr = 0; nodeNr < nLeaves; nodeNr++)
                leafTimes[nodeNr] = origin;
        }

        // Draw samples

        Color[] sampleCols = new Color[] {Color.RED, Color.ORANGE};

        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {


            int rowNum;
            if (isScalar) {
                rowNum = 0;
            } else {
                Node node = tree.getNode(nodeNr);
                String typeName = typeTraitSet.getStringValue(node.getID());
                rowNum = param.getNTypes()-1-typeSet.getTypeIndex(typeName);
            }

            g2.setColor(sampleCols[getEpoch(leafTimes[nodeNr]) % sampleCols.length]);

            int circleRad = fm.getHeight()/4;
            g2.fillOval(getHorizontalPixel(leafTimes[nodeNr]) - circleRad,
                    axisBaseY - fm.getHeight()*HEIGHT_PER_TYPE/2 - fm.getHeight()*HEIGHT_PER_TYPE*rowNum - circleRad,
                    circleRad*2, circleRad*2);
        }
    }

    int getHorizontalPixel(double time) {
        boolean useAges = param.timesAreAgesInput.get();

        int charHeight = getFontMetrics(getFont()).getHeight();
        int axisXStart = charHeight*2;
        int axisXEnd = getWidth() - charHeight*2;

        int scaledTime = (int)Math.round((axisXEnd-axisXStart)
                *time/param.processLengthInput.get().getArrayValue());

        return useAges
                ? axisXEnd - scaledTime
                : axisXStart + scaledTime;
    }

    int getEpoch(double time) {
        int epoch=0;

        while (epoch < param.getChangeCount() && time > param.getChangeTimes()[epoch])
            epoch ++;

        if (param.timesAreAgesInput.get())
            epoch = param.getChangeCount() - epoch;

        return epoch;
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
