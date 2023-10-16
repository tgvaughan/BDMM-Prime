/*
 * Copyright (C) 2023 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.beauti;

import bdmmprime.parameterization.SkylineParameter;
import bdmmprime.parameterization.TypeSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.layout.BorderPane;
import javafx.scene.paint.Color;
import javafx.scene.text.Text;

import java.awt.*;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;

public class EpochVisualizerPane extends BorderPane {

    Tree tree;
    TraitSet typeTraitSet;
    SkylineParameter param;

    Canvas canvas;

    boolean isScalar = false;

    final int HEADER_HEIGHT = 1;
    final int FOOTER_HEIGHT = 3;
    final int HEIGHT_PER_TYPE = 1;
    final int MARGIN_WIDTH = 2;

    public EpochVisualizerPane(Tree tree, TraitSet typeTraitSet, SkylineParameter param) {
        this.tree = tree;
        this.typeTraitSet = typeTraitSet;
        this.param = param;

        canvas = new Canvas(getWidth(), getHeight());
        setCenter(canvas);

        widthProperty().addListener(e -> canvas.setWidth(getWidth()));
        heightProperty().addListener(e -> canvas.setHeight(getHeight()));

        lines.setValue(HEADER_HEIGHT + FOOTER_HEIGHT +
                HEIGHT_PER_TYPE*(isScalar ? 1 : param.getNTypes()));
    }

    @Override
    protected void layoutChildren() {
        super.layoutChildren();

        GraphicsContext gc = canvas.getGraphicsContext2D();

        gc.clearRect(0,0,getWidth(),getHeight());

        gc.setFill(Color.ROYALBLUE);
        gc.fillRect(5, 5, 5, 5);

        boolean useAges = param.timesAreAgesInput.get();
        double origin = param.processLengthInput.get().getArrayValue();
        TypeSet typeSet = param.typeSetInput.get();

        if (origin <= 0.0) {
            gc.setStroke(Color.BLACK);
            gc.strokeText("Can't visualize: invalid origin value.", 0, getHeight()/2);
            return;
        }

        // Draw rectangles highlighting type bands

        double fontHeight = gc.getFont().getSize();

        if (param.getNTypes()>1 && !isScalar) {
            for (int rowNum=0; rowNum<param.getNTypes(); rowNum++) {

                // Draw bar

                if (rowNum % 2 == 1) {
                    gc.setFill(Color.LIGHTGRAY);
                    gc.fillRect(fontHeight*MARGIN_WIDTH,
                            getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE*(rowNum+1))*fontHeight,
                            getWidth()-2*MARGIN_WIDTH*fontHeight,
                            HEIGHT_PER_TYPE*fontHeight);
                }

                int typeIdx = param.getNTypes()-1-rowNum;

                // Draw label

                String typeName = param.typeSetInput.get().getTypeName(typeIdx);

                gc.setFill(Color.DARKGRAY);

                gc.fillText(typeName,
                        fontHeight*MARGIN_WIDTH,
                        getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE * rowNum) * fontHeight);

                gc.fillText(typeName,
                        getWidth() - fontHeight*MARGIN_WIDTH - getStringWidth(gc, typeName),
                        getHeight() - (FOOTER_HEIGHT + HEIGHT_PER_TYPE * rowNum) * fontHeight);
            }
        }


        // Draw axis

        double axisBaseY = getHeight() - fontHeight*FOOTER_HEIGHT;
        double boundaryLabelY = fontHeight;

        gc.setStroke(Color.BLACK);
        gc.strokeLine(getHorizontalPixel(gc, 0.0), axisBaseY,
                getHorizontalPixel(gc, origin), axisBaseY);

        String axisLabel = useAges
                ? "Age before most recent sample"
                : "Time relative to start of BD process";

        gc.setFill(Color.BLACK);
        gc.fillText(axisLabel,
                (getWidth()-getStringWidth(gc, axisLabel))/2.0,
                getHeight()-0.5*fontHeight);


        double delta = Math.pow(10, Math.ceil(Math.log10(origin/10)));

        gc.setStroke(Color.BLACK);
        gc.setFill(Color.BLACK);
        for (double t=0; t<=origin; t += delta) {
            gc.strokeLine(getHorizontalPixel(gc, t), axisBaseY,
                    getHorizontalPixel(gc, t),
                    axisBaseY+fontHeight/4);

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

            gc.fillText(tickLabel,
                    getHorizontalPixel(gc, t) - fontHeight/2,
                    axisBaseY+fontHeight/4+fontHeight);
        }

        // Mark Origin

        String originLabel = "Origin";
        double originPosition = getHorizontalPixel(gc, 0.0);
        gc.strokeLine(originPosition, axisBaseY, originPosition, fontHeight*HEADER_HEIGHT);
        gc.fillText(originLabel,
                originPosition-getStringWidth(gc, originLabel)/2,
                boundaryLabelY);

        // Mark Epochs

        gc.setStroke(Color.BLUE);
        gc.setFill(Color.BLUE);

        for (int epoch=0; epoch<param.getChangeCount(); epoch++) {
            String boundaryLabel = (epoch+1) + " -> " + (epoch+2);

            double changeTime = useAges
                    ? param.getChangeTimes()[param.getChangeCount()-epoch-1]
                    : param.getChangeTimes()[epoch];

            double boundaryPosition = getHorizontalPixel(gc, changeTime);
            gc.strokeLine(boundaryPosition, axisBaseY, boundaryPosition,
                    fontHeight*HEADER_HEIGHT);

            gc.fillText(boundaryLabel,
                    boundaryPosition - getStringWidth(gc, boundaryLabel)/2,
                    fontHeight*HEADER_HEIGHT);
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
            Arrays.fill(leafTimes, origin);
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

            gc.setFill(sampleCols[getEpoch(leafTimes[nodeNr]) % sampleCols.length]);

            double circleRad = fontHeight/4;
            gc.fillOval(getHorizontalPixel(gc, leafTimes[nodeNr]) - circleRad,
                    axisBaseY - fontHeight*HEIGHT_PER_TYPE/2 - fontHeight*HEIGHT_PER_TYPE*rowNum - circleRad,
                    circleRad*2, circleRad*2);
        }
    }

    int getEpoch(double time) {
        int epoch=0;

        while (epoch < param.getChangeCount() && time > param.getChangeTimes()[epoch])
            epoch ++;

        if (param.timesAreAgesInput.get())
            epoch = param.getChangeCount() - epoch;

        return epoch;
    }


    public void setScalar(boolean scalar) {
        isScalar = scalar;
    }

    private final SimpleIntegerProperty lines = new SimpleIntegerProperty(
            HEADER_HEIGHT + FOOTER_HEIGHT + HEIGHT_PER_TYPE);

    public int getLines() {
        return lines.get();
    }
    public void setLines(int value) {
        lines.set(value);
    }
    public IntegerProperty linesProperty() {
        return lines;
    }

    double getHorizontalPixel(GraphicsContext gc, double time) {
        boolean useAges = param.timesAreAgesInput.get();

        double charHeight = gc.getFont().getSize();
        double axisXStart = charHeight*2;
        double axisXEnd = getWidth() - charHeight*2;

        int scaledTime = (int)Math.round((axisXEnd-axisXStart)
                *time/param.processLengthInput.get().getArrayValue());

        return useAges
                ? axisXEnd - scaledTime
                : axisXStart + scaledTime;
    }

    private double getStringWidth(GraphicsContext gc, String string) {
        Text theText = new Text(string);
        theText.setFont(gc.getFont());
        return theText.getBoundsInLocal().getWidth();
    }
}
