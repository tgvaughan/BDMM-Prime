/*
 * Copyright (C) 2023 ETH Zurich
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
import bdmmprime.util.ProcessLength;
import beast.base.evolution.tree.*;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.inference.StateNodeInitialiser;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.Text;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.Optional;

public class EpochVisualizerPane extends Canvas {

    Tree tree;
    TraitSet typeTraitSet;
    SkylineParameter param;

    boolean isScalar = false;

    final int HEADER_HEIGHT = 1;
    final int FOOTER_HEIGHT = 3;
    final int HEIGHT_PER_TYPE = 1;
    final int MARGIN_WIDTH = 2;

    public EpochVisualizerPane(Tree tree, TraitSet typeTraitSet, SkylineParameter param) {
        this.tree = tree;
        this.typeTraitSet = typeTraitSet;
        this.param = param;

        widthProperty().addListener(e -> repaintCanvas());
    }

    public void repaintCanvas() {
        setHeight((HEADER_HEIGHT + FOOTER_HEIGHT
                + (isScalar ? 1 : param.getNTypes()) * HEIGHT_PER_TYPE)
                * Font.getDefault().getSize());

        GraphicsContext gc = getGraphicsContext2D();

        gc.clearRect(0,0,getWidth(),getHeight());

        reinitTree(); // Ensure tree matches any constraints

        boolean useAges = param.timesAreAgesInput.get();
        double processLength = param.processLengthInput.get().getArrayValue();
        TypeSet typeSet = param.typeSetInput.get();

        if (processLength <= 0.0) {
            gc.setStroke(Color.BLACK);
            gc.strokeText("Can't visualize: invalid process length value.", 0, getHeight()/2);
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
        gc.strokeLine(getHorizontalPixel(gc, 0.0, processLength), axisBaseY,
                getHorizontalPixel(gc, processLength, processLength), axisBaseY);

        String axisLabel = useAges
                ? "Age before most recent sample"
                : "Time relative to start of BD process";

        gc.setFill(Color.BLACK);
        gc.fillText(axisLabel,
                (getWidth()-getStringWidth(gc, axisLabel))/2.0,
                getHeight()-0.5*fontHeight);


        double delta = Math.pow(10, Math.ceil(Math.log10(processLength/10)));

        gc.setStroke(Color.BLACK);
        gc.setFill(Color.BLACK);
        for (double t=0; t<=processLength; t += delta) {
            gc.strokeLine(getHorizontalPixel(gc, t, processLength), axisBaseY,
                    getHorizontalPixel(gc, t, processLength),
                    axisBaseY+fontHeight/4);

            double val = useAges
                    ? processLength - t
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
                    getHorizontalPixel(gc, t, processLength) - fontHeight/2,
                    axisBaseY+fontHeight/4+fontHeight);
        }

        // Mark Origin

        String processStartLabel = ((ProcessLength)(param.processLengthInput.get())).treeInput.get() != null
                ? "Root" : "Origin";

        double originPosition = getHorizontalPixel(gc, 0.0, processLength);
        gc.strokeLine(originPosition, axisBaseY, originPosition, fontHeight*HEADER_HEIGHT);
        gc.fillText(processStartLabel,
                originPosition-getStringWidth(gc, processStartLabel)/2,
                boundaryLabelY);

        // Mark Epochs

        gc.setStroke(Color.BLUE);
        gc.setFill(Color.BLUE);

        for (int epoch=0; epoch<param.getChangeCount(); epoch++) {
            String boundaryLabel = (epoch+1) + " -> " + (epoch+2);

            double changeTime = useAges
                    ? param.getChangeTimes()[param.getChangeCount()-epoch-1]
                    : param.getChangeTimes()[epoch];

            double boundaryPosition = getHorizontalPixel(gc, changeTime, processLength);
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
                leafTimes[nodeNr] = processLength - tree.getDateTrait().getValue(tree.getNode(nodeNr).getID());
        } else {
            Arrays.fill(leafTimes, processLength);
        }

        // Draw samples

        Color[] sampleCols = new Color[] {Color.RED, Color.ORANGE};

        for (int nodeNr=0; nodeNr<nLeaves; nodeNr++) {

            int rowNum;
            if (isScalar) {
                rowNum = 0;
            } else {
                Node node = tree.getNode(nodeNr);
                String typeName = typeTraitSet.getStringValue(tree.getTaxonId(node));
                rowNum = param.getNTypes()-1-typeSet.getTypeIndex(typeName);
            }

            gc.setFill(sampleCols[getEpoch(leafTimes[nodeNr]) % sampleCols.length]);

            double circleRad = fontHeight/4;
            gc.fillOval(getHorizontalPixel(gc, leafTimes[nodeNr], processLength) - circleRad,
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
        repaintCanvas();
    }

    double getHorizontalPixel(GraphicsContext gc, double time, double processLength) {
        boolean useAges = param.timesAreAgesInput.get();

        double charHeight = gc.getFont().getSize();
        double axisXStart = charHeight*2;
        double axisXEnd = getWidth() - charHeight*2;

        int scaledTime = (int)Math.round((axisXEnd-axisXStart)*time/processLength);

        return useAges
                ? axisXEnd - scaledTime
                : axisXStart + scaledTime;
    }

    private double getStringWidth(GraphicsContext gc, String string) {
        Text theText = new Text(string);
        theText.setFont(gc.getFont());
        return theText.getBoundsInLocal().getWidth();
    }

    /**
     * Ensure any initialisers are applied to the tree before using its root
     * age to determine process length.
     */
    private void reinitTree() {
        ProcessLength procLengthObj = (ProcessLength) param.processLengthInput.get();
        if (procLengthObj.treeInput.get() != null) {
            Tree tree = procLengthObj.treeInput.get();
            Optional<StateNodeInitialiser> maybeInitialiser =
                    tree.getOutputs().stream()
                            .filter(bo -> bo instanceof StateNodeInitialiser)
                            .map(bo -> (StateNodeInitialiser)bo)
                            .findFirst();

            if (maybeInitialiser.isPresent()) {
                StateNodeInitialiser initialiser = maybeInitialiser.get();

                if (initialiser instanceof RandomTree ||
                        initialiser instanceof ClusterTree ||
                        initialiser instanceof TreeParser) {
                    initialiser.initStateNodes();
                }

            }
        }
    }
}
