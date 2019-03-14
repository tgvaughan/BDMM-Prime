package beast.app.bdmmprime.beauti;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.parameterization.SkylineVectorParameter;
import beast.app.beauti.BeautiDoc;
import beast.app.beauti.BeautiSubTemplate;
import beast.app.draw.BEASTObjectInputEditor;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.util.List;

public class SkylineVectorInputEditor extends InputEditor.Base {

    SkylineVectorParameter skylineVectorParameter;


    public SkylineVectorInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return SkylineVectorParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {

        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr= itemNr;

        skylineVectorParameter = (SkylineVectorParameter)input.get();

        addInputLabel();

        Box boxVert, boxHoriz;

        boxVert = Box.createVerticalBox();
        boxVert.setBorder(new EtchedBorder());

        boxHoriz = Box.createHorizontalBox();
        JLabel changePointLabel = new JLabel("Number of change points:");
        JSpinner changePointsSpinner = new JSpinner();
        changePointsSpinner.setModel(
                new SpinnerNumberModel(0, 0, Integer.MAX_VALUE,1));
        boxHoriz.add(changePointLabel);
        boxHoriz.add(changePointsSpinner);

        boxVert.add(boxHoriz);

        boxHoriz = Box.createHorizontalBox();
        JCheckBox isScalar = new JCheckBox("Scalar rates");
        boxHoriz.add(isScalar);

        boxVert.add(boxHoriz);

        add(boxVert);

//        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    }
}
