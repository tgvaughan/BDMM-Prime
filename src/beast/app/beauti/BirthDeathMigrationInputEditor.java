package beast.app.beauti;

import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.speciation.BirthDeathMigrationModel;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ItemEvent;
import java.util.ArrayList;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.TableModelEvent;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
/**
 * Created by Denise on 04.07.16.
 */
public class BirthDeathMigrationInputEditor extends InputEditor.Base {

    DefaultTableModel R0Model, deltaModel, samplingModel, samplingTimesModel, rateMatrixModel;
    SpinnerNumberModel nTypesModel;
    BirthDeathMigrationModel bdmm;

    JCheckBox R0EstCheckBox, deltaEstCheckBox, samplingEstCheckBox, rateMatrixEstCheckBox;

    boolean dimChangeInProgress = false;

    public BirthDeathMigrationInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return BirthDeathMigrationModel.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
        ExpandOption bExpandOption, boolean bAddButtons) {

        // Set up fields
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_beastObject = beastObject;
		this.itemNr = itemNr;

        // Adds label to left of input editor
        addInputLabel();

        // Create component models and fill them with data from input
        bdmm = (BirthDeathMigrationModel) ((ArrayList) input.get()).get(0);
        nTypesModel = new SpinnerNumberModel(2, 2, Short.MAX_VALUE, 1);
        R0Model = new DefaultTableModel();
        deltaModel = new DefaultTableModel();
        samplingModel = new DefaultTableModel();
        samplingTimesModel = new DefaultTableModel();
        rateMatrixModel = new DefaultTableModel() {
            @Override
            public boolean isCellEditable(int row, int column) {
                return row != column;
            }
        };
        R0EstCheckBox = new JCheckBox("estimate");
        deltaEstCheckBox = new JCheckBox("estimate");
        samplingEstCheckBox = new JCheckBox("estimate");
        rateMatrixEstCheckBox = new JCheckBox("estimate");
        loadFromBDMM();

        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(new EtchedBorder());

        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(3, 3, 3, 3);
        c.weighty = 0.5;

        // Deme count spinner:
        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Number of demes: "), c);
        JSpinner dimSpinner = new JSpinner(nTypesModel);
        dimSpinner.setMaximumSize(new Dimension(100, Short.MAX_VALUE));
        c.gridx = 1;
        c.gridy = 0;
        c.weightx = 1.0;
        c.anchor = GridBagConstraints.LINE_START;
        panel.add(dimSpinner, c);

        // Reproduction number table
        c.gridx = 0;
        c.gridy = 1;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Reproduction number per type: "), c);
        JTable R0Table = new JTable(R0Model);
        R0Table.setShowVerticalLines(true);
        R0Table.setCellSelectionEnabled(true);
        R0Table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        R0Table.setMaximumSize(new Dimension(100, Short.MAX_VALUE));

        c.gridx = 1;
        c.gridy = 1;
        c.weightx = 1.0;
        c.anchor = GridBagConstraints.LINE_START;
        panel.add(R0Table, c);
        R0EstCheckBox.setSelected(bdmm.R0.get().isEstimatedInput.get());
        c.gridx = 2;
        c.gridy = 1;
        c.anchor = GridBagConstraints.LINE_END;
        c.weightx = 1.0;
        panel.add(R0EstCheckBox, c);

        // becomeUninfectiousRate table
        c.gridx = 0;
        c.gridy = 2;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("BecomeUninfectiousRate per type: "), c);
        JTable deltaTable = new JTable(deltaModel);
        deltaTable.setShowVerticalLines(true);
        deltaTable.setCellSelectionEnabled(true);
        deltaTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        deltaTable.setMaximumSize(new Dimension(100, Short.MAX_VALUE));

        c.gridx = 1;
        c.gridy = 2;
        c.weightx = 1.0;
        c.anchor = GridBagConstraints.LINE_START;
        panel.add(deltaTable, c);
        deltaEstCheckBox.setSelected(bdmm.becomeUninfectiousRate.get().isEstimatedInput.get());
        c.gridx = 2;
        c.gridy = 2;
        c.anchor = GridBagConstraints.LINE_END;
        c.weightx = 1.0;
        panel.add(deltaEstCheckBox, c);


        // Sampling proportion table
        c.gridx = 0;
        c.gridy = 3;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("SamplingProportion per type: "), c);
        JTable samplingTable = new JTable(samplingModel);
        samplingTable.setShowVerticalLines(true);
        samplingTable.setCellSelectionEnabled(true);
        samplingTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        samplingTable.setMaximumSize(new Dimension(100, Short.MAX_VALUE));

        c.gridx = 1;
        c.gridy = 3;
        c.weightx = 1.0;
        c.anchor = GridBagConstraints.LINE_START;
        panel.add(samplingTable, c);
        samplingEstCheckBox.setSelected(bdmm.samplingProportion.get().isEstimatedInput.get());
        c.gridx = 2;
        c.gridy = 3;
        c.anchor = GridBagConstraints.LINE_END;
        c.weightx = 1.0;
        panel.add(samplingEstCheckBox, c);

        // Sampling change times table
        c.gridx = 0;
        c.gridy = 4;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Sampling change time: "), c);
        JTable samplingTimesTable = new JTable(samplingTimesModel);
        samplingTimesTable.setShowVerticalLines(true);
        samplingTimesTable.setCellSelectionEnabled(true);
        samplingTimesTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        samplingTimesTable.setMaximumSize(new Dimension(100, Short.MAX_VALUE));

        c.gridx = 1;
        c.gridy = 4;
        c.weightx = 1.0;
        c.anchor = GridBagConstraints.LINE_START;
        panel.add(samplingTimesTable, c);

        // Migration rate table
        // (Uses custom cell renderer to grey out diagonal elements.)
        c.gridx = 0;
        c.gridy = 5;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Migration rates: "), c);
        JTable rateMatrixTable = new JTable(rateMatrixModel) {
            @Override
            public TableCellRenderer getCellRenderer(int row, int column) {
                if (row != column)
                    return super.getCellRenderer(row, column);
                else
                    return new DefaultTableCellRenderer() {
                        @Override
                        public Component getTableCellRendererComponent(
                            JTable table, Object value, boolean isSelected,
                            boolean hasFocus, int row, int column) {
                            JLabel label = new JLabel();
                            label.setOpaque(true);
                            label.setBackground(Color.GRAY);
                            return label;
                        }
                    };
            }
        };
        rateMatrixTable.setShowGrid(true);
        rateMatrixTable.setCellSelectionEnabled(true);
        rateMatrixTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        rateMatrixTable.setMaximumSize(new Dimension(100, Short.MAX_VALUE));

        c.gridx = 1;
        c.gridy = 5;
        c.anchor = GridBagConstraints.LINE_START;
        c.weightx = 1.0;
        panel.add(rateMatrixTable, c);
        rateMatrixEstCheckBox.setSelected(bdmm.migrationMatrix.get().isEstimatedInput.get());

        c.gridx = 2;
        c.gridy = 5;
        c.anchor = GridBagConstraints.LINE_END;
        c.weightx = 1.0;
        panel.add(rateMatrixEstCheckBox, c);

        add(panel);


        // Event handlers

        dimSpinner.addChangeListener((ChangeEvent e) -> {
            JSpinner spinner = (JSpinner)e.getSource();
            int newDim = (int)spinner.getValue();
            int oldDim = bdmm.stateNumber.get();

            dimChangeInProgress = true;

            System.out.println("Dimension change starting.");

            int samplingIntervals = samplingModel.getColumnCount()/oldDim;

            R0Model.setColumnCount(R0Model.getColumnCount()/oldDim*newDim);
            bdmm.R0.get().setDimension(R0Model.getColumnCount());
            deltaModel.setColumnCount(deltaModel.getColumnCount()/oldDim*newDim);
            bdmm.becomeUninfectiousRate.get().setDimension(deltaModel.getColumnCount());
            samplingModel.setColumnCount(samplingIntervals*newDim);
            bdmm.samplingProportion.get().setDimension(samplingModel.getColumnCount());
            bdmm.setInputValue("stateNumber",newDim);

            StringBuilder sbfreqs = new StringBuilder();
            double fr = Math.round(100./newDim)/100.;
            for (int i=0; i<newDim; i++) {
                if (i>0)
                    sbfreqs.append(" ");

                if (i==0)   // make sure frequencies add up to 1
                    sbfreqs.append(Double.toString(Math.round(100*(1-(newDim-1)*fr))/100.));
                else
                    sbfreqs.append(Double.toString(fr));

            }
            bdmm.frequencies.get().valuesInput.setValue(
                    sbfreqs.toString(),
                    bdmm.frequencies.get());

            bdmm.setInputValue("frequencies",sbfreqs.toString());

            rateMatrixModel.setColumnCount(newDim);
            rateMatrixModel.setRowCount(newDim);
            bdmm.migrationMatrix.get().setDimension(newDim*newDim);
            for (int i=0; i<newDim; i++) {
                if (R0Model.getValueAt(0, i) == null) {
                    R0Model.setValueAt(2.0, 0, i);
                }
                if (deltaModel.getValueAt(0, i) == null) {
                    deltaModel.setValueAt(1.0, 0, i);
                }
                for (int j=0; j<samplingIntervals; j++) {
                    int k = i * samplingIntervals + j;
                    if (samplingModel.getValueAt(0, k) == null) {
                        samplingModel.setValueAt(1.0, 0, k);
                    }
                }
                for (int j=0; j<newDim; j++) {
                    if (i==j)
                        continue;
                    if (rateMatrixModel.getValueAt(j, i) == null) {
                        rateMatrixModel.setValueAt(0.1, j, i);
                    }
                }
            }

            System.out.println("Dimension change finishing.");

            dimChangeInProgress = false;

            saveToBDMM();
        });

        R0Model.addTableModelListener((TableModelEvent e) -> {
            if (e.getType() != TableModelEvent.UPDATE)
                return;

            if (!dimChangeInProgress)
                saveToBDMM();
        });

        R0EstCheckBox.addItemListener((ItemEvent e) -> {
            saveToBDMM();
        });

        deltaModel.addTableModelListener((TableModelEvent e) -> {
            if (e.getType() != TableModelEvent.UPDATE)
                return;

            if (!dimChangeInProgress)
                saveToBDMM();
        });

        deltaEstCheckBox.addItemListener((ItemEvent e) -> {
            saveToBDMM();
        });

        samplingModel.addTableModelListener((TableModelEvent e) -> {
            if (e.getType() != TableModelEvent.UPDATE)
                return;

            if (!dimChangeInProgress)
                saveToBDMM();
        });

        samplingEstCheckBox.addItemListener((ItemEvent e) -> {
            saveToBDMM();
        });

        samplingTimesModel.addTableModelListener((TableModelEvent e) -> {
            if (e.getType() != TableModelEvent.UPDATE)
                return;

            if (!dimChangeInProgress)
                saveToBDMM();
        });

        rateMatrixModel.addTableModelListener((TableModelEvent e) -> {
            if (e.getType() != TableModelEvent.UPDATE)
                return;

            if (!dimChangeInProgress)
                saveToBDMM();
        });

        rateMatrixEstCheckBox.addItemListener((ItemEvent e) -> {
            saveToBDMM();
        });
    }

    public void loadFromBDMM() {
        nTypesModel.setValue(bdmm.stateNumber.get());
        R0Model.setRowCount(1);
        R0Model.setColumnCount(bdmm.R0.get().getDimension());
        deltaModel.setRowCount(1);
        deltaModel.setColumnCount(bdmm.becomeUninfectiousRate.get().getDimension());
        samplingModel.setRowCount(1);
        samplingModel.setColumnCount(bdmm.samplingProportion.get().getDimension());
        samplingTimesModel.setRowCount(1);
        samplingTimesModel.setColumnCount(bdmm.samplingRateChangeTimesInput.get().getDimension()-1);
        rateMatrixModel.setRowCount(bdmm.stateNumber.get());   // todo: allow changes in ratMatrix as well!
        rateMatrixModel.setColumnCount(bdmm.stateNumber.get());

        for (int i=0; i<bdmm.R0.get().getDimension(); i++) {
            R0Model.setValueAt(bdmm.R0.get().getValue(i), 0, i);
        }
        for (int i=0; i<bdmm.becomeUninfectiousRate.get().getDimension(); i++) {
            deltaModel.setValueAt(bdmm.becomeUninfectiousRate.get().getValue(i), 0, i);
        }
        for (int i=0; i<bdmm.samplingProportion.get().getDimension(); i++) {
            samplingModel.setValueAt(bdmm.samplingProportion.get().getValue(i), 0, i);
        }

        for (int i=1; i<bdmm.samplingRateChangeTimesInput.get().getDimension(); i++) {
            samplingTimesModel.setValueAt(bdmm.samplingRateChangeTimesInput.get().getValue(i), 0, i-1);
        }

        for (int i=0; i<bdmm.stateNumber.get(); i++) {
            for (int j=0; j<bdmm.stateNumber.get(); j++) {
                if (i == j)
                    continue;
                rateMatrixModel.setValueAt(bdmm.getNbyNRate(i, j), i, j);
            }
        }

        R0EstCheckBox.setSelected(bdmm.R0.get().isEstimatedInput.get());
        deltaEstCheckBox.setSelected(bdmm.becomeUninfectiousRate.get().isEstimatedInput.get());
        samplingEstCheckBox.setSelected(bdmm.samplingProportion.get().isEstimatedInput.get());
        rateMatrixEstCheckBox.setSelected(bdmm.migrationMatrix.get().isEstimatedInput.get());
    }

    public void saveToBDMM() {
        StringBuilder sbR0 = new StringBuilder();
        for (int i=0; i<R0Model.getColumnCount(); i++) {
            if (i>0)
                sbR0.append(" ");

            if (R0Model.getValueAt(0, i) != null)
                sbR0.append(R0Model.getValueAt(0, i));
            else
                sbR0.append("2.0");
        }
        bdmm.R0.get().setDimension(R0Model.getColumnCount());
        bdmm.R0.get().valuesInput.setValue(
                sbR0.toString(),
                bdmm.R0.get());

        StringBuilder sbdelta = new StringBuilder();
        for (int i=0; i<deltaModel.getColumnCount(); i++) {
            if (i>0)
                sbdelta.append(" ");

            if (deltaModel.getValueAt(0, i) != null)
                sbdelta.append(deltaModel.getValueAt(0, i));
            else
                sbdelta.append("1.0");
        }
        bdmm.becomeUninfectiousRate.get().setDimension(deltaModel.getColumnCount());
        bdmm.becomeUninfectiousRate.get().valuesInput.setValue(
                sbdelta.toString(),
                bdmm.becomeUninfectiousRate.get());

        StringBuilder sbsampling = new StringBuilder();
        for (int i=0; i<samplingModel.getColumnCount(); i++) {
            if (i>0)
                sbsampling.append(" ");

            if (samplingModel.getValueAt(0, i) != null)
                sbsampling.append(samplingModel.getValueAt(0, i));
            else
                sbsampling.append("0.0 0.01");
        }
        bdmm.samplingProportion.get().setDimension(samplingModel.getColumnCount());
        bdmm.samplingProportion.get().valuesInput.setValue(
                sbsampling.toString(),
                bdmm.samplingProportion.get());

        StringBuilder sbsamplingtimes = new StringBuilder();
        for (int i=0; i<samplingTimesModel.getColumnCount(); i++) {
            if (i==0)
                sbsamplingtimes.append("0.0 ");

            if (i>0)
                sbsamplingtimes.append(" ");

            if (samplingTimesModel.getValueAt(0, i) != null)
                sbsamplingtimes.append(samplingTimesModel.getValueAt(0, i));
        }
        bdmm.samplingRateChangeTimesInput.get().setDimension(samplingTimesModel.getColumnCount()+1);
        bdmm.samplingRateChangeTimesInput.get().valuesInput.setValue(
                sbsamplingtimes.toString(),
                bdmm.samplingRateChangeTimesInput.get());


        StringBuilder sbRateMatrix = new StringBuilder();
        boolean first = true;
        for (int i=0; i<rateMatrixModel.getRowCount(); i++) {
            for (int j=0; j<rateMatrixModel.getColumnCount(); j++) {
                if (i == j)
                    continue;

                if (first)
                    first = false;
                else
                    sbRateMatrix.append(" ");

                if (rateMatrixModel.getValueAt(i, j) != null)
                    sbRateMatrix.append(rateMatrixModel.getValueAt(i, j));
                else
                    sbRateMatrix.append("0.1");
            }
        }
        bdmm.migrationMatrix.get().setDimension(
            R0Model.getColumnCount()*(R0Model.getColumnCount()-1));
        bdmm.migrationMatrix.get().valuesInput.setValue(
            sbRateMatrix.toString(),
            bdmm.migrationMatrix.get());

        bdmm.R0.get().isEstimatedInput.setValue(
            R0EstCheckBox.isSelected(), bdmm.R0.get());
        bdmm.migrationMatrix.get().isEstimatedInput.setValue(
            rateMatrixEstCheckBox.isSelected(), bdmm.migrationMatrix.get());
        bdmm.samplingProportion.get().isEstimatedInput.setValue(
                samplingEstCheckBox.isSelected(), bdmm.samplingProportion.get());

        try {
            bdmm.R0.get().initAndValidate();
            bdmm.samplingProportion.get().initAndValidate();
            bdmm.samplingRateChangeTimesInput.get().initAndValidate();
            bdmm.becomeUninfectiousRate.get().initAndValidate();
            bdmm.migrationMatrix.get().initAndValidate();
            bdmm.initAndValidate();
        } catch (Exception ex) {
            System.err.println(ex.getCause());
            System.err.println("Error updating tree prior.");
        }

        refreshPanel();
    }
}
