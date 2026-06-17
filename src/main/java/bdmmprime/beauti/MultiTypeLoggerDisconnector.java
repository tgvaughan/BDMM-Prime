/*
 * Copyright (C) 2019-2025 ETH Zurich
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

import bdmmprime.mapping.TypedTreeStatsLogger;
import bdmmprime.parameterization.TypeSet;
import bdmmprime.util.OptionalLogger;
import beast.base.core.BEASTInterface;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beastfx.app.inputeditor.BeautiDoc;

import java.util.List;

public class MultiTypeLoggerDisconnector {

    public static void customConnector(BeautiDoc doc) {
        List<BEASTInterface> treePartitions = doc.getPartitions("TreeModel");
        MCMC mcmc = (MCMC) doc.pluginmap.get("mcmc");

        System.out.println("MultiTypeLoggerDisconnector.customConnector:\n" +
                "Checking for multi-type loggers applied to single-type trees.");
        for (BEASTInterface treePartition : treePartitions) {
            System.out.println("Processing partition " + treePartition.getID() + ":");

            String partitionID = BeautiDoc.parsePartition(treePartition.getID());
            if (doc.pluginmap.containsKey("typeSetBDMMPrime.t:" + partitionID)) {
                TypeSet typeSet = (TypeSet) doc.pluginmap.get("typeSetBDMMPrime.t:" + partitionID);

                if (typeSet != null && typeSet.getNTypes()>1) {
                    System.out.println("Multi-type model detected, skipping.");
                    continue;
                } else {
                    System.out.println("Single-type model detected, attempting disconnections.");
                }

                // Disconnect typed tree and tree stats loggers
                if (doc.pluginmap.containsKey("typedTreeLogger.t:" + partitionID)) {
                    OptionalLogger logger = (OptionalLogger) doc.pluginmap.get("typedTreeLogger.t:" + partitionID);
                    mcmc.loggersInput.get().remove(logger);
                    System.out.println("Disconnected typedTreeLogger.t:" + partitionID);
                }
                if (doc.pluginmap.containsKey("nodeTypedTreeLogger.t:" + partitionID)) {
                    OptionalLogger logger = (OptionalLogger) doc.pluginmap.get("nodeTypedTreeLogger.t:" + partitionID);
                    mcmc.loggersInput.get().remove(logger);
                    System.out.println("Disconnected nodeTypedTreeLogger.t:" + partitionID);
                }
                if (doc.pluginmap.containsKey("typedTreeStats.t:" + partitionID)) {
                    TypedTreeStatsLogger logger = (TypedTreeStatsLogger) doc.pluginmap.get("typedTreeStats.t:" + partitionID);
                    Logger tracelog = (Logger) doc.pluginmap.get("tracelog");
                    tracelog.loggersInput.get().remove(logger);
                    System.out.println("Disconnected typedTreeStats.t:" + partitionID);
                }
            }
        }

        System.out.println("MultiTypeLoggerDisconnector.customConnector done.");
    }
}
