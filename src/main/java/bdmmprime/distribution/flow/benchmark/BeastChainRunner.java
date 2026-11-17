/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution.flow.benchmark;

import beast.base.core.ProgramStatus;
import beast.base.inference.MCMC;
import beast.base.parser.XMLParser;
import beast.base.util.Randomizer;
import beast.pkgmgmt.PackageManager;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Shared machinery for the benchmarks that run short MCMC chains of the BEAST analyses in the
 * BDMM-Flow-Supplementary repository. The tree prior of each analysis is reconfigured by rewriting
 * its XML element before parsing, so that the model code itself stays untouched.
 */
public class BeastChainRunner {

    /**
     * A BEAST analysis to run, with its XML file relative to the benchmark directory and the length
     * of the measured chains. The chain lengths differ because a step of the H3N2 analysis is far
     * more expensive than one of the MERS-CoV analysis.
     */
    record Dataset(String name, String xmlPath, long chainLength) {

        /**
         * Returns the length of the warm-up chains that let the JIT compiler optimize the hot code
         * before anything is measured.
         */
        long getWarmupChainLength() {
            return this.chainLength / 10;
        }
    }

    static final List<Dataset> DATASETS = List.of(
            new Dataset("MERS-CoV", "mers/bdmm-flow.xml", 20_000),
            new Dataset("H3N2", "h3n2/bdmm-flow.xml", 5_000)
    );

    static final String DEFAULT_BENCHMARK_DIRECTORY = "../BDMM-Flow-Supplementary/beast_benchmarks";

    /**
     * Tree prior inputs shared by all runs: the flow engine with the inverse flow and the average
     * inverse preconditioner ("IF AIP"), parallelized and with all caching enabled.
     */
    static final Map<String, String> DEFAULT_TREE_PRIOR_INPUTS = Map.of(
            "engine", "flow",
            "initialMatrixStrategy", "average_inverse",
            "useInverseFlow", "true",
            "maxConditioningNumber", "1e8",
            "useLoucaPennellIntervals", "false",
            "parallelize", "true",
            "useSubtreeCaching", "true",
            "useFlowCaching", "true"
    );

    static final String TREE_PRIOR_SPEC = "spec=\"bdmmprime.distribution.BirthDeathMigrationDistribution\"";
    static final Pattern TREE_PRIOR_ELEMENT = Pattern.compile("<\\w+[^>]*" + Pattern.quote(TREE_PRIOR_SPEC) + "[^>]*>");

    /**
     * Loads the installed BEAST packages the analyses depend on (e.g. ORC and feast) and gives
     * BEAST's own thread pool one thread per core, as {@code beast -threads -1} does.
     */
    static void setUpBeast() throws IOException {
        PackageManager.loadExternalJars();

        ProgramStatus.m_nThreads = Runtime.getRuntime().availableProcessors();
        ProgramStatus.g_exec = Executors.newFixedThreadPool(ProgramStatus.m_nThreads);
    }

    /**
     * Shuts down BEAST's thread pool so that the JVM can exit.
     */
    static void tearDownBeast() {
        ProgramStatus.g_exec.shutdown();
    }

    /**
     * Returns the directory containing the BEAST analyses, taken from the first command line
     * argument if given.
     */
    static Path getBenchmarkDirectory(String[] args) {
        return Path.of(args.length > 0 ? args[0] : DEFAULT_BENCHMARK_DIRECTORY);
    }

    /**
     * Parses the analysis with the given tree prior inputs and prepares it to run a chain of the
     * given length. Loggers are removed and the state file is written to a temporary location so
     * that the run leaves no files behind in the analysis directory.
     */
    static MCMC loadChain(Path xmlFile, Map<String, String> treePriorInputs, long chainLength, long seed) throws Exception {
        String xml = Files.readString(xmlFile);
        Path configuredXmlFile = Files.createTempFile("bdmm-benchmark", ".xml");
        Path stateFile = Files.createTempFile("bdmm-benchmark", ".state");
        configuredXmlFile.toFile().deleteOnExit();
        stateFile.toFile().deleteOnExit();
        Files.writeString(configuredXmlFile, configureTreePrior(xml, treePriorInputs));

        // the variables only appear in logger file names, which are never written
        Map<String, String> parserDefinitions = new HashMap<>(Map.of(
                "filebase", "benchmark",
                "seed", Long.toString(seed)
        ));

        Randomizer.setSeed(seed);
        MCMC mcmc = (MCMC) new XMLParser(parserDefinitions).parseFile(configuredXmlFile.toFile());

        mcmc.chainLengthInput.setValue(chainLength, mcmc);
        mcmc.loggersInput.get().clear();
        mcmc.setStateFile(stateFile.toString(), false);

        return mcmc;
    }

    /**
     * Replaces the given attributes of the tree prior element, keeping all others.
     */
    static String configureTreePrior(String xml, Map<String, String> treePriorInputs) {
        Matcher matcher = TREE_PRIOR_ELEMENT.matcher(xml);
        if (!matcher.find())
            throw new IllegalArgumentException("The XML contains no element with " + TREE_PRIOR_SPEC + ".");

        int start = matcher.start();
        int end = matcher.end();
        if (matcher.find())
            throw new IllegalArgumentException("The XML contains more than one element with " + TREE_PRIOR_SPEC + ".");

        String element = xml.substring(start, end);
        StringBuilder attributes = new StringBuilder();
        for (Map.Entry<String, String> input : treePriorInputs.entrySet()) {
            element = element.replaceAll("\\s" + Pattern.quote(input.getKey()) + "=\"[^\"]*\"", "");
            attributes.append(" %s=\"%s\"".formatted(input.getKey(), input.getValue()));
        }
        element = element.replace(TREE_PRIOR_SPEC, TREE_PRIOR_SPEC + attributes);

        return xml.substring(0, start) + element + xml.substring(end);
    }

    /**
     * Runs the chain and returns its wall-clock duration in nanoseconds, excluding parsing.
     */
    static long runChain(MCMC mcmc) throws Exception {
        long start = System.nanoTime();
        mcmc.run();
        return System.nanoTime() - start;
    }

    /**
     * Overwrites the results file with the given header line.
     */
    static void createResultsFile(Path resultsFile, String header) throws IOException {
        Files.writeString(resultsFile, header + "\n");
    }

    /**
     * Appends one line to the results file, so that completed runs survive an interruption.
     */
    static void appendResult(Path resultsFile, String line) throws IOException {
        Files.writeString(resultsFile, line + "\n", StandardOpenOption.APPEND);
    }

}
