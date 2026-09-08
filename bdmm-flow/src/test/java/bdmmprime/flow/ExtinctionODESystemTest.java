package bdmmflow.flow;


import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.extinctionSystem.ExtinctionProbabilitiesODESystem;
import bdmmflow.intervals.IntervalODESystem;
import bdmmprime.distribution.P0GeSystem;
import bdmmprime.parameterization.*;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;

public class ExtinctionODESystemTest {

//    private void testExtinctionODESolver(
//            Parameterization parameterization,
//            double endTime
//    ) {
//        // setup extinction ODE
//
//        double[] initialState = new double[parameterization.getNTypes()];
//        Arrays.fill(initialState, 1.0);
//
//        IntervalODESystem extinctionSystem = new ExtinctionProbabilitiesODESystem(parameterization, 1e-100, 1e-20);
//        ExtinctionProbabilities extinctionProbabilities = new ExtinctionProbabilities(extinctionSystem.integrateBackwards(
//                initialState.clone()
//        ));
//        int intervalEndTime = parameterization.getIntervalIndex(endTime);
//        double[] flowIntegral = extinctionProbabilities.getProbability(endTime);
//
//        // setup extinction ODE in the distribution package
//
//        P0GeSystem oldSystem = new P0GeSystem(parameterization, 1e-100, 1e-20);
//
//        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
//                parameterization.getTotalProcessLength() * 1e-100,
//                parameterization.getTotalProcessLength() / 10,
//                1e-100, 1e-20
//        );
//        double[] oldIntegral = new double[]{
//                initialState[0],
//                initialState[1],
//                0, // we don't care about the actual likelihood state
//                0
//        };
//        integrator.integrate(oldSystem, parameterization.getTotalProcessLength(), oldIntegral, endTime, oldIntegral);
//        oldIntegral = Arrays.copyOfRange(oldIntegral, 0, 2);
//
//        System.out.println(Arrays.toString(oldIntegral));
//        System.out.println(Arrays.toString(flowIntegral));
//
//        for (int i = 0; i < flowIntegral.length; i++) {
//            assertEquals(oldIntegral[i], flowIntegral[i], 1e-10);
//        }
//    }

//    void testDerivatives(
//            Parameterization parameterization,
//            double t,
//            double[] initialState
//    ) {
//        // setup extinction ODE
//
//        IntervalODESystem extinctionSystem = new ExtinctionProbabilitiesODESystem(parameterization, 1e-100, 1e-20);
//
//        // get flow ODE derivatives
//
//        double[] flowDerivatives = new double[initialState.length];
//        extinctionSystem.computeDerivatives(t, initialState.clone(), flowDerivatives);
//
//        // get old ODE derivatives
//
//        P0GeSystem oldSystem = new P0GeSystem(parameterization, 1e-100, 1e-20);
//
//        double[] oldInitialState = new double[]{
//                initialState[0],
//                initialState[1],
//                0,
//                0
//        };
//        double[] oldDerivatives = new double[2 * initialState.length];
//        oldSystem.computeDerivatives(t, oldInitialState, oldDerivatives);
//
//        System.out.println(Arrays.toString(flowDerivatives));
//        System.out.println(Arrays.toString(oldDerivatives));
//
//        for (int i = 0; i < initialState.length; i++) {
//            assertEquals(oldDerivatives[i], flowDerivatives[i], 1e-10);
//        }
//    }

//    @Test
//    public void testFlowODESolver3() {
//        Parameterization parameterization = new EpiParameterization();
//        parameterization.initByName(
//                "processLength", new RealParameter("2.5"),
//                "typeSet", new TypeSet(2),
//                "R0", new SkylineVectorParameter(
//                        null,
//                        new RealParameter(Double.toString(40.0 / 3.0) + " " + Double.toString(1.0 / 3.0)), 2),
//                "becomeUninfectiousRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.5"), 2),
//                "samplingProportion", new SkylineVectorParameter(
//                        null,
//                        new RealParameter(Double.toString(1.0 / 3.0)), 2),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.0"), 2),
//                "migrationRate", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.1"),
//                        2)
//        );
//
//        this.testDerivatives(
//                parameterization, 0.5, new double[]{0.25, 0.25}
//        );
//        this.testDerivatives(
//                parameterization, 0.25, new double[]{0.4, 0.6}
//        );
//
//        this.testExtinctionODESolver(
//                parameterization,
//                1.0
//        );
//        this.testExtinctionODESolver(
//                parameterization,
//                0.0
//        );
//    }

//    @Test
//    public void testFlowODESolver4() {
//        Parameterization parameterization = new EpiParameterization();
//        parameterization.initByName(
//                "typeSet", new TypeSet(2),
//                "processLength", new RealParameter("6.0"),
//                "R0", new SkylineVectorParameter(
//                        null,
//                        new RealParameter((4.0/3.0) + " 1.1"),
//                        2),
//                "becomeUninfectiousRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.5 1.4"),
//                        2),
//                "R0AmongDemes", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.0"),
//                        2),
//                "migrationRate", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.2 0.3"),
//                        2),
//                "samplingProportion", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.33"),
//                        2),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.3 0.4"))
//        );
//
//
//        double[] initialState = new double[]{1.0, 1.0};
//
//        this.testExtinctionODESolver(
//                parameterization,
//                1.0
//        );
//    }


}