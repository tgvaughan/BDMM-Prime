package bdmm.mapping;

import bdmm.distributions.BirthDeathMigrationDistribution;
import bdmm.parameterization.*;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class TypeMappedTreeTest {

    @Test
    public void testBackwardIntegration() {

        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "origin", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0 / 3.0) + " " + 5.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.25")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter((1.0 / 3.0) + " " + (1.0/2.0))),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1")),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        RealParameter frequencies = new RealParameter("0.5 0.5");

        // Compute density using regular BDMM phylodynamic likelihood

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", frequencies,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        double logProbTrue = density.calculateLogP();

        TypeMappedTree typeMappedTree = new TypeMappedTree();
        typeMappedTree.initByName(
                "parameterization", parameterization,
                "frequencies", frequencies,
                "untypedTree", tree,
                "typeLabel", "type");

        // Compute density using result of backwards integration method

        double[] y = typeMappedTree.backwardsIntegrateSubtree(tree.getRoot(), 0.0);
        double logScaleFactor = typeMappedTree.geScaleFactors[tree.getRoot().getNr()];

        double logProb = 0.0;
        for (int type=0; type<parameterization.getNTypes(); type++) {
            logProb += y[type+parameterization.getNTypes()]*Math.exp(logScaleFactor)*frequencies.getValue(type);
        }

        logProb = Math.log(logProb);

        assertEquals(logProbTrue, logProb, 1e-5);
    }

    @Test
    public void testBackwardIntegrationWithRhoSampling() {

        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "origin", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0 / 3.0) + " " + 5.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.25")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter((1.0 / 3.0) + " " + (1.0/2.0))),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "rhoSampling", new TimedParameter(
                        new RealParameter("1.3 6.0"),
                        new RealParameter("0.1 0.2 0.25 0.5")));

        RealParameter frequencies = new RealParameter("0.5 0.5");

        // Compute density using regular BDMM phylodynamic likelihood

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", frequencies,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        double logProbTrue = density.calculateLogP();

        TypeMappedTree typeMappedTree = new TypeMappedTree();
        typeMappedTree.initByName(
                "parameterization", parameterization,
                "frequencies", frequencies,
                "untypedTree", tree,
                "typeLabel", "type");

        // Compute density using result of backwards integration method

        double[] y = typeMappedTree.backwardsIntegrateSubtree(tree.getRoot(), 0.0);
        double logScaleFactor = typeMappedTree.geScaleFactors[tree.getRoot().getNr()];


        double logProb = 0.0;
        for (int type=0; type<parameterization.getNTypes(); type++) {
            logProb += y[type+parameterization.getNTypes()]*frequencies.getValue(type);
        }
        logProb = Math.log(logProb) + logScaleFactor;

        assertEquals(logProbTrue, logProb, 1e-5);

    }

    @Test
    public void testForwardSimulation () {
        // uncoloured tree, 291 tips

        Tree tree = new TreeParser(
                "(((((((t1[&type=1]:0.9803361397,t2[&type=0]:0.9035540882):0.0532383481,t3[&type=0]:0.2637392259):0.6273536528,(t4[&type=1]:0.8624112266,t5[&type=0]:0.3278892266):0.2606245542):0.2941323873,(t6[&type=0]:0.09820114588,t7[&type=0]:0.533115675):0.8625875909):0.7040311908,(((t8[&type=0]:0.8696136218,t9[&type=0]:0.08719484485):0.4204288905,(t10[&type=1]:0.102143287,(t11[&type=1]:0.9850614571,t12[&type=1]:0.7407912319):0.8715072596):0.5182644848):0.524062254,(((((((t13[&type=1]:0.3981794417,(t14[&type=1]:0.03889928572,t15[&type=1]:0.5187105467):0.1127638209):0.3431177251,((t16[&type=1]:0.4239511855,t17[&type=1]:0.001895790454):0.690600364,t18[&type=1]:0.6283850113):0.4073564562):0.6862231812,(((t19[&type=1]:0.9947085041,t20[&type=0]:0.4739363373):0.1873670686,t21[&type=0]:0.151270482):0.803061039,((t22[&type=0]:0.8899249982,((t23[&type=0]:0.1329096023,t24[&type=0]:0.84205155):0.8838408566,(t25[&type=0]:0.7541888549,t26[&type=0]:0.8602364615):0.8912267659):0.771449636):0.1022819551,(((t27[&type=0]:0.3134289116,(t28[&type=0]:0.2446750235,t29[&type=0]:0.8565168788):0.8277210968):0.4307989818,((t30[&type=0]:0.2330717787,t31[&type=0]:0.4438336496):0.6521712865,(t32[&type=0]:0.2534400895,t33[&type=0]:0.7885409284):0.3051449039):0.1196702593):0.4061951274,t34[&type=0]:0.8415271267):0.4365981282):0.753448925):0.1580670979):0.04210642632,(((t35[&type=0]:0.7504386581,t36[&type=0]:0.6328390085):0.9047614154,t37[&type=0]:0.4946133171):0.2264722914,((((t38[&type=0]:0.06683212146,t39[&type=0]:0.479845396):0.9424520086,t40[&type=0]:0.894530142):0.3844042511,(((t41[&type=0]:0.5215392481,t42[&type=0]:0.2366602973):0.8142298241,(t43[&type=0]:0.2968777204,(t44[&type=0]:0.655541793,t45[&type=0]:0.8608812049):0.3564132168):0.04912991729):0.1511388237,t46[&type=0]:0.9031036345):0.1874918914):0.9690212663,(t47[&type=0]:0.07753491728,(t48[&type=0]:0.8349514075,(t49[&type=0]:0.9689748741,t50[&type=0]:0.925813166):0.4534903264):0.3571097804):0.1324767114):0.5515443345):0.3330309158):0.7202291801,((t51[&type=0]:0.6977306763,((t52[&type=0]:0.9157640305,t53[&type=0]:0.4226291834):0.5872618856,t54[&type=0]:0.2063144948):0.1422286083):0.7182746637,t55[&type=0]:0.759545143):0.7437628019):0.2425582204,((t56[&type=0]:0.4614429038,(t57[&type=0]:0.9092229386,((t58[&type=0]:0.1049408391,t59[&type=0]:0.6328130178):0.642241966,((t60[&type=0]:0.264340204,t61[&type=0]:0.5904771155):0.7333205172,(t62[&type=0]:0.9183179205,t63[&type=0]:0.1090340314):0.3010568973):0.3240860389):0.3192155454):0.1835780439):0.5942421539,t64[&type=0]:0.7931551472):0.967891278):0.06263663713,(t65[&type=0]:0.5774453548,((t66[&type=0]:0.07208712469,((t67[&type=0]:0.8918803469,t68[&type=0]:0.5110983853):0.1491188321,t69[&type=0]:0.2471361952):0.9591872343):0.3133718621,(t70[&type=0]:0.944087367,t71[&type=0]:0.7830825299):0.2284035049):0.5492361034):0.1136150162):0.002181729767):0.4548798562):0.4258609388,((((((t72[&type=0]:0.27679418,t73[&type=0]:0.5398862793):0.8871422287,(((((t74[&type=0]:0.2531923286,t75[&type=0]:0.3796772889):0.4489221217,t76[&type=0]:0.2554209188):0.3248268673,t77[&type=0]:0.5372577759):0.5699883625,t78[&type=0]:0.1656995732):0.957750936,(t79[&type=0]:0.1301121258,t80[&type=0]:0.8925942327):0.2838441601):0.5258686764):0.47825964,(t81[&type=0]:0.5749240227,((t82[&type=0]:0.9574132746,(t83[&type=0]:0.00485483068,t84[&type=0]:0.8091488208):0.1985368489):0.3703975577,(((t85[&type=0]:0.3991035291,(t86[&type=0]:0.03201846033,t87[&type=0]:0.8380640063):0.05616304209):0.8414494572,t88[&type=0]:0.6844437125):0.2426782607,((t89[&type=0]:0.7543559887,t90[&type=0]:0.7162597755):0.8230077426,t91[&type=0]:0.08967904118):0.4460245941):0.8679371702):0.51572948):0.4362259945):0.2631344711,(((t92[&type=0]:0.3353162925,((t93[&type=0]:0.4025212794,t94[&type=0]:0.0281926766):0.7965471447,t95[&type=0]:0.1145715592):0.5993301494):0.08854756854,(t96[&type=0]:0.1461353719,((t97[&type=0]:0.3158547124,t98[&type=0]:0.06653800653):0.5634025722,t99[&type=0]:0.9711292514):0.9727503664):0.7684133062):0.4824229684,((t100[&type=1]:0.06834940333,t101[&type=1]:0.7794982188):0.3453287922,(t102[&type=1]:0.627945075,t103[&type=1]:0.1914187325):0.9974814849):0.6312927424):0.04858242651):0.2845227425,((t104[&type=1]:0.6782600286,(t105[&type=1]:0.03190574702,t106[&type=1]:0.5840284519):0.03041352634):0.725893975,(((t107[&type=1]:0.9885271091,t108[&type=1]:0.07126446022):0.8419693699,t109[&type=1]:0.1546431775):0.898004594,t110[&type=1]:0.2500803664):0.1493327522):0.4266726137):0.5946582041,(t111[&type=1]:0.1395377244,(((t112[&type=1]:0.7170655408,(t113[&type=1]:0.976886861,t114[&type=1]:0.9406369971):0.7471234254):0.8065501407,((t115[&type=1]:0.1713845057,(t116[&type=1]:0.7861330248,t117[&type=1]:0.6082276558):0.8413775554):0.3245444677,t118[&type=1]:0.3892389825):0.5992471091):0.7592411407,(((t119[&type=1]:0.535931844,t120[&type=1]:0.09058958571):0.4227561057,(t121[&type=1]:0.5531579193,t122[&type=1]:0.8276180199):0.6653355309):0.0941624688,t123[&type=1]:0.3623022255):0.1494971744):0.3526274569):0.9720881658):0.8149677955):0.6065687414,((((((t124[&type=1]:0.5406888947,t125[&type=1]:0.8892341822):0.06211395678,((t126[&type=1]:0.8203180477,(t127[&type=1]:0.8536844573,t128[&type=1]:0.360511546):0.9030223228):0.9095590916,((t129[&type=1]:0.9110714826,(t130[&type=1]:0.2346256471,t131[&type=1]:0.6523390864):0.1288849309):0.7077432328,(t132[&type=1]:0.4060195235,t133[&type=1]:0.1661393729):0.3910941551):0.205704404):0.8609933471):0.3724007562,((t134[&type=1]:0.1731842053,(t135[&type=1]:0.7232482471,(t136[&type=1]:0.3883952193,((t137[&type=1]:0.6709475764,t138[&type=1]:0.0372075201):0.5473196667,(t139[&type=1]:0.8092764446,t140[&type=1]:0.4123262055):0.2000603897):0.55258787):0.2654263263):0.745555162):0.2956101163,((t141[&type=1]:0.52147611,(t142[&type=1]:0.9462005703,t143[&type=1]:0.5671354234):0.6887917654):0.362258781,t144[&type=1]:0.4798202242):0.8242726682):0.6072624433):0.695287361,((((t145[&type=1]:0.03793937969,t146[&type=1]:0.07275558705):0.3482963489,t147[&type=1]:0.1457363514):0.1479936559,(t148[&type=1]:0.7158309214,((t149[&type=1]:0.2174433649,t150[&type=1]:0.04072828358):0.4112026501,t151[&type=1]:0.6422409331):0.3413406226):0.1693999742):0.6631712937,(((t152[&type=1]:0.2706006162,t153[&type=1]:0.9267972289):0.1387761638,((((t154[&type=1]:0.2563392594,t155[&type=1]:0.3058371837):0.5946117372,t156[&type=1]:0.6161190302):0.6970871226,(t157[&type=1]:0.2388902532,(t158[&type=1]:0.9486316761,t159[&type=1]:0.215360787):0.168830334):0.03888285463):0.1640696453,t160[&type=1]:0.6803096831):0.1418975852):0.4218000816,(((t161[&type=1]:0.8702562298,t162[&type=1]:0.9289729816):0.05807372741,t163[&type=1]:0.3533785399):0.5012762842,(((t164[&type=1]:0.8666574673,t165[&type=1]:0.9603798252):0.7887994377,t166[&type=1]:0.857058729):0.4139410679,(t167[&type=1]:0.5900272813,t168[&type=1]:0.3345388798):0.06017537019):0.9609203783):0.7103463742):0.696603697):0.6451920038):0.1909481271,((((t169[&type=1]:0.9171597108,t170[&type=1]:0.9479122513):0.7170342554,(t171[&type=1]:0.2722596873,((t172[&type=1]:0.1194724559,(t173[&type=1]:0.03922236571,t174[&type=1]:0.6290624789):0.07739861775):0.8598598302,(t175[&type=1]:0.2009421999,(t176[&type=1]:0.06154947914,t177[&type=1]:8.997193072E-4):0.04738179315):0.3235510678):0.3443877005):0.6351028818):0.5525081949,((((t178[&type=1]:0.7599076207,t179[&type=1]:0.2997759853):0.5921433992,t180[&type=1]:0.7098581635):0.3725496214,(t181[&type=1]:0.5053773888,(t182[&type=1]:0.5991492711,(t183[&type=1]:0.5036820578,t184[&type=1]:0.6361607853):0.510631816):0.9604382808):0.2464167587):0.6073093358,(((t185[&type=1]:0.03128415369,(t186[&type=1]:0.5260852403,(t187[&type=1]:0.878767435,t188[&type=1]:0.4992109234):0.5333148066):0.00347468094):0.5590308013,t189[&type=1]:0.3710992143):0.5034162949,(t190[&type=1]:0.778916508,((t191[&type=1]:0.3069154553,(((t192[&type=1]:0.9946115273,t193[&type=1]:0.9138687006):0.5209144899,t194[&type=1]:0.5152770842):0.9462409306,t195[&type=1]:0.7395236609):0.4110851623):0.930918345,(((t196[&type=1]:0.7895439987,((t197[&type=1]:0.4697002599,t198[&type=1]:0.1383787312):0.6911794308,(t199[&type=1]:0.8664436699,t200[&type=0]:0.1959039853):0.8656513852):0.3620497067):0.2839249384,(t201[&type=0]:0.6558795469,t202[&type=0]:0.2103423763):0.969477433):0.9058840063,(t203[&type=0]:0.0856692954,t204[&type=0]:0.4175976661):0.820434629):0.5355881769):0.2263581599):0.4512835185):0.7323478526):0.2479199937):0.1964542414,((t205[&type=0]:0.7537573762,(t206[&type=0]:0.1392466244,(t207[&type=0]:0.5136175761,(t208[&type=0]:0.7852529553,t209[&type=0]:0.07355738804):0.1220811389):0.7572090242):0.1422528555):0.5948274662,(((((t210[&type=0]:0.3068353184,(t211[&type=0]:0.3314456891,((t212[&type=0]:0.5265486804,t213[&type=0]:0.1382007354):0.1814086549,t214[&type=0]:0.9276472756):0.07718444197):0.03486835537):0.1617580003,(t215[&type=0]:0.3328830956,t216[&type=0]:0.8558843595):0.8366736979):0.347376487,t217[&type=0]:0.8222538356):0.2337225529,(t218[&type=0]:0.06199815008,t219[&type=0]:0.45975962):0.179990889):0.0635867205,(t220[&type=0]:0.3214025751,(t221[&type=0]:0.5022090652,t222[&type=0]:0.6454557138):0.6956466341):0.2711792416):0.1847200533):0.1051658324):0.4945860899):0.936143348,(((t223[&type=0]:0.06268779701,((t224[&type=0]:0.3337278806,t225[&type=0]:0.1570303424):0.3089733059,(t226[&type=0]:0.5069784883,t227[&type=0]:0.1434204187):0.2001587199):0.04750720505):0.3600859912,((((t228[&type=0]:0.9994731578,(t229[&type=0]:0.8934116936,t230[&type=0]:0.03698333143):0.8173468311):0.3089058488,((((t231[&type=0]:0.3216121283,t232[&type=0]:0.5232846253):0.8687884973,(t233[&type=0]:0.6280638413,((t234[&type=0]:0.6543256822,t235[&type=0]:0.8677638234):0.8895299246,t236[&type=0]:0.4047793006):0.7147388768):0.3533478715):0.9470084386,t237[&type=0]:0.7769409856):0.4955915695,((t238[&type=0]:0.2772087415,(t239[&type=0]:0.4904922615,(t240[&type=0]:0.05356206303,t241[&type=0]:0.08998329984):0.8154862223):0.5610961432):0.1617916438,(t242[&type=0]:0.5707751412,(t243[&type=0]:0.9836868793,t244[&type=0]:0.1984052949):0.6953297216):0.05552111682):0.9476150468):0.2473166997):0.9623488116,((t245[&type=0]:0.7935025664,t246[&type=0]:0.08509867964):0.3953444003,(t247[&type=0]:0.09163277131,(t248[&type=0]:0.5201428954,t249[&type=0]:0.8055520628):0.7452739514):0.3989078877):0.07581191277):0.9779064963,(((t250[&type=0]:0.943611098,(t251[&type=0]:0.33392801,t252[&type=0]:0.5996331484):0.4291575127):0.4906436009,((((t253[&type=0]:0.7749450852,(t254[&type=0]:0.8616885878,t255[&type=0]:0.585028409):0.06060880423):0.1238881133,((t256[&type=0]:0.7451687793,t257[&type=0]:0.6925335305):0.05338745634,t258[&type=0]:0.3357626374):0.2069296469):0.09644073155,((((t259[&type=0]:0.2258843291,t260[&type=0]:0.2671526412):0.3940743534,(t261[&type=0]:0.5022506947,(t262[&type=0]:0.9498897423,t263[&type=0]:0.1406114365):0.2847759123):0.04320593993):0.6982026948,t264[&type=0]:0.2693712024):0.959781138,(((t265[&type=0]:0.6035173486,t266[&type=0]:0.5529949202):0.9900399651,(t267[&type=0]:0.5455351078,t268[&type=0]:0.3530619899):0.4626278321):0.2735997427,(t269[&type=0]:0.9580646451,(t270[&type=0]:0.3280033092,t271[&type=0]:0.7206294278):0.03739526332):0.4967516926):0.9350089293):0.4371789068):0.1014483059,t272[&type=0]:0.2867298371):0.07522285799):0.06352435821,((t273[&type=0]:0.4001782183,t274[&type=0]:0.7190070178):0.1696753846,(t275[&type=0]:0.5535608665,t276[&type=0]:0.01324651297):0.2691543309):0.8676247413):0.8461736294):0.1769516913):0.344365149,(((t277[&type=0]:0.3245107541,(t278[&type=0]:0.4142541443,t279[&type=0]:0.5857141651):0.819547887):0.0867733527,(t280[&type=0]:0.4938162852,(t281[&type=0]:0.2444119717,t282[&type=0]:0.08141433029):0.05381231918):0.8375963389):0.176160393,((t283[&type=0]:0.4199601968,t284[&type=0]:0.8354801824):0.3150380594,(((t285[&type=0]:0.9818797186,(t286[&type=0]:0.8971825438,((t287[&type=0]:0.5155417006,t288[&type=0]:0.8260786769):0.7060374152,t289[&type=0]:0.6001661876):0.4120474763):0.9949228324):0.8038698458,t290[&type=0]:0.1939124272):0.6380942846,t291[&type=0]:0.3665255161):0.459349304):0.482901911):0.4833473735):0.5903116504):0.9973697898);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "origin", new RealParameter(Double.toString(tree.getRoot().getHeight()+0.1)),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0 / 3.0) + " " + 5.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.25")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter((1.0 / 3.0) + " " + 0.5)),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1")),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2)

        );

        RealParameter frequencies = new RealParameter("0.5 0.5");

        // Compute density using regular BDMM phylodynamic likelihood

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", frequencies,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);

        double logProbTrue = density.calculateLogP();

        TypeMappedTree typeMappedTree = new TypeMappedTree();
        typeMappedTree.initByName(
                "parameterization", parameterization,
                "frequencies", frequencies,
                "untypedTree", tree,
                "typeLabel", "type");

        // Compute density using result of backwards integration method

        double[] y = typeMappedTree.backwardsIntegrateSubtree(tree.getRoot(), 0.0);
        double logScaleFactor = typeMappedTree.geScaleFactors[tree.getRoot().getNr()];

        double logProb = 0.0;
        for (int type=0; type<parameterization.getNTypes(); type++) {
            logProb += y[type+parameterization.getNTypes()]*frequencies.getValue(type);
        }
        logProb = Math.log(logProb) + logScaleFactor;


        assertEquals(logProbTrue, logProb, 1e-3);
    }
}
