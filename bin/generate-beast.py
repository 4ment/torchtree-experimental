#!/usr/bin/env python

import argparse
import re

from dendropy import Tree

parser = argparse.ArgumentParser(description="script generating BEAST XML")
parser.add_argument(
    "--input", type=argparse.FileType("r"), required=True, help="""XML template file"""
)
parser.add_argument(
    "--output",
    type=argparse.FileType("w"),
    required=True,
    help="""beast xml file""",
)
parser.add_argument("--stem", required=True, help="""stem file""")
parser.add_argument("--model", required=True, choices=["JC69", "HKY", "GTR"])
parser.add_argument(
    "--coalescent",
    required=True,
    choices=["constant", "exponential", "skyride", "skygrid"],
)
parser.add_argument("--clock", required=True, choices=["strict", "ucln"])
parser.add_argument("--site", required=True, choices=["gamma", "weibull", "invariant"])
parser.add_argument(
    "--ngen",
    default=1000000,
    type=int,
    help="""number of iterations per chain""",
)
parser.add_argument(
    "--categories",
    default=1,
    type=int,
    help="""number of categories for gamma distribution""",
)
parser.add_argument(
    "--samplefreq",
    default=1000,
    type=int,
)
parser.add_argument(
    "--grid",
    type=int,
)
parser.add_argument(
    "--cutoff",
    type=float,
)
parser.add_argument(
    "--inv",
    action="store_true",
    help="""use proportion of invariant site model""",
)
parser.add_argument(
    "--topology",
    help="""path to containing tree topology for fixed topology inference""",
)

args = parser.parse_args()

CONSTANT_SIZE = """
	<constantSize id="constant" units="years">
		<populationSize>
			<parameter id="constant.popSize" value="1.0" lower="0.0"/>
		</populationSize>
	</constantSize>
"""

STRICT_CLOCK = """
	<strictClockBranchRates id="branchRates">
		<rate>
			<parameter id="clock.rate" value="1.0" lower="0.0"/>
		</rate>
	</strictClockBranchRates>
	
	<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<strictClockBranchRates idref="branchRates"/>
	</rateStatistic>
"""
STRICT_OPERATOR = """
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="clock.rate"/>
		</scaleOperator>
		<upDownOperator scaleFactor="0.75" weight="3">
			<up>
				<parameter idref="treeModel.allInternalNodeHeights"/>
			</up>
			<down>
				<parameter idref="clock.rate"/>
			</down>
		</upDownOperator>
"""
STRICT_PRIOR = """
				<ctmcScalePrior>
					<ctmcScale>
						<parameter idref="clock.rate"/>
					</ctmcScale>
					<treeModel idref="treeModel"/>
				</ctmcScalePrior>
"""
STRICT_LOG = '<parameter idref="clock.rate"/>'
STRICT_LOG_TREE = """
			<trait name="rate" tag="rate">
				<strictClockBranchRates idref="branchRates"/>
			</trait>		
"""

UCLN_CLOCK = """
	<!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->
	<discretizedBranchRates id="branchRates">
		<treeModel idref="treeModel"/>
		<distribution>
			<logNormalDistributionModel meanInRealSpace="true">
				<mean>
					<parameter id="ucld.mean" value="1.0" lower="0.0"/>
				</mean>
				<stdev>
					<parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0"/>
				</stdev>
			</logNormalDistributionModel>
		</distribution>
		<rateCategories>
			<parameter id="branchRates.categories"/>
		</rateCategories>
	</discretizedBranchRates>
	
	<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
	
	<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
	
	<rateCovarianceStatistic id="covariance" name="covariance">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateCovarianceStatistic>
"""

UCLN_OPERATOR = """
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="ucld.mean"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="ucld.stdev"/>
		</scaleOperator>
		<upDownOperator scaleFactor="0.75" weight="3">
			<up>
				<parameter idref="treeModel.allInternalNodeHeights"/>
			</up>
			<down>
				<parameter idref="ucld.mean"/>
			</down>
		</upDownOperator>
		<swapOperator size="1" weight="10" autoOptimize="false">
			<parameter idref="branchRates.categories"/>
		</swapOperator>
		<uniformIntegerOperator weight="10">
			<parameter idref="branchRates.categories"/>
		</uniformIntegerOperator>
"""

UCLN_LOG = """
			<rateStatistic idref="meanRate"/>
			<parameter idref="ucld.mean"/>
			<parameter idref="ucld.stdev"/>
			<rateStatistic idref="coefficientOfVariation"/>
			<rateCovarianceStatistic idref="covariance"/>			
"""

UCLN_LOG_TREE = '<discretizedBranchRates idref="branchRates"/>'


HKY = """
	<HKYModel id="substmodel">
		<frequencies>
			<frequencyModel dataType="nucleotide">
				<frequencies>
					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>
				</frequencies>
			</frequencyModel>
		</frequencies>
		<kappa>
			<parameter id="kappa" value="1.0" lower="0.0"/>
		</kappa>
	</HKYModel>
"""
HKY_OPERATOR = """
		<scaleOperator scaleFactor="0.75" weight="1">
			<parameter idref="kappa"/>
		</scaleOperator>
		<deltaExchange delta="0.01" weight="1">
			<parameter idref="frequencies"/>
		</deltaExchange>
"""

HKY_PRIOR = """
				<logNormalPrior mu="1.0" sigma="1.25" offset="0.0">
					<parameter idref="kappa"/>
				</logNormalPrior>				
				<dirichletPrior alpha="1.0" sumsTo="1.0">
					<parameter idref="frequencies"/>
				</dirichletPrior>
"""
HKY_LOG = """
			<parameter idref="kappa"/>
			<parameter idref="frequencies"/>
"""
GTR = """
	<!-- The general time reversible (GTR) substitution model                    -->
	<gtrModel id="substmodel">
		<frequencies>
			<frequencyModel dataType="nucleotide">
				<frequencies>
					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>
				</frequencies>
			</frequencyModel>
		</frequencies>
		<rates>
			<parameter id="gtr.rates" dimension="6" value="1.0" lower="0.0"/>
		</rates>
	</gtrModel>
"""
GTR_OPERATOR = """
		<deltaExchange delta="0.01" weight="1">
			<parameter idref="gtr.rates"/>
		</deltaExchange>
		<deltaExchange delta="0.01" weight="1">
			<parameter idref="frequencies"/>
		</deltaExchange>
"""
GTR_PRIOR = """
				<dirichletPrior alpha="1.0" sumsTo="6.0">
					<parameter idref="gtr.rates"/>
				</dirichletPrior>
				<dirichletPrior alpha="1.0" sumsTo="1.0">
					<parameter idref="frequencies"/>
				</dirichletPrior>
"""
GTR_LOG = """
			<parameter idref="gtr.rates"/>
			<parameter idref="frequencies"/>
"""

CONSTANT_SITE = """
	<siteModel id="siteModel">
		<substitutionModel>
			<HKYModel idref="substmodel"/>
		</substitutionModel>
	</siteModel>
"""

WEIBULL_SITE = """
	<siteModelWeibull id="siteModel">
		<substitutionModel>
			<HKYModel idref="substmodel"/>
		</substitutionModel>
		<gammaShape weibullCategories="CATEGORIES">
			<parameter id="alpha" value="0.5" lower="0.0"/>
		</gammaShape>
	</siteModelWeibull>
"""

GAMMA_SITE = """
	<siteModel id="siteModel" gamma="true">
		<substitutionModel>
			<HKYModel idref="substmodel"/>
		</substitutionModel>
		<gammaShape gammaCategories="CATEGORIES">
			<parameter id="alpha" value="0.5" lower="0.0"/>
		</gammaShape>
	</siteModel>
"""
GAMMA_SITE_OPERATOR = """

		<scaleOperator scaleFactor="0.75" weight="1">
			<parameter idref="alpha"/>
		</scaleOperator>
"""
GAMMA_SITE_PRIOR = """
				<exponentialPrior mean="0.5" offset="0.0">
					<parameter idref="alpha"/>
				</exponentialPrior>
"""
GAMMA_SITE_LOG = '			<parameter idref="alpha"/>'

CONSTANT_COALESCENT = """
	<coalescentLikelihood id="coalescent">
		<model>
			<constantSize idref="constant"/>
		</model>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</coalescentLikelihood>
"""
CONSTANT_COALESCENT_PRIOR = """
				<oneOnXPrior>
					<parameter idref="constant.popSize"/>
				</oneOnXPrior>
				<coalescentLikelihood idref="coalescent"/>
"""
CONSTANT_COALESCENT_OPERATOR = """
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="constant.popSize"/>
		</scaleOperator>
"""
CONSTANT_COALESCENT_LOG = """
			<parameter idref="constant.popSize"/>
			<coalescentLikelihood idref="coalescent"/>
"""

EXPONENTIAL_COALESCENT = """
	<coalescentLikelihood id="coalescent">
		<model>
			<exponentialGrowth idref="exponential"/>
		</model>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</coalescentLikelihood>
"""
EXPONENTIAL_COALESCENT_PRIOR = """
				<oneOnXPrior>
					<parameter idref="exponential.popSize"/>
				</oneOnXPrior>
				<laplacePrior mean="0.0" scale="1.0">
					<parameter idref="exponential.growthRate"/>
				</laplacePrior>
				<coalescentLikelihood idref="coalescent"/>
"""
EXPONENTIAL_COALESCENT_OPERATOR = """
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="exponential.popSize"/>
		</scaleOperator>
		<randomWalkOperator windowSize="1.0" weight="3">
			<parameter idref="exponential.growthRate"/>
		</randomWalkOperator>
"""
EXPONENTIAL_COALESCENT_LOG = """
			<parameter idref="exponential.popSize"/>
			<parameter idref="exponential.growthRate"/>
			<coalescentLikelihood idref="coalescent"/>
"""

SKYRIDE = """
	<!-- Generate a gmrfSkyrideLikelihood for GMRF Bayesian Skyride process      -->
	<gmrfSkyrideLikelihood id="skyride" timeAwareSmoothing="true" randomizeTree="false">
		<populationSizes>

			<!-- skyride.logPopSize is in log units unlike other popSize                 -->
			<parameter id="skyride.logPopSize" dimension="NTAX" value="1.0"/>
		</populationSizes>
		<groupSizes>
			<parameter id="skyride.groupSize" dimension="NTAX"/>
		</groupSizes>
		<precisionParameter>
			<parameter id="skyride.precision" value="1.0" lower="0.0"/>
		</precisionParameter>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</gmrfSkyrideLikelihood>
"""
SKYRIDE_OPERATOR = """
		<gmrfBlockUpdateOperator scaleFactor="2.0" weight="2">
			<gmrfSkyrideLikelihood idref="skyride"/>
		</gmrfBlockUpdateOperator>
"""
SKYRIDE_PRIOR = """
				<gammaPrior shape="0.001" scale="1000.0" offset="0.0">
					<parameter idref="skyride.precision"/>
				</gammaPrior>
				<gmrfSkyrideLikelihood idref="skyride"/>
"""
SKYRIDE_LOG = """
			<gmrfSkyrideLikelihood idref="skyride"/>
			<parameter idref="skyride.precision"/>
			<parameter idref="skyride.logPopSize"/>
			<parameter idref="skyride.groupSize"/>

"""

SKYGRID = """
	<!-- Generate a gmrfSkyGridLikelihood for the Bayesian SkyGrid process       -->
	<gmrfSkyGridLikelihood id="skygrid">
		<populationSizes>

			<!-- skygrid.logPopSize is in log units unlike other popSize                 -->
			<parameter id="skygrid.logPopSize" dimension="POP_DIM" value="1.0"/>
		</populationSizes>
		<precisionParameter>
			<parameter id="skygrid.precision" value="0.1" lower="0.0"/>
		</precisionParameter>
		<numGridPoints>
			<parameter id="skygrid.numGridPoints" value="GRID_POINTS"/>
		</numGridPoints>
		<cutOff>
			<parameter id="skygrid.cutOff" value="GRID_CUTOFF"/>
		</cutOff>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</gmrfSkyGridLikelihood>
"""

SKYGRID_OPERATOR = """
		<gmrfGridBlockUpdateOperator scaleFactor="1.0" weight="2">
			<gmrfSkyrideLikelihood idref="skygrid"/>
		</gmrfGridBlockUpdateOperator>
		<scaleOperator scaleFactor="0.75" weight="1">
			<parameter idref="skygrid.precision"/>
		</scaleOperator>
"""
SKYGRID_PRIOR = """
				<gammaPrior shape="0.001" scale="1000.0" offset="0.0">
					<parameter idref="skygrid.precision"/>
				</gammaPrior>
				
				<gmrfSkyGridLikelihood idref="skygrid"/>
"""
SKYGRID_LOG = """
			<parameter idref="skygrid.precision"/>
			<parameter idref="skygrid.logPopSize"/>
			<parameter idref="skygrid.cutOff"/>
			<gmrfSkyGridLikelihood idref="skygrid"/>
"""

TOPOLOGY_FIXED = """
<newick id="startingTree">
    TREE_NEWICK
</newick>
"""

TOPOLOGY_SIMULATOR = """
<coalescentSimulator id="startingTree">
		<taxa idref="taxa"/>
		<constantSize idref="constant"/>
</coalescentSimulator>
"""
TOPOLOGY_OPERATOR = """
<subtreeSlide size="1.0" gaussian="true" weight="30">
	<treeModel idref="treeModel"/>
</subtreeSlide>
<narrowExchange weight="30">
	<treeModel idref="treeModel"/>
</narrowExchange>
<wideExchange weight="3">
	<treeModel idref="treeModel"/>
</wideExchange>
<wilsonBalding weight="3">
	<treeModel idref="treeModel"/>
</wilsonBalding>
"""

template = args.input.read()

ntax = int(re.findall(r"ntax=(\d+)", template)[0])

if args.categories == 1:
    template = (
        template.replace("SITE_MODEL", CONSTANT_SITE)
        .replace("SITE_PRIOR", "")
        .replace("SITE_OPERATOR", "")
        .replace("SITE_LOG", "")
    )
else:
    if args.site == "weibull":
        template = (
            template.replace(
                "SITE_MODEL", GAMMA_SITE.replace("CATEGORIES", str(args.categories))
            )
            .replace("SITE_PRIOR", GAMMA_SITE_PRIOR)
            .replace("SITE_OPERATOR", GAMMA_SITE_OPERATOR)
            .replace("SITE_LOG", GAMMA_SITE_LOG)
            .replace('gamma="true"', 'gamma="false"')
        )
    elif args.site == "invariant":
        template = (
            template.replace("SITE_MODEL", CONSTANT_SITE)
            .replace("SITE_PRIOR", "")
            .replace("SITE_OPERATOR", "")
            .replace("SITE_LOG", "")
        )
    else:
        template = (
            template.replace(
                "SITE_MODEL", GAMMA_SITE.replace("CATEGORIES", str(args.categories))
            )
            .replace("SITE_PRIOR", GAMMA_SITE_PRIOR)
            .replace("SITE_OPERATOR", GAMMA_SITE_OPERATOR)
            .replace("SITE_LOG", GAMMA_SITE_LOG)
        )


if args.model == "JC69":
    template = (
        template.replace("SUBSTITUTION_MODEL", HKY)
        .replace("SUBSTITUTION_PRIOR", "")
        .replace("SUBSTITUTION_OPERATOR", "")
        .replace("SUBSTITUTION_LOG", "")
    )
elif args.model == "HKY":
    template = (
        template.replace("SUBSTITUTION_MODEL", HKY)
        .replace("SUBSTITUTION_PRIOR", HKY_PRIOR)
        .replace("SUBSTITUTION_OPERATOR", HKY_OPERATOR)
        .replace("SUBSTITUTION_LOG", HKY_LOG)
    )
elif args.model == "GTR":
    template = (
        template.replace("SUBSTITUTION_MODEL", GTR)
        .replace("SUBSTITUTION_PRIOR", GTR_PRIOR)
        .replace("SUBSTITUTION_OPERATOR", GTR_OPERATOR)
        .replace("HKYModel", "gtrModel")
        .replace("SUBSTITUTION_LOG", GTR_LOG)
    )

if args.coalescent == "constant":
    template = (
        template.replace("COALESCENT_MODEL", CONSTANT_COALESCENT)
        .replace("COALESCENT_PRIOR", CONSTANT_COALESCENT_PRIOR)
        .replace("COALESCENT_OPERATOR", CONSTANT_COALESCENT_OPERATOR)
        .replace("COALESCENT_LOG", CONSTANT_COALESCENT_LOG)
    )
elif args.coalescent == "exponential":
    template = (
        template.replace(
            "COALESCENT_MODEL", EXPONENTIAL_COALESCENT.replace("NTAX", str(ntax - 1))
        )
        .replace("COALESCENT_PRIOR", EXPONENTIAL_COALESCENT_PRIOR)
        .replace("COALESCENT_OPERATOR", EXPONENTIAL_COALESCENT_OPERATOR)
        .replace("COALESCENT_LOG", EXPONENTIAL_COALESCENT_LOG)
    )
elif args.coalescent == "skyride":
    template = (
        template.replace("COALESCENT_MODEL", SKYRIDE.replace("NTAX", str(ntax - 1)))
        .replace("COALESCENT_PRIOR", SKYRIDE_PRIOR)
        .replace("COALESCENT_OPERATOR", SKYRIDE_OPERATOR)
        .replace("COALESCENT_LOG", SKYRIDE_LOG)
    )
elif args.coalescent == "skygrid":
    skygrid = (
        SKYGRID.replace("POP_DIM", str(args.grid))
        .replace("GRID_POINTS", str(args.grid - 1))
        .replace("GRID_CUTOFF", str(args.cutoff))
    )
    template = (
        template.replace("COALESCENT_MODEL", skygrid)
        .replace("COALESCENT_PRIOR", SKYGRID_PRIOR)
        .replace("COALESCENT_OPERATOR", SKYGRID_OPERATOR)
        .replace("COALESCENT_LOG", SKYGRID_LOG)
    )

if args.clock == "strict":
    template = (
        template.replace("CLOCK_MODEL", STRICT_CLOCK)
        .replace("CLOCK_PRIOR", STRICT_PRIOR)
        .replace("CLOCK_OPERATOR", STRICT_OPERATOR)
        .replace("CLOCK_LOG", STRICT_LOG)
        .replace("CLOCK_TREE_LOG", STRICT_LOG_TREE)
    )
elif args.clock == "ucln":
    template = (
        template.replace("CLOCK_MODEL", UCLN_CLOCK)
        .replace("CLOCK_PRIOR", "")
        .replace("CLOCK_OPERATOR", UCLN_OPERATOR)
        .replace("CLOCK_LOG", UCLN_LOG)
        .replace("CLOCK_TREE_LOG", UCLN_LOG_TREE)
    )

if args.topology:
    template = template.replace("TOPOLOGY_OPERATOR", "")
    tree = Tree.get(
        path=args.topology,
        schema="nexus",
        tree_offset=0,
        preserve_underscores=True,
        rooting="force-rooted",
    )
    newick = tree.as_string(unquoted_underscores=True, schema="newick")
    newick = newick[newick.index("(") :]
    template = template.replace(
        "STARTING_TOPOLOGY", TOPOLOGY_FIXED.replace("TREE_NEWICK", newick)
    )

else:
    template = template.replace("TOPOLOGY_OPERATOR", TOPOLOGY_OPERATOR).replace(
        "STARTING_TOPOLOGY", TOPOLOGY_SIMULATOR
    )

template = (
    template.replace("ITERATIONS", str(args.ngen))
    .replace("SAMPLE_FREQUENCY", str(args.samplefreq))
    .replace("STEM", args.stem)
)

args.output.write(template)
