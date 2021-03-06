#!/usr/bin/env Rscript

##
## Example V3.2: Evolving codon sequences
## See also the package vignette (vignette("PhyloSim",package="phylosim")).
##

library(phylosim)

# Enable "fast & careless" mode:
PSIM_FAST<-TRUE;

# Construct a GY94 codon substitution model:
p<-GY94();

# Set the transition/transverion rate ratio:
p$kappa=2

# Sample codon frequencies from a normal distribution:
p$equDist<-abs(rnorm(61,mean=10,sd=3))

# Get object summary for p:
summary(p)

# Get a bubble plot of p:
plot(p,scale=0.8)

# Construct a discrete deletor process:
d<-DiscreteDeletor(
	rate=1,
	sizes=1:4,
	probs=c(4,3,2,1)/10
);

# Construct a discrete insertor process inserting neutrally evolving sites:
i<-DiscreteInsertor(
	rate=1.5,
	sizes=1:4,
	probs=c(4,3,2,1)/10,
	template.seq=CodonSequence(length=4,processes=list(list(p)))
);

# Construct root sequence and attach process p:

s<-CodonSequence(length=30,processes=list(list(p)))

# Sample omegas from a discrete model:
omegaVarM3(s,p,omegas=c(0,1,2),probs=c(2/4,1/4,1/4))

# Plot the omega values across sites:
plotParametersAtSites(s,p,"omega");

# Sample states:

sampleStates(s)

# Construct the simulation object:
sim<-PhyloSim(
	root.seq=s,
	phylo=read.tree("data/smalldemotree.nwk")
);

# Create a node hook function and attach to node 9:
node.hook<-function(seq){

	# Set all omegas to 1 (neutral):
	setOmegas(seq,p,1);
	# attach the deletion process:
	attachProcess(seq,d)
	# attach the insertion process:
	attachProcess(seq,i)

        return(seq);
}

attachHookToNode(
                sim,                    # PhyloSim object.
                node=9,                 # the node
                fun=node.hook           # the node hook function
);

# Disable fast mode just before simulation in order to preserve branch statistics:
rm(PSIM_FAST)

# Run the simulation:
Simulate(sim)

# Plot the resulting alingment alongside the tree:
plot(sim);

# Export the nonsynonymous substitution counts as a phylo object:
nsyn.subst<-exportStatTree(sim,"nr.nsyn.subst")

# Plot the exported phylo object:
plot(nsyn.subst)
nodelabels()

# Save the resulting alignment:
saveAlignment(
                sim,                            # the phylo object      
                file="example_V3.2_aln.fas",      # filename for alignment
);

