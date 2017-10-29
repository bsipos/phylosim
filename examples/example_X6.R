#!/usr/bin/env Rscript

##
## Example X6: Simulating amino acid sequences under a mixture model.
##

# load PhyloSim
library("phylosim");

# Enable the "fast & careless mode":
PSIM_FAST <- TRUE

lg<-LG();	# Create a LG substitution process.

#
# Setting up the substitution processes:
#

# Create LG substitution processes with different equilibrium frequencies:
nr.processes<- 5

processes<-list()
for (i in 1:nr.processes) {
    tmp <- clone(lg)
    tmp$equDist<-abs(rnorm(20)) # Sample the equilibrium distribution.
    processes<-c(processes,list(tmp))
    # plot(tmp)
}

# Create a vector of relative rates:
rel.rates <- numeric(nr.processes)
for (i in 1:nr.processes) {
    rel.rates[i]<-abs(rnorm(1))  # Here the rates can be set in a more complicated way.  
}
rel.rates<-rel.rates/sum(rel.rates)

# Create a continous deletor process:
cont.del<-ContinuousDeletor(
			rate=0.01,	# global rate for this deletion process
			max.length=10,	# the maximum allowed deletion length
			dist=expression(rnorm(1,mean=5,sd=3))	# length sampling expression
		);

templ.seq.lg<-AminoAcidSequence(length=10); # this is just a sequence with length 10.

# Note that the template sequence state is undefined, so the states
# will be sampled from the equlibrium distribution of the substitution process(es). 

# Create a continous insertor process object:
cont.ins.lg<-ContinuousInsertor(
			rate=0.01,	# global rate for this insertion process
			max.length=10,	# the maximum allowed insertion length
			dist=expression(rnorm(1,mean=5,sd=3)) # length sampling expression
		);

# Setting up the template sequence for the insertion process:
setProcesses(templ.seq.lg, list(processes));

# Set mixture proportions for the template sequence:
for (i in 1:nr.processes){
    setRateMultipliers(templ.seq.lg, processes[[i]], rel.rates[i])
}

# Attach insertion and deletion processes to the template sequence:
attachProcess(templ.seq.lg, cont.del)
attachProcess(templ.seq.lg, cont.ins.lg)

#print(templ.seq.lg$processes)

# Disabling write protection for the insertion processes:
cont.ins.lg$writeProtected<-FALSE;

# Setting the template sequence for the insertion processes:
cont.ins.lg$templateSeq<-templ.seq.lg;

#
# Setting up the root sequence:
#

seq<-AminoAcidSequence(length=200); # Create a sequence of length 200.

# Attach the substitution processes:
setProcesses(seq, list(processes));

# Set mixture proportions for the root sequence:
for (i in 1:nr.processes){
    setRateMultipliers(seq, processes[[i]], rel.rates[i])
}

# Attach deletion and insertion processes to the root sequence:
attachProcess(seq, cont.del)
attachProcess(seq, cont.ins.lg)

# Sample the states from the attached substitution process(es):
sampleStates(seq);

print(seq); # Print the actual sequence.

plot(seq);  # Plot the "rate landscape".

# Read in the tree using APE:
tree<-read.tree(
		file="data/smalldemotree.nwk"	# the path to the tree file
	);

# Create the simulation object:
sim<-PhyloSim(
			phylo=tree,	# the tree as an APE phylo object
			root.seq=seq	# the root sequence.
		);

# Run the simulation:
Simulate(sim)

# Plot the resulting alingment alongside the tree:
plot(sim)

# Save the resulting alignment, skip internal nodes:
saveAlignment(
		sim,				# the phylo object	
		file="example_X6_aln.fas",
		skip.internal=FALSE		# filename for alignment
);

# Disable fast mode:
rm(PSIM_FAST)
