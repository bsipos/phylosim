#!/usr/bin/env Rscript

##
## Example X7: Simulating amino acid sequences under a site-hetereogeneous model.
##

# load PhyloSim
library("phylosim");

# Enable the "fast & careless mode":
PSIM_FAST <- TRUE

lg<-LG();	# Create a LG substitution process.

# Specify root sequence length:
root.len<-1000

#
# Setting up the substitution processes:
#

# Create LG substitution processes with different equilibrium frequencies:
nr.processes<- root.len

processes<-list()
for (i in 1:nr.processes) {
    tmp <- clone(lg)
    tmp$name <- paste("root site", i)
    tmp$equDist<-abs(rnorm(20)) # One could set the equilibrium distribution here based on some other logic.
    processes<-c(processes,list(tmp))

}

# Create a continous deletor process:
cont.del<-ContinuousDeletor(
			rate=0.05,	# global rate for this deletion process
			max.length=10,	# the maximum allowed deletion length
			dist=expression(rnorm(1,mean=5,sd=3))	# length sampling expression
		);

# Set up a template sequence for the insertion process:
# We attach the substitution processes in the insert hook, so the site states
# must be set. We will clear them in the insert hook:
templ.seq.lg<-AminoAcidSequence(string="CCCCCCCCCC"); # this is just a sequence with length 10.

# Note that the template sequence state is undefined, so the states
# will be sampled from the equlibrium distribution of the substitution process(es). 

# Create a continous insertor process object:
cont.ins.lg<-ContinuousInsertor(
			rate=0.05,	# global rate for this insertion process
			max.length=10,	# the maximum allowed insertion length
			dist=expression(rnorm(1,mean=5,sd=3)) # length sampling expression
		);

# Setting up the template sequence for the insertion process:

# Attach insertion and deletion processes to the template sequence:
attachProcess(templ.seq.lg, cont.del)
attachProcess(templ.seq.lg, cont.ins.lg)

#print(templ.seq.lg$processes)

# Disabling write protection for the insertion processes:
cont.ins.lg$writeProtected<-FALSE;

# Setting the template sequence for the insertion processes:
cont.ins.lg$templateSeq<-templ.seq.lg;


# Setting up the insert hook for the insertion processes:
cont.ins.lg$insertHook<-function(seq,target.seq,event.pos,insert.pos){
        # Clear sequence states:
        clearStates(seq)
	# Sample a substitution process for each site from the list of
	# available processes:
	for (i in 1:seq$length) {
		p<-sample(processes,1)[[1]]
		attachProcess(seq$sites[[i]], p)
	}
        # Sample states from the subsitution processes:
        sampleStates(seq)
	return(seq)
}

#
# Setting up the root sequence:
#

seq<-AminoAcidSequence(length=root.len); # Create a sequence of specified length.

# Attach the substitution processes:
for(i in 1:seq$length) {
    attachProcess(seq$sites[[i]], processes[[i]])
}
# print(seq$processes)


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
