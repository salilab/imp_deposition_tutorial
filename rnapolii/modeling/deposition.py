#!/usr/bin/env python3



# Imports needed to use ProtocolOutput
import IMP.pmi.mmcif
import ihm
import ihm.location
import ihm.reference
import ihm.model

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import ihm.cross_linkers

import os
import sys


#---------------------------
# Define Input Files
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology.txt"
target_gmm_file = datadirectory+'emd_1883.map.mrc.gmm.50.txt'

# Initialize model
m = IMP.Model()

# Read in the topology file.
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=datadirectory,
                                  fasta_dir=datadirectory,
                                  gmm_dir=datadirectory)

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)

# Record the modeling protocol to an mmCIF file
po = IMP.pmi.mmcif.ProtocolOutput()
bs.system.add_protocol_output(po)
po.system.title = "Modeling of RNA Pol II"
# Add publication
po.system.citations.append(ihm.Citation.from_pubmed_id(25161197))

# Each state can be specified by a topology file.
bs.add_state(topology)

# Build the system representation and degrees of freedom
root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
                                  max_rb_rot=0.3,
                                  max_bead_trans=4.0,
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)

# Fix all rigid bodies but not Rpb4 and Rpb7 (the stalk)
# First select and gather all particles to fix.
fixed_particles=[]
for prot in ["Rpb1","Rpb2","Rpb3","Rpb5","Rpb6","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12"]:
    fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot).get_selected_particles()

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,
                                          IMP.pmi.TransformMover])

# Randomize the initial configuration before sampling, of only the molecules
# we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

outputobjects = [] # reporter objects (for stat files)

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files

# Connectivity keeps things connected along the backbone (ignores if inside
# same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
        mol,scale=2.0, label=molname)
    cr.add_to_model()
    outputobjects.append(cr)

# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.add_to_model()
outputobjects.append(ev)


# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging components together)
xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("pep1.accession")
xldbkwc.set_protein2_key("pep2.accession")
xldbkwc.set_residue1_key("pep1.xlinked_aa")
xldbkwc.set_residue2_key("pep2.xlinked_aa")

xl1db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl1db.create_set_from_file(datadirectory+'polii_xlinks.csv')

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xl1db,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="Trnka",
                                   linker=ihm.cross_linkers.dss,
                                   weight=1.)

xl1.add_to_model()             # crosslink must be added to the model
outputobjects.append(xl1)


# Crosslinks - dataset 2
#  We can easily add a second set of crosslinks.
#  These have a different format and label, but other settings are the same
xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("prot1")
xldbkwc.set_protein2_key("prot2")
xldbkwc.set_residue1_key("res1")
xldbkwc.set_residue2_key("res2")

xl2db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl2db.create_set_from_file(datadirectory+'polii_juri.csv')

xl2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xl2db,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="Chen",
                                   linker=ihm.cross_linkers.bs3,
                                   weight=1.)
xl2.add_to_model()
outputobjects.append(xl2)


# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed because the EM restraint is quasi-Bayesian
em_components = IMP.pmi.tools.get_densities(root_hier)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                                 target_gmm_file,
                                                 scale_target_to_mass=True,
                                                 slope=0.000001,
                                                 weight=80.0)
gemt.add_to_model()
outputobjects.append(gemt)

#--------------------------
# Monte-Carlo Sampling
#--------------------------

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 20000
if '--test' in sys.argv: num_frames=100
num_mc_steps = 10

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange(m,
                                   root_hier=root_hier,
                                   monte_carlo_sample_objects=dof.get_movers(),
                                   output_objects=outputobjects,
                                   monte_carlo_temperature=1.0,
                                   simulated_annealing=True,
                                   simulated_annealing_minimum_temperature=1.0,
                                   simulated_annealing_maximum_temperature=2.5,
                                   simulated_annealing_minimum_temperature_nframes=200,
                                   simulated_annealing_maximum_temperature_nframes=20,
                                   replica_exchange_minimum_temperature=1.0,
                                   replica_exchange_maximum_temperature=2.5,
                                   number_of_best_scoring_models=10,
                                   monte_carlo_steps=num_mc_steps,
                                   number_of_frames=num_frames,
                                   global_output_directory="output",
                                   test_mode=True)

# Start Sampling
mc1.execute_macro()

po.finalize()

s = po.system

print("all subunits:",
      [a.details for a in s.asym_units])

print("first subunit sequence:",
      "".join(r.code for r in s.asym_units[0].entity.sequence))

print("all restraints on the system:", s.restraints)

import ihm.dumper
with open('initial.cif', 'w') as fh:
    ihm.dumper.write(fh, [s])

print("restraint datasets:", [r.dataset for r in s.restraints])

# Datasets for XL-MS restraint
for r in s.restraints:
    if isinstance(r, ihm.restraint.CrossLinkRestraint):
        print("XL-MS dataset at:", r.dataset.location.path)
        print("Details:", r.dataset.location.details)

# Dataset for EM restraint
em, = [r for r in s.restraints
       if isinstance(r, ihm.restraint.EM3DRestraint)]
d = em.dataset
print("GMM file at", d.location.path)
print("is derived from EMDB entry", d.parents[0].location.access_code)

# Definitions of some common crosslinkers
import ihm.cross_linkers

# Set linker types to match those used (DSS and BS3)
trnka, chen = [r for r in s.restraints
               if isinstance(r, ihm.restraint.CrossLinkRestraint)]
trnka.linker = ihm.cross_linkers.dss
chen.linker = ihm.cross_linkers.bs3

# Get last step of last protocol (protocol is an 'orphan' because
# we haven't used it for a model yet)
last_step = s.orphan_protocols[-1].steps[-1]
print(last_step.num_models_end)

# Correct number of output models to account for multiple runs
last_step.num_models_end = 200000

for subunit, accession in [('Rpb1', 'P04050'),
                           ('Rpb2', 'P08518')]:
    ref = ihm.reference.UniProtSequence.from_accession(accession)
    po.entities[subunit].references.append(ref)

# Get last protocol in the file
protocol = po.system.orphan_protocols[-1]
# State that we filtered the 200000 frames down to one cluster of
# 100 models:
analysis = ihm.analysis.Analysis()
protocol.analyses.append(analysis)
analysis.steps.append(ihm.analysis.ClusterStep(
                      feature='RMSD', num_models_begin=200000,
                      num_models_end=100))

mg = ihm.model.ModelGroup(name="Cluster 0")

# Add to last state
po.system.state_groups[-1][-1].append(mg)

e = ihm.model.Ensemble(model_group=mg,
                       num_models=100,
                       post_process=analysis.steps[-1],
                       name="Cluster 0")
po.system.ensembles.append(e)

import RMF
import IMP.rmf

# Add the model from RMF (let's say it's frame 42 from the output RMF file)
rh = RMF.open_rmf_file_read_only('output/rmfs/0.rmf3')
IMP.rmf.link_hierarchies(rh, [root_hier])
IMP.rmf.load_frame(rh, RMF.FrameID(42))
del rh
model = po.add_model(e.model_group)

# Look up the ihm.AsymUnit corresponding to a PMI component name
asym = po.asym_units['Rpb4.0']
# Add path to a local output file
loc = ihm.location.OutputFileLocation('output/cluster0.Rpb4.mrc')
den = ihm.model.LocalizationDensity(file=loc, asym_unit=asym)
# Add to ensemble
e.densities.append(den)

# Replace local links with DOIs
repo = ihm.location.Repository(doi="10.5281/zenodo.2598760", root="../..",
                  top_directory="salilab-imp_deposition_tutorial-1ad5919",
                  url="https://zenodo.org/record/2598760/files/salilab/"
                      "imp_deposition_tutorial-v0.2.zip")
s.update_locations_in_repositories([repo])

# Datasets for XL-MS restraint
trnka, chen = [r for r in s.restraints
               if isinstance(r, ihm.restraint.CrossLinkRestraint)]
d = trnka.dataset
print("XL-MS dataset now at %s/%s inside %s"
      % (d.location.repo.top_directory,
         d.location.path, d.location.repo.url))

po.finalize()
with open('rnapolii.cif', 'w') as fh:
    ihm.dumper.write(fh, [po.system])

#with open('rnapoliii.bcif', 'wb') as fh:
#    ihm.dumper.write(fh, [s], format='BCIF')

import ihm.reader
with open('rnapolii.cif') as fh:
    s, = ihm.reader.read(fh)
print(s.title, s.restraints, s.ensembles, s.state_groups)

import urllib.request

fh = urllib.request.urlopen('https://pdb-ihm.org/cif/8zze.cif')
s, = ihm.reader.read(fh)
print(s.title, s.restraints, s.ensembles, s.state_groups)
