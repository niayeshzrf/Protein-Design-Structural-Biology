bg_color white
set cartoon_transparency,0.6
set cartoon_fancy_sheets, 1
set cartoon_fancy_helices, 1
set ambient, 0.5
set ray_trace_mode, 1
set cartoon_gap_cutoff, 0
set cartoon_side_chain_helper,1
set ray_shadow,0
set valence,0

directory='../ensembles/'

load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/ENSEMBLE_ref_tx_1_ptls_0.6_tbath_5/1A53-core_TSA-coot-3_refine_001_refine_001.updated_ensemble.pdb", 1A53_C_TSA_ensemble_1
#load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/ENSEMBLE_ref_tx_1_ptls_0.7_tbath_5/1A53-core_TSA-coot-3_refine_001_refine_001.updated_ensemble.pdb", 1A53_C_TSA_ensemble_2
#load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/ENSEMBLE_ref_tx_1_ptls_0.9_tbath_5/1A53-core_TSA-coot-3_refine_001_refine_001.updated_ensemble.pdb", 1A53_C_TSA_ensemble_3
#load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/ENSEMBLE_ref_tx_1_ptls_1.0_tbath_5/1A53-core_TSA-coot-3_refine_001_refine_001.updated_ensemble.pdb", 1A53_C_TSA_ensemble_4


#load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Designed/withTSA/ENSEMBLE_ref_tx_1_ptls_0.6_tbath_5/3NZ1.updated_ensemble.pdb", 1A53_D_TSA_ensemble_1
load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Designed/withTSA/ENSEMBLE_ref_tx_1_ptls_0.7_tbath_5/3NZ1.updated_ensemble.pdb", 1A53_D_TSA_ensemble_2
#load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Designed/withTSA/ENSEMBLE_ref_tx_1_ptls_0.9_tbath_5/3NZ1.updated_ensemble.pdb", 1A53_D_TSA_ensemble_3
#load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Designed/withTSA/ENSEMBLE_ref_tx_1_ptls_1.0_tbath_5/3NZ1.updated_ensemble.pdb", 1A53_D_TSA_ensemble_4


fetch 1A53, 1A53


remove chain B+C+D+E+F+G
remove solvent
extra_fit *,1A53

color palegreen, name C*
color palegreen, *_ensemble and name C*
color wheat, resn 6NT and elem C
color wheat, resn 3NY and elem C
hide everything
select active_site, resi 178+210+110+ resn 3NY
show spheres, active_site and name ca
set sphere_scale,0.1
set line_width, 2
show lines, active_site
hide lines, hydrogens
hide lines, name c+n+o
zoom active_site

set_view (\
     0.381290466,   -0.181602970,   -0.906443179,\
     0.177895308,   -0.947776616,    0.264713943,\
    -0.907175124,   -0.262185067,   -0.329071313,\
     0.000121589,    0.000052478,  -53.814769745,\
    -2.961925507,  -18.382286072,   34.152565002,\
    42.426811218,   65.199783325,  -20.000000000 )



set state,0



disable all
enable 1A53_D_TSA_ensemble_1
scene 1A53_D_TSA, store




disable all
enable 1A53_C_TSA_ensemble_1
scene 1A53_C_TSA, store

