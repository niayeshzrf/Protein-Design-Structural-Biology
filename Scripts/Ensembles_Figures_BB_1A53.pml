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


load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/ENSEMBLE_ref_tx_1_ptls_0.6_tbath_5/1A53-core_TSA-coot-3_refine_001_refine_001.updated_ensemble.pdb", 1A53_C_TSA_ensemble
load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Designed/withTSA/ENSEMBLE_ref_tx_1_ptls_0.7_tbath_5/3NZ1.updated_ensemble.pdb", 1A53_D_TSA_ensemble


fetch 3NZ1, 1A53_D_TSA_Xtal
load "/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/1A53-core_TSA-coot-3_refine_001_refine_001.pdb"
fetch 1A53, 1A53


remove chain B+C+D+E+F+G
remove solvent
extra_fit *,1A53


color gray60, *_ensemble and name C*

color warmpink, *Xtal

hide everything
show ribbon
orientD

set_view (\
     0.907538593,   -0.419851094,    0.009940351,\
     0.017975040,    0.062480159,    0.997884333,\
    -0.419583887,   -0.905439854,    0.064249992,\
     0.000000000,    0.000000000, -173.359161377,\
    15.930103302,    3.870658875,   15.548716545,\
   136.677673340,  210.040649414,  -20.000000000 )

set state,0





disable all
enable 1A53_D_TSA_ensemble
enable 1A53_D_TSA_Xtal
scene 1A53_D_TSA, store




disable all
enable 1A53_C_TSA_ensemble
enable 1A53_C_TSA_Xtal
scene 1A53_C_TSA, store

