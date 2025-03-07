from pymol import cmd
import os


def count_Hbond_GLU(structure,res1,atom1,atom1b,res2,atom2):
    struct = os.path.basename(structure)[:-4]
    cmd.load(structure,struct)
    cmd.remove("chain B")
    states = cmd.count_states(struct)
    hbond = 0
    for state in range(0,states):
        dist1 = cmd.get_distance(struct + " and resi " + res1 + " and name " + atom1, struct +" and resi " + res2     + " and name " + atom2,state)
        dist2 = cmd.get_distance(struct + " and resi " + res1 + " and name " + atom1b, struct +" and resi " + res2     + " and name " + atom2,state)
        dist = min(dist1,dist2)
        if dist <= 3.5:
            print(dist)
            hbond+=1
    cmd.delete('*')
    return [str(states),str(hbond),str(round((hbond/states)*100,1))]

Designed = count_Hbond_GLU("/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Designed/withTSA/ENSEMBLE_ref_tx_1_ptls_0.7_tbath_5/3NZ1.updated_ensemble.pdb",'178','OE1','OE2','262','N1')
Core = count_Hbond_GLU("/Volumes/Nia_HardDrive/ensemble_refinement_KEs/1A53-2/Core/withTSA/ENSEMBLE_ref_tx_1_ptls_0.6_tbath_5/1A53-core_TSA-coot-3_refine_001_refine_001.updated_ensemble.pdb",'178','OE1','OE2','401','N1')

with open('hbonds_1A53.csv','w') as infile:
    infile.write('Structure, No. of models, No. of Hbonds residue 178 to TSA, % of structures with hbond to resi 178'+'\n')
    infile.write('Designed'+ ',' + Designed[0] + ',' + Designed[1] + ',' + Designed[2]+'\n')
    infile.write('Core' + ',' + Core[0] + ',' + Core[1] + ',' + Core[2]+'\n')

