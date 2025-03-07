# Rojo's script for creating a figure containing the 
# electron density with PyMol 
# usage: 
# open PyMol, make sure you are in the directory containing
# the script (do pwd in the PyMol command line window and cd to your directory) 
# in the command line type: run create_map_for_figure.py


from pymol import cmd, stored, math 

# List of the structures 
struct=["HG3_TSA","3NYD"]

# open PDB, make sure you put the correct path where your PDBs are
for pdb in struct: 
    cmd.load("../Xray/structures_for_review/"+pdb+".pdb", pdb)

# open Map file, here the format in MTZ
for pdb in struct: 
    cmd.load_mtz("../Xray/structures_for_review/"+pdb+".mtz", pdb+".map")

    # --- prepare PDB

    # Remove H and chains that you don't need 
    cmd.remove('chain B')
    cmd.remove('h.')

    # If you are comparing 2 structures with different unit cell
    # you need to align the structures first (here 3NYD aligned to HG3_TSA)
    # and then you copy 3NYD matrix to the new coordinates of 3NYD
    cmd.extra_fit()
    cmd.align("3NYD", "HG3_TSA") 
    cmd.matrix_copy("3NYD", "3NYD.map.2fofc") 

    # Use this command for coloring all C atoms to grey80
    cmd.color("grey80", " n. C*") 

    # --- prepare map 

    # Since in this script the density map will be represented 
    # in volume, you need to specify a range of color that will match 
    # the sigma of the density (for more info go to pymolwiki) 
    cmd.volume_ramp_new('ramp537', [\
      0.51, 0.00, 1.00, 1.00, 0.00, \
      0.60, 0.33, 1.00, 1.00, 0.17, \
      0.65, 0.33, 1.00, 1.00, 0.00, \
      1.32, 0.00, 0.00, 1.00, 0.00, \
      1.46, 0.00, 0.00, 1.00, 0.06, \
      1.62, 0.00, 0.00, 1.00, 0.00, \
      1.72, 0.00, 0.00, 1.00, 0.00, \
	  ])

    # to create the map, first select the residues of interest 
    cmd.select("site_"+pdb, pdb+" and chain A")

    # Create the map: here we are using the volume representation,
    # color the volume with coloration option we set up previously
    # if you want to visualize in mesh, you can use isomesh instead 
    # the command is similar (more info about isomesh: pymolwiki)
    cmd.volume("site_map_"+pdb, pdb+".map.2fofc", "ramp537", "site_"+pdb,"0.0","1","1.6", "0","1")

    # if you have many selections to analyse, you can group them with the 
    # group command  
    cmd.group(pdb+"group", " site_"+pdb+" site_map_"+pdb)

    # ---- set up the view

    # if you are using isomesh, a width of 0.5 is better
    cmd.set("mesh_width","0.5")

    # a simple command to color mesh 
    cmd.set("mesh_color","lightblue")

    # Hide the representations that you don't need
    cmd.hide("sticks") 
    cmd.hide("spheres") 
    cmd.hide("cartoon")

    # a command to remove double bonds from the representation
    cmd.set("valence", "0")

    # For a better visualization, remove the shadow 
    cmd.set("ray_shadow", "0")

# In your pymol GUI, set up the orientation needed and render
# the image 
# Final note: unfortunately, it is not yet possible to render 
# images with "ray" when you are using the volume representation 
# Instead, use: draw 4000, 3000, antialias=2
# png highres.png
