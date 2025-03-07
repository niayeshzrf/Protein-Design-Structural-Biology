from pymol import cgo
from pymol import cmd

# Unit sphere at the origin
#plain_sphere = [cgo.SPHERE, 0, 0, 0, 1]
#cmd.load_cgo(plain_sphere, "plain_sphere")

# Red unit sphere at (3, 4, 5)
#red_sphere = [cgo.COLOR, 1, 0, 0, cgo.SPHERE, 3, 4, 5, 1]
#cmd.load_cgo(red_sphere, "red_sphere")


cmd.fetch("3NZ1")

cmd.pseudoatom("PS1","resn 3NY")

xyz = cmd.get_coords("name PS1")[0]
Sphere = [cgo.COLOR,1,0.8,0.5, cgo.SPHERE, xyz[0], xyz[1], xyz[2], 4]
cmd.load_cgo(Sphere, "plain_sphere")

Sphere2 = [cgo.COLOR,0.4,0.1,0.7, cgo.SPHERE, xyz[0], xyz[1], xyz[2], 11]
cmd.load_cgo(Sphere2, "plain_sphere2")

cmd.remove("solvent")
cmd.remove("inorganic")
cmd.remove("resn TLA")
cmd.set("cgo_transparency", 0.3 ,"plain_sphere")
cmd.set("cgo_transparency", 0.5 ,"plain_sphere2")
cmd.set("ray_trace_mode",1)
cmd.set("ambient",0.8)
#cmd.set("bg_color","white")
cmd.color("white","elem C")
cmd.color("orange","resn 3NY and elem C")
cmd.run("/Users/Nia/Documents/PymolWD/Circle.py")
cmd.circleSelection("resn 3NY",4)

