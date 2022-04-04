import yt

ds = yt.load("BR2D_stripes_0040.dat")

yt.SlicePlot(ds, "x", ("amrvac", "v"), width=(800.0,"kpc")).save()
