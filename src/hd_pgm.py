#!/usr/bin/env python

from matplotlib import rc
from daft import PGM, Node, Plate
rc("font", family="serif", size=10)
rc("text", usetex=True)


pgm = PGM([10.5, 7.5], origin=[0., 0.2], observed_style='inner')

pgm.add_node(Node('G',r"$G$", 3,1))
pgm.add_node(Node('Coords',r"$\mathit{RA/Dec}_i$", 2,2,scale=1.5,aspect=1.2,fixed=True))
pgm.add_node(Node('G_i',r"$g_i$", 3,2))

pgm.add_node(Node('mu',r"$\theta_\mu$", 1,4))
pgm.add_node(Node('HD',r"$\mu_i$", 4,4))
pgm.add_node(Node('theta_T',r"$\theta_T^{Ia}$, $\theta_T^{non-Ia}$", 3,7,scale=1.5,aspect=1.5))
pgm.add_node(Node('theta_Ti',r"$\theta_{Ti}^{Ia}$, $\theta_{Ti}^{non-Ia}$", 2,6,scale=1.5,aspect=1.5))


#pgm.add_node(Node('theta_G',r"$\theta_G$", 2,1))



#pgm.add_node(Node('theta_Gi',r"$\theta_{Gi}$", 3,4))


pgm.add_node(Node('Host',r"$\theta_{gi}$", 2, 3))
pgm.add_node(Node('z',r"$z_i$", 3, 3))


pgm.add_node(Node('Type',r"$T_i$", 2, 5))
pgm.add_node(Node('Luminosity',r"$L_i(t,\lambda)$", 4, 6, scale=1.2))
pgm.add_node(Node('Flux',r"$n_i(t,\lambda)$", 4, 5, scale=1.2))
pgm.add_node(Node('Flux_g',r"$n_{gi}(\lambda)$", 4, 2))
pgm.add_node(Node('Transmission',r"$\phi(\lambda)$", 5, 7))
pgm.add_node(Node('Counts',r"$\mathit{ADU}_i$", 6, 5,scale=1.2))
pgm.add_node(Node('Counts_g',r"$\mathit{ADU}_{gi}$", 6, 2,scale=1.2))
pgm.add_node(Node('Zeropoints',r"$\hat{Z}$", 7, 7, observed=True))
pgm.add_node(Node('Counts_V',r"$\mathit{ADU}_{Vi}$", 6, 4, scale=1.3))
pgm.add_node(Node('^Counts',r"$\hat{\mathit{ADU}_i}$,$\hat{\mathit{ADU}_{gi}}$", 8, 4, observed=True,scale=2,aspect=1.2))

pgm.add_node(Node('Spars',r"$\hat{T}_{Si}$, $\hat{z}_{Si},\hat{\theta}_{Si}$", 8, 6, scale=1.8,aspect=1.2, observed=True))
#pgm.add_node(Node('^Host',r"$\hat{z}_{Hi},\hat{\theta}_{Hi}$", 7, 2, scale=1.5,observed=True))
pgm.add_node(Node('Galaxies',r"$\hat{\mathit{Gals.}}$", 8, 1, observed=True))

pgm.add_node(Node('Detected',r"$\hat{D}_i$", 7, 5, observed=True))


pgm.add_node(Node('^Type',r"$\hat{T}_i$", 9, 4, observed=True))

pgm.add_edge("mu","HD")
pgm.add_edge("Spars","^Type")
#pgm.add_edge("^Host","^Type")
pgm.add_edge("^Counts","^Type")

pgm.add_edge("Host","Luminosity")


pgm.add_edge("z","HD")


pgm.add_edge("Coords","G_i")
pgm.add_edge("G_i","Host")

pgm.add_edge("G_i","z")
pgm.add_edge("HD","Flux")
pgm.add_edge("z","Flux")
pgm.add_edge("HD","Flux_g")
pgm.add_edge("z","Flux_g")
pgm.add_edge("^Counts","Galaxies")
#pgm.add_edge("HD","^Host")


pgm.add_edge("Host","Type")

pgm.add_edge("Type","Luminosity")
pgm.add_edge("Luminosity","Flux")
pgm.add_edge("G","G_i")


pgm.add_edge("theta_T","Luminosity")
#pgm.add_edge("theta_G","Flux_g")
pgm.add_edge("G","Galaxies")

pgm.add_edge("theta_Ti","Luminosity")
pgm.add_edge("Host","Flux_g")


pgm.add_edge("Flux","Counts")
pgm.add_edge("Flux_g","Counts_g")

pgm.add_edge("Counts","Counts_V")
pgm.add_edge("Counts_V","^Counts")
pgm.add_edge("Counts_g","^Counts")


pgm.add_edge("Transmission","Counts")
pgm.add_edge("Transmission","Counts_g")
pgm.add_edge("Transmission","Zeropoints")

pgm.add_edge("Flux","Spars")

#pgm.add_edge("Galaxies","^Host")

pgm.add_edge("Counts_g","Counts_V")

pgm.add_edge("Counts_V","Detected")
pgm.add_edge("Detected","^Counts")
pgm.add_edge("Detected","Spars")
#pgm.add_edge("Detected","^Host")




# pgm.add_edge("Galaxies", "gal_zi")


# Per-Galaxy parameters: third line in the plate


# Observed photometry

# pgm.add_edge("gal_ni","^gal_zi")
# pgm.add_edge("gal_ni","^gal_thetai")
# pgm.add_edge("gal_ni","^gal_ni")


# pgm.add_edge("gal_D","^gal_zi")
# pgm.add_edge("gal_D","^gal_thetai")
# pgm.add_edge("gal_D","^gal_ni")

# pgm.add_edge("^ni","^Type")


# pgm.add_edge("Standard","Z")


#pgm.add_node(Node("t0true", r"$t_0^{\mathrm{true}}$", 7, 1))


# Big Plate: Galaxy
pgm.add_plate(Plate([1.5, 1.5, 7, 5.],
                    label=r"SNe $\i = 1, \cdots, N_{SN}$",
                    shift=-0.1))

# # Big Plate: SNe
# pgm.add_plate(Plate([5.5, 0.5, 4, 6.],
#                     label=r"SNe $i = 1, \cdots, N_{SN}$",
#                     shift=-0.1))


# pgm.add_plate(Plate([2.6, 1.6, 1.8, 3.8],
#                     label=r"Photo-$z$",
#                     shift=-0.1))

# pgm.add_plate(Plate([3.6, 5.6, 2.8, 1.8],
#                     label=r"Host match",
#                     shift=-0.1))

# pgm.add_plate(Plate([5.6, 1.6, 1.8, 3.8],
#                     label=r"LC Fit",
#                     shift=-0.1))

# pgm.add_plate(Plate([5.6, 2.6, 2.8, 0.8],
#                     label=r"Efficiency",
#                     shift=-0.1))


# Render and save.
pgm.render()
pgm.figure.savefig("../results/hdpgm.eps")
