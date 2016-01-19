#!/usr/bin/env python

from matplotlib import rc
import matplotlib
from daft import PGM, Node, Plate
rc("font", family="serif", size=8)
rc("text", usetex=True)


pgm = PGM([9.5, 8.5], origin=[0., 0.2], observed_style='inner')

#pgm.add_node(Node('dispersion',r"\center{$\sigma_{Ia}$ \newline $\sigma_{!Ia}$}", 1,6,scale=1.2,aspect=1.8))
#pgm.add_node(Node('rate',r"{\center Relative \newline Rates}", 2,8,scale=1.6,aspect=1.2))
#pgm.add_node(Node('theta_T',r"\center{$\alpha_{Ia}$, $\alpha_{!Ia}$ \newline $\beta_{Ia}$, $\beta_{!Ia}$}", 4,6,scale=1.4,aspect=1.8))
pgm.add_node(Node('theta_T',r"\center{SNe~Ia, Non-Ia Populations \newline Rates}", 3,8,scale=1.8,aspect=3))
pgm.add_node(Node('Global Transmission',r"\centering{Global \newline Throughput}", 8, 8,scale=1.6,aspect=1.2))
pgm.add_node(Node('Transmission',r"Throughput", 8, 4,scale=1.6,aspect=1.2))
#pgm.add_node(Node('theta_T2',r"{Non-Ia}", 4,8,scale=1.6,aspect=1.2))
pgm.add_node(Node('mu',r"{\center Cosmology}", 7,8, scale=1.6,aspect=1.2))

#pgm.add_node(Node('G_Ni',r"$g_{Ni}$", 8,5))


#pgm.add_node(Node('Flux_g',r"$f_{gi}(\lambda)$", 8, 3))

#pgm.add_node(Node('z',r"$z_i$", 3, 3, fixed=True,offset=(-10,-5)))


#pgm.add_node(Node('Counts',r"$\overline{n}_i$", 7, 2,scale=1.2,fixed=True,offset=(-15,0)))
#pgm.add_node(Node('Counts_g',r"$\overline{\mathit{ADU}}_{gi}$", 8, 2,scale=1.2,fixed=True,offset=(20,0)))
#pgm.add_node(Node("Summary",r"$S({\mathit{ADU}_i})$", 9, 1, fixed=True, offset=(0,-25)))

pgm.add_node(Node('G_i',r"{\centering Host Galaxy \newline (Redshift)}", 3,7, scale=1.6,aspect=1.2))


pgm.add_node(Node('Type',r"\centering{Type \newline Subtype}", 2, 6, scale=1.6,aspect=1.2))


#pgm.add_node(Node('theta_Ti',r"\centering{Transient \newline Parameters}", 1,5,scale=1.5,aspect=1.2))
pgm.add_node(Node('HD',r"Distance", 7,6,fixed=True,offset=(-25,-15)))


pgm.add_node(Node('Luminosity',r"Luminosity", 5, 4,scale=1.6,aspect=1.2))


pgm.add_node(Node('Flux',r"Flux", 7, 3, scale=1.2,fixed=True,offset=(-20,-20)))




#pgm.add_node(Node('Stheta',r"Indicators", 1, 1,  observed=True,scale=1.6,aspect=1.2))
pgm.add_node(Node('ST',r"\centering{Type \newline Subtype$_o$}", 2, 1, observed=True,scale=1.6,aspect=1.2))
pgm.add_node(Node('Sz',r"\centering{Redshift$_o$\newline Phot \& Spec}", 3, 1, observed=True,scale=1.6,aspect=1.2))
pgm.add_node(Node('^Counts',r"Photometry$_o$", 8, 1, observed=True,scale=1.6,aspect=1.2))
pgm.add_node(Node('Gals',r"\centering{Host Properties$_o$\newline Phot \& Spec}", 4,1, observed=True,scale=1.6,aspect=1.2))


#pgm.add_node(Node('Soectroscopy',r"Transient Spectroscpy", 1, 1,  observed=True))


#pgm.add_edge("^Counts","Summary")
pgm.add_edge("G_i","Gals")
#pgm.add_edge("G_Ni","Gals")
pgm.add_edge("theta_T","Type")
#pgm.add_edge("Type", "theta_Ti")
pgm.add_edge("mu","HD")
#pgm.add_edge("theta_T", "theta_Ti")
pgm.add_edge("G_i","HD")
pgm.add_edge("G_i","Sz")
pgm.add_edge("G_i","Luminosity")
pgm.add_edge("G_i","Flux")
pgm.add_edge("HD","Flux")
#pgm.add_edge("z","Flux")
#pgm.add_edge("G_Ni","Flux_g")
#pgm.add_edge("HD","Flux_g")
#pgm.add_edge("z","Flux_g")
pgm.add_edge("G_i","Type")
#pgm.add_edge("G_i","theta_Ti")
pgm.add_edge("Type","Luminosity")
pgm.add_edge("Luminosity","Flux")
#pgm.add_edge("Type","Stheta")
pgm.add_edge("theta_T","Luminosity")
#pgm.add_edge("theta_Ti","Luminosity")
pgm.add_edge("Flux","^Counts")
#pgm.add_edge("Flux_g","Counts_g")
#pgm.add_edge("Counts","^Counts")
#pgm.add_edge("Counts_g","^Counts")
pgm.add_edge("Global Transmission","Transmission")
pgm.add_edge("Transmission","^Counts")
#pgm.add_edge("Transmission","Counts_g")
pgm.add_edge("Type","ST")
#pgm.add_edge("theta_Ti","Stheta")
#pgm.add_edge("Counts_g","^Counts")

# Big Plate: Galaxy
pgm.add_plate(Plate([1.4, 0.5, 7.2, 7.],
                    label=r"SNe $i = 1, \cdots, N_{SN}$",
                    shift=-0.2,label_offset=[20,2]))


pgm.add_plate(Plate([4.5, 0.6, 4, 4.],
                    label=r"{\centering LC Point $j = 1, \cdots, N_{band}$, $k = 1, \cdots, N_{date}$}",
                    shift=-0.2,label_offset=[2,2]))


# Render and save.
pgm.render()
pgm.figure.text(0.01,0.9,r'\underline{UNIVERSAL}',size='large')
pgm.figure.text(0.01,0.55,r'{\centering \underline{INDIVIDUAL} \newline \underline{SN}}',size='large')
pgm.figure.text(0.01,0.2,r'\underline{OBSERVATORY}',size='large')
pgm.figure.text(0.01,0.1,r'\underline{DATA}',size='large')


pgm.figure.savefig("../results/hdpgm.pdf")
