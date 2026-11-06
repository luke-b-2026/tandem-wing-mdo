#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================
from base import get_mass
import numpy as np


if True: #proteus geometry
    fwng_xpos0	=-6
    fwng_zpos0	=0
    fwng_chrd0	=1.32
    fwng_inci0	=0
    fwng_span0	=3.78
    fwng_dihe0	=5
    fwng_swep0	=0
    fwng_chrd1	=1.132
    fwng_inci1	=0
    fwng_span1	=2.63
    fwng_dihe1	=5
    fwng_swep1	=0
    fwng_chrd2	=0.943
    fwng_inci2	=0
    fwng_span2	=2.08
    fwng_dihe2	=5
    fwng_swep2	=0
    fwng_chrd3	=0.755
    fwng_inci3	=0
    rwng_xpos0	= 1
    rwng_zpos0	=0
    rwng_chrd0	=1.51
    rwng_inci0	=0
    rwng_span0	=3.9
    rwng_dihe0	=0
    rwng_swep0	=0
    rwng_chrd1	=1.51
    rwng_inci1	=0
    rwng_span1	=4.63
    rwng_dihe1	=8
    rwng_swep1	=0
    rwng_chrd2	=0.943
    rwng_inci2	=0
    rwng_span2	=3.45
    rwng_dihe2	=-8
    rwng_swep2	=0
    rwng_chrd3	=0.755
    rwng_inci3	= 0
    tail_xpos	= 8
    tail_zpos	= -2
    tail_span	= 4.15
    tail_chrd0	= 1.89
    tail_chrd1  = 1.8
    tail_swep   = 0

geom = [fwng_xpos0,	fwng_zpos0,	fwng_chrd0,	fwng_inci0,	fwng_span0,	fwng_dihe0,	fwng_swep0,	fwng_chrd1,	fwng_inci1,	fwng_span1,	fwng_dihe1,	fwng_swep1,	fwng_chrd2,	fwng_inci2,	fwng_span2,	fwng_dihe2,	fwng_swep2,	fwng_chrd3,	fwng_inci3,	rwng_xpos0,	rwng_zpos0,	rwng_chrd0,	rwng_inci0,	rwng_span0,	rwng_dihe0,	rwng_swep0,	rwng_chrd1,	rwng_inci1,	rwng_span1,	rwng_dihe1,	rwng_swep1,	rwng_chrd2,	rwng_inci2,	rwng_span2,	rwng_dihe2,	rwng_swep2,	rwng_chrd3,	rwng_inci3,	tail_xpos,	tail_zpos,	tail_span,	tail_chrd0, tail_chrd1, tail_swep]
#proteus geom


mass_payload = 1000 #kg
aircraft_range = 3000*1e3 #m
aircraft_endurance = 20*3600 #s

mass, stability = get_mass(geom, mass_payload, aircraft_range, aircraft_endurance)
print(mass, stability)
#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================