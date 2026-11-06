#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================

import json
from math import radians, sqrt, tan
import avlwrapper as avl
from avlwrapper import default_config
import numpy as np
import matplotlib.pyplot as plt

def get_mass(geom, mass_payload, aircraft_range, aircraft_endurance):

    rho = 1.225 #kg/m3

    # print(geom)
    fwng_xpos0	=	geom	[	0	]
    fwng_zpos0	=	geom	[	1	]
    fwng_chrd0	=	geom	[	2	]
    fwng_inci0	=	geom	[	3	]
    fwng_span0	=	geom	[	4	]
    fwng_dihe0	=	geom	[	5	]
    fwng_swep0	=	geom	[	6	]
    fwng_chrd1	=	geom	[	7	]
    fwng_inci1	=	geom	[	8	]
    fwng_span1	=	geom	[	9	]
    fwng_dihe1	=	geom	[	10	]
    fwng_swep1	=	geom	[	11	]
    fwng_chrd2	=	geom	[	12	]
    fwng_inci2	=	geom	[	13	]
    fwng_span2	=	geom	[	14	]
    fwng_dihe2	=	geom	[	15	]
    fwng_swep2	=	geom	[	16	]
    fwng_chrd3	=	geom	[	17	]
    fwng_inci3	=	geom	[	18	]
    rwng_xpos0	=	geom	[	19	]
    rwng_zpos0	=	geom	[	20	]
    rwng_chrd0	=	geom	[	21	]
    rwng_inci0	=	geom	[	22	]
    rwng_span0	=	geom	[	23	]
    rwng_dihe0	=	geom	[	24	]
    rwng_swep0	=	geom	[	25	]
    rwng_chrd1	=	geom	[	26	]
    rwng_inci1	=	geom	[	27	]
    rwng_span1	=	geom	[	28	]
    rwng_dihe1	=	geom	[	29	]
    rwng_swep1	=	geom	[	30	]
    rwng_chrd2	=	geom	[	31	]
    rwng_inci2	=	geom	[	32	]
    rwng_span2	=	geom	[	33	]
    rwng_dihe2	=	geom	[	34	]
    rwng_swep2	=	geom	[	35	]
    rwng_chrd3	=	geom	[	36	]
    rwng_inci3	=	geom	[	37	]
    tail_xpos	=	geom	[	38	]
    tail_zpos	=	geom	[	39	]
    tail_span	=	geom	[	40	]
    tail_chrd0	=	geom	[	41	]
    tail_chrd1	=	geom	[	42	]
    tail_swep	=	geom	[	43	]

    fwng_span = 2*(fwng_span0+fwng_span1+fwng_span2)
    rwng_span = 2*(rwng_span0+rwng_span1+rwng_span2)

    fuselage_width = 0.06*max(fwng_span,rwng_span) #6% of wingspan
    fuselage_length = 3 + max(rwng_xpos0,tail_xpos) - fwng_xpos0 #3 metres more than distance between fwng and rwng/tail
    fuselage_length1= max(fwng_span,rwng_span)#another fuselage length based on fineness ratio of 1
    # print("span = {}".format(fuselage_length1))
    imag_fuse_length= max(fuselage_length,fuselage_length1)
    fuselage_vol = imag_fuse_length*np.pi*(0.5*fuselage_width)**2 #imaginary fuselage volume
    # print(fuselage_vol)
    fuselage_width = 2*np.sqrt(fuselage_vol/(fuselage_length*np.pi)) # new fuselage width based on fixed volume
    # print(fuselage_width)
    # print("fuselage_length = {}".format(fuselage_length))

    #convert all angles to radians
    fwng_inci0	=fwng_inci0*np.pi/180
    fwng_dihe0	=fwng_dihe0*np.pi/180
    fwng_swep0	=fwng_swep0*np.pi/180
    fwng_inci1	=fwng_inci1*np.pi/180
    fwng_dihe1	=fwng_dihe1*np.pi/180
    fwng_swep1	=fwng_swep1*np.pi/180
    fwng_inci2	=fwng_inci2*np.pi/180
    fwng_dihe2	=fwng_dihe2*np.pi/180
    fwng_swep2	=fwng_swep2*np.pi/180
    fwng_inci3	=fwng_inci3*np.pi/180
    rwng_inci0	=rwng_inci0*np.pi/180
    rwng_dihe0	=rwng_dihe0*np.pi/180
    rwng_swep0	=rwng_swep0*np.pi/180
    rwng_inci1	=rwng_inci1*np.pi/180
    rwng_dihe1	=rwng_dihe1*np.pi/180
    rwng_swep1	=rwng_swep1*np.pi/180
    rwng_inci2	=rwng_inci2*np.pi/180
    rwng_dihe2	=rwng_dihe2*np.pi/180
    rwng_swep2	=rwng_swep2*np.pi/180
    rwng_inci3	=rwng_inci3*np.pi/180
    tail_swep   =tail_swep*np.pi/180

    strips = 32

    cbar = 0.25 #geometry refers to this chordwise position, but avl needs leading edge

    fwng_0_le_pnt_x = fwng_xpos0-cbar*fwng_chrd0*np.cos(fwng_inci0)
    fwng_0_le_pnt_y = 0
    fwng_0_le_pnt_z = fwng_zpos0+cbar*fwng_chrd0*np.sin(fwng_inci0)
    fwng_1_le_pnt_x = fwng_0_le_pnt_x+cbar*fwng_chrd0*np.cos(fwng_inci0)+fwng_span0*np.tan(fwng_swep0)-cbar*fwng_chrd1*np.cos(fwng_inci1)
    fwng_1_le_pnt_y = fwng_0_le_pnt_y + fwng_span0
    fwng_1_le_pnt_z = fwng_0_le_pnt_z-cbar*fwng_chrd0*np.sin(fwng_inci0)+fwng_span0*np.tan(fwng_dihe0)+cbar*fwng_chrd1*np.sin(fwng_inci1)
    fwng_2_le_pnt_x = fwng_1_le_pnt_x+cbar*fwng_chrd1*np.cos(fwng_inci1)+fwng_span1*np.tan(fwng_swep1)-cbar*fwng_chrd2*np.cos(fwng_inci2)
    fwng_2_le_pnt_y = fwng_1_le_pnt_y + fwng_span1
    fwng_2_le_pnt_z = fwng_1_le_pnt_z-cbar*fwng_chrd1*np.sin(fwng_inci1)+fwng_span1*np.tan(fwng_dihe1)+cbar*fwng_chrd2*np.sin(fwng_inci2)
    fwng_3_le_pnt_x = fwng_2_le_pnt_x+cbar*fwng_chrd2*np.cos(fwng_inci2)+fwng_span2*np.tan(fwng_swep2)-cbar*fwng_chrd3*np.cos(fwng_inci3)
    fwng_3_le_pnt_y = fwng_2_le_pnt_y + fwng_span2
    fwng_3_le_pnt_z = fwng_2_le_pnt_z-cbar*fwng_chrd2*np.sin(fwng_inci2)+fwng_span2*np.tan(fwng_dihe2)+cbar*fwng_chrd3*np.sin(fwng_inci3)

    fwng_0_le_pnt = avl.Point(fwng_0_le_pnt_x,fwng_0_le_pnt_y,fwng_0_le_pnt_z)
    fwng_1_le_pnt = avl.Point(fwng_1_le_pnt_x,fwng_1_le_pnt_y,fwng_1_le_pnt_z)
    fwng_2_le_pnt = avl.Point(fwng_2_le_pnt_x,fwng_2_le_pnt_y,fwng_2_le_pnt_z)
    fwng_3_le_pnt = avl.Point(fwng_3_le_pnt_x,fwng_3_le_pnt_y,fwng_3_le_pnt_z)

    fwng_inci0 = np.rad2deg(fwng_inci0)
    fwng_inci1 = np.rad2deg(fwng_inci1)
    fwng_inci2 = np.rad2deg(fwng_inci2)
    fwng_inci3 = np.rad2deg(fwng_inci3)

    fwng_0 = avl.Section(leading_edge_point=fwng_0_le_pnt,
                            chord=fwng_chrd0,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=fwng_inci0)
    fwng_1 = avl.Section(leading_edge_point=fwng_1_le_pnt,
                            chord=fwng_chrd1,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=fwng_inci1)
    fwng_2 = avl.Section(leading_edge_point=fwng_2_le_pnt,
                            chord=fwng_chrd2,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=fwng_inci2)
    fwng_3 = avl.Section(leading_edge_point=fwng_3_le_pnt,
                            chord=fwng_chrd3,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=fwng_inci3)


    # y_duplicate=0.0 duplicates the wing over a XZ-plane at Y=0.0
    forward_wing = avl.Surface(name='forward_wing',
                    n_chordwise=8,
                    chord_spacing=avl.Spacing.cosine,
                    n_spanwise=strips,
                    span_spacing=avl.Spacing.cosine,
                    y_duplicate=0.0,
                    sections=[fwng_0, fwng_1, fwng_2, fwng_3])

    ##define rear wing

    rwng_0_le_pnt_x = rwng_xpos0-cbar*rwng_chrd0*np.cos(rwng_inci0)
    rwng_0_le_pnt_y = 0
    rwng_0_le_pnt_z = rwng_zpos0+cbar*rwng_chrd0*np.sin(rwng_inci0)
    rwng_1_le_pnt_x = rwng_0_le_pnt_x+cbar*rwng_chrd0*np.cos(rwng_inci0)+rwng_span0*np.tan(rwng_swep0)-cbar*rwng_chrd1*np.cos(rwng_inci1)
    rwng_1_le_pnt_y = rwng_0_le_pnt_y + rwng_span0
    rwng_1_le_pnt_z = rwng_0_le_pnt_z-cbar*rwng_chrd0*np.sin(rwng_inci0)+rwng_span0*np.tan(rwng_dihe0)+cbar*rwng_chrd1*np.sin(rwng_inci1)
    rwng_2_le_pnt_x = rwng_1_le_pnt_x+cbar*rwng_chrd1*np.cos(rwng_inci1)+rwng_span1*np.tan(rwng_swep1)-cbar*rwng_chrd2*np.cos(rwng_inci2)
    rwng_2_le_pnt_y = rwng_1_le_pnt_y + rwng_span1
    rwng_2_le_pnt_z = rwng_1_le_pnt_z-cbar*rwng_chrd1*np.sin(rwng_inci1)+rwng_span1*np.tan(rwng_dihe1)+cbar*rwng_chrd2*np.sin(rwng_inci2)
    rwng_3_le_pnt_x = rwng_2_le_pnt_x+cbar*rwng_chrd2*np.cos(rwng_inci2)+rwng_span2*np.tan(rwng_swep2)-cbar*rwng_chrd3*np.cos(rwng_inci3)
    rwng_3_le_pnt_y = rwng_2_le_pnt_y + rwng_span2
    rwng_3_le_pnt_z = rwng_2_le_pnt_z-cbar*rwng_chrd2*np.sin(rwng_inci2)+rwng_span2*np.tan(rwng_dihe2)+cbar*rwng_chrd3*np.sin(rwng_inci3)

    rwng_0_le_pnt = avl.Point(rwng_0_le_pnt_x,rwng_0_le_pnt_y,rwng_0_le_pnt_z)
    rwng_1_le_pnt = avl.Point(rwng_1_le_pnt_x,rwng_1_le_pnt_y,rwng_1_le_pnt_z)
    rwng_2_le_pnt = avl.Point(rwng_2_le_pnt_x,rwng_2_le_pnt_y,rwng_2_le_pnt_z)
    rwng_3_le_pnt = avl.Point(rwng_3_le_pnt_x,rwng_3_le_pnt_y,rwng_3_le_pnt_z)

    rwng_inci0 = np.rad2deg(rwng_inci0)
    rwng_inci1 = np.rad2deg(rwng_inci1)
    rwng_inci2 = np.rad2deg(rwng_inci2)
    rwng_inci3 = np.rad2deg(rwng_inci3)

    rwng_0 = avl.Section(leading_edge_point=rwng_0_le_pnt,
                            chord=rwng_chrd0,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=rwng_inci0)
    rwng_1 = avl.Section(leading_edge_point=rwng_1_le_pnt,
                            chord=rwng_chrd1,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=rwng_inci1)
    rwng_2 = avl.Section(leading_edge_point=rwng_2_le_pnt,
                            chord=rwng_chrd2,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=rwng_inci2)
    rwng_3 = avl.Section(leading_edge_point=rwng_3_le_pnt,
                            chord=rwng_chrd3,
                            airfoil=avl.NacaAirfoil('2414'),
                            angle=rwng_inci3)

    # y_duplicate=0.0 duplicates the wing over a XZ-plane at Y=0.0
    rear_wing = avl.Surface(name='rear_wing',
                    n_chordwise=8,
                    chord_spacing=avl.Spacing.cosine,
                    n_spanwise=strips,
                    span_spacing=avl.Spacing.cosine,
                    y_duplicate=0.0,
                    sections=[rwng_0, rwng_1, rwng_2, rwng_3])

    ##define tail wing

    wing_span = tail_span
    wing_root_chord = tail_chrd0
    wing_tip_chord = tail_chrd1

    wing_root_le_pnt = avl.Point(tail_xpos-cbar*tail_chrd0, 0, tail_zpos)
    wing_tip_le_pnt = avl.Point(tail_xpos+tail_span*np.tan(tail_swep)-cbar*tail_chrd0,
                                0,
                                tail_zpos+tail_span)

    root_section = avl.Section(leading_edge_point=wing_root_le_pnt,
                            chord=wing_root_chord,
                            airfoil=avl.NacaAirfoil('0012'))
    tip_section = avl.Section(leading_edge_point=wing_tip_le_pnt,
                            chord=wing_tip_chord,
                            airfoil=avl.NacaAirfoil('0010'))

    # y_duplicate=0.0 duplicates the wing over a XZ-plane at Y=0.0
    tail_wing = avl.Surface(name='tail_wing',
                    n_chordwise=8,
                    chord_spacing=avl.Spacing.cosine,
                    n_spanwise=10,
                    span_spacing=avl.Spacing.cosine,
                    sections=[root_section, tip_section])

    #define fuselage surfaces

    # fuselage_centreline_le_pnt = avl.Point(fwng_xpos0-3, 0, 0)
    fuselage_side_le_pnt0= avl.Point(fwng_xpos0-3,-fuselage_width/2,0)
    fuselage_side_le_pnt = avl.Point(fwng_xpos0-3,fuselage_width/2,0)

    # centreline = avl.Section(leading_edge_point=fuselage_centreline_le_pnt,
    #                         chord=fuselage_length,
    #                         airfoil=avl.NacaAirfoil('0000'))
    fuselage_side0 = avl.Section(leading_edge_point=fuselage_side_le_pnt0,
                            chord=fuselage_length,
                            airfoil=avl.NacaAirfoil('0000'))
    fuselage_side = avl.Section(leading_edge_point=fuselage_side_le_pnt,
                            chord=fuselage_length,
                            airfoil=avl.NacaAirfoil('0000'))

    # y_duplicate=0.0 duplicates the wing over a XZ-plane at Y=0.0
    fuselage_horizontal = avl.Surface(name='fuselage_horizontal',
                    n_chordwise=25,
                    chord_spacing=avl.Spacing.cosine,
                    n_spanwise=20,
                    span_spacing=avl.Spacing.cosine,
                    # y_duplicate=0.0,
                    sections=[fuselage_side0,fuselage_side])


    fuselage_bottom_le_pnt = avl.Point(fwng_xpos0-3, 0, -fuselage_width/2)
    fuselage_top_le_pnt = avl.Point(fwng_xpos0-3,0,fuselage_width/2)

    fuse_bottom = avl.Section(leading_edge_point=fuselage_bottom_le_pnt,
                            chord=fuselage_length,
                            airfoil=avl.NacaAirfoil('0000'))
    fuse_top = avl.Section(leading_edge_point=fuselage_top_le_pnt,
                            chord=fuselage_length,
                            airfoil=avl.NacaAirfoil('0000'))

    # y_duplicate=0.0 duplicates the wing over a XZ-plane at Y=0.0
    fuselage_vertical = avl.Surface(name='fuselage_vertical',
                    n_chordwise=25,
                    chord_spacing=avl.Spacing.cosine,
                    n_spanwise=5,
                    span_spacing=avl.Spacing.cosine,
                    sections=[fuse_bottom,fuse_top])

    ## define operating point and reference values

    mach = 0.3

    fwng_area = 2*(
        (fwng_span0*0.5*(fwng_chrd0+fwng_chrd1))*np.cos(fwng_dihe0)+
        (fwng_span1*0.5*(fwng_chrd1+fwng_chrd2))*np.cos(fwng_dihe1)+
        (fwng_span2*0.5*(fwng_chrd2+fwng_chrd3))*np.cos(fwng_dihe2)
    )
    rwng_area = 2*(
        (rwng_span0*0.5*(rwng_chrd0+rwng_chrd1))*np.cos(rwng_dihe0)+
        (rwng_span1*0.5*(rwng_chrd1+rwng_chrd2))*np.cos(rwng_dihe1)+
        (rwng_span2*0.5*(rwng_chrd2+rwng_chrd3))*np.cos(rwng_dihe2)
    )

    fwng_area_not_projected = 2*(
        (fwng_span0*0.5*(fwng_chrd0+fwng_chrd1))+
        (fwng_span1*0.5*(fwng_chrd1+fwng_chrd2))+
        (fwng_span2*0.5*(fwng_chrd2+fwng_chrd3))
    )
    rwng_area_not_projected = 2*(
        (rwng_span0*0.5*(rwng_chrd0+rwng_chrd1))+
        (rwng_span1*0.5*(rwng_chrd1+rwng_chrd2))+
        (rwng_span2*0.5*(rwng_chrd2+rwng_chrd3))
    )
    wing_area = fwng_area + rwng_area

    # print("S = {}".format(wing_area))

    wing_mac = wing_area/(max(fwng_span,rwng_span)) #need better MAC calulation
    # print(wing_mac)

    # calculate the m.a.c. leading edge location
    def mac_le_pnt(root_chord, tip_chord, root_pnt, tip_pnt):
        pnt = ((2*root_chord*root_pnt[dim] +
                root_chord*tip_pnt[dim] + 
                tip_chord*root_pnt[dim] +
                2*tip_chord*tip_pnt[dim]) / 
            (3*(root_chord+tip_chord))
            for dim in range(3))
        return avl.Point(*pnt)

    le_pnt = mac_le_pnt(wing_root_chord, wing_tip_chord,
                        wing_root_le_pnt, wing_tip_le_pnt)

    # ref_pnt = avl.Point(x=le_pnt.x + 0.25*wing_mac,
    #                     y=le_pnt.y, z=le_pnt.z)

    ref_pnt = avl.Point(x = 0.5*(fwng_xpos0+rwng_xpos0),
                        y=0, z=0)

    ## put it into avl format

    aircraft = avl.Aircraft(name='aircraft',
                            reference_area=wing_area,
                            reference_chord=wing_mac,
                            reference_span=wing_span,
                            reference_point=ref_pnt,
                            mach=mach,
                            surfaces=[forward_wing, rear_wing, tail_wing, fuselage_horizontal, fuselage_vertical])

    ##start session with specified geometry
    # default_config.local_copy() 
    # session = avl.Session(geometry=aircraft)
    



    # create a function for showing the Trefftz plot, since we'll be using it more often
    def show_treffz(session):
        if 'gs_bin' in session.config.settings:
            images = session.save_trefftz_plots()
            for img in images:
                avl.show_image(img)
        else:
            for idx, _ in enumerate(session.cases):
                session.show_trefftz_plot(idx+1) # cases start from 1


    base_case = avl.Case(name='sweep', mach=0.3)

    alphas = np.linspace(-2,8,11)#sweep through alphas change this to -2 8 11
    all_cases = avl.create_sweep_cases(base_case=base_case,
                                    parameters=[{'name': 'alpha',
                                                    'values': alphas}
                                                    ])
    session = avl.Session(geometry=aircraft, cases=all_cases)
    # session.show_geometry()
    # show_treffz(session)

    partitions = avl.partitioned_cases(all_cases)
    results = {}
    for partition in partitions:
        session = avl.Session(geometry=aircraft, cases=partition)
        results.update(session.run_all_cases())
        
    # Write everything to json
    with open('all_cases.json', 'w') as f:
        f.write(json.dumps(results))
        f.close()

    



    ## drag buildup
    # 0 = fuselage
    # 1 = forward wing
    # 2 = rear wing
    # 3 = tail fin


    laminar_ratio = [0.0,0,0,0]#proportion of laminar flow over components
    length = [fuselage_length,fwng_chrd0,rwng_chrd0,tail_chrd0]#characteristic lengths of components (fuselage length, wing section chords)
    dist = np.sqrt(np.square(fwng_xpos0-rwng_xpos0)+np.square(fwng_zpos0-rwng_zpos0))#distance between origins of wings
    dist_c_forward = dist/rwng_chrd0#effective distance of rear wing from forward wing, normalised to chord of rear wing
    dist_c_rear = dist/fwng_chrd0


    interference_factor = [1,1,1,1.03]#Q, from Raymer
    interference_factor[1] = max(1,1.5-0.25*dist_c_forward)#interference factor adjustments for wings close together
    interference_factor[2] = max(1,1.5-0.25*dist_c_rear)#source: adapted from raymer
    # print(interference_factor)
    S_wet = [np.pi*fuselage_width*fuselage_length+2*(np.pi*(fuselage_width/2)**2),2.2*fwng_area_not_projected,2.2*rwng_area_not_projected,2.2*(tail_chrd0+tail_chrd1)*0.5*tail_span]#wetted surface area
    # print(sum(S_wet))
    fineness = fuselage_length/fuselage_width#fineness ratio
    # print(fineness)
    FF_fuse = 0.9 + (5/(fineness**1.5)) + (fineness/400)#form factor, from Raymer
    thickness = [0.14,0.14,0.12]#thickness of aerofoils of wing sections
    FF_forwardwing = 1 + (0.6/0.25)*thickness[0] + 100*(thickness[0]**4)
    FF_rearwing = 1 + (0.6/0.25)*thickness[1] + 100*(thickness[1]**4)
    FF_tail = 1 + (0.6/0.25)*thickness[2] + 100*(thickness[2]**4)
    FF = [FF_fuse,FF_forwardwing,FF_rearwing,FF_tail]

    CD0_components = 0
    for n in [0,1,2,3]:
        # print(n)
        # print(length[n])
        Re = rho*40*length[n]/(1.789e-5)#component Reynolds
        Re_cutoff = 38.21*((length[n]/0.634e-5)**1.053)#cutoff reynolds, from raymer, not really useful tbh
        Re = min(Re,Re_cutoff)
        Cfc = (laminar_ratio[n]*1.328 + (1-laminar_ratio[n])*4.5)/sqrt(Re)#skin friction coeff, apportioning between laminar and turbulent
        # print("Cfc = {}".format(Cfc))
        CD0_comp = (Cfc*FF[n]*interference_factor[n]*S_wet[n])/wing_area#CD0 contribution from component
        # print(CD0_comp)
        if n == 0: #if fuselage
            # print(np.pi*(fuselage_width/2)**2)
            CD0_comp = CD0_comp + (0.3 * (np.pi*(fuselage_width/2)**2))/wing_area #form drag of fuselage
        # else:
        #     a = 1
        # print("CD0_comp = {}".format(CD0_comp))
        CD0_components = CD0_components + CD0_comp#add them up

    # print(CD0_components)



    alpha = []
    CL = []
    CL_fwng = []
    CL_rwng = []
    CD = []
    LD = []
    CD_ind21 = []
    # print(results[1]["SurfaceForces"]["forward_wing"]["CL"].keys())
    for n in results.keys():#for each run case (from alpha sweep)
        alpha.append(results[n]['Totals']['Alpha'])#find values for lift curve
        CL.append(results[n]["Totals"]["CLtot"])
        CL_fwng.append(results[n]["SurfaceForces"]["forward_wing"]["CL"])
        CL_rwng.append(results[n]["SurfaceForces"]["rear_wing"]["CL"])
        # alpha_induced_drag = 5 #the alpha at which the mutual induced drag is calculated
        CL_fwng_ind = results[n]["SurfaceForces"]["forward_wing"]["CL"]
        # print(CL_fwng_ind)
        CL_rwng_ind = results[n]["SurfaceForces"]["rear_wing"]["CL"]
        # print(CL_rwng_ind)
        delta = np.atan((rwng_zpos0-fwng_zpos0)/max(0.001,(rwng_xpos0-fwng_xpos0)))#relative wing angle with reference to fuselage
        h = dist*np.sin(delta-np.deg2rad(results[n]['Totals']['Alpha']))
        # h = fwng_zpos0-rwng_zpos0#height diff         
        c = (fwng_span-rwng_span)/2#wing exceedance
        b = (fwng_span+rwng_span)/2#avg wingspan
        sigma = 0.125*np.log((b**2+h**2)/max(0.00001,(c**2+h**2)))
        # print(sigma)
        CD_ind21.append(2*sigma*((CL_fwng_ind*CL_rwng_ind*fwng_area_not_projected*rwng_area_not_projected)/(np.pi*(fwng_area_not_projected+rwng_area_not_projected)*fwng_span*rwng_span)))
    # print(CD_ind21)
    # CD0_components = CD0_components + CD_ind21
    
    for n in results.keys():
        CD.append(results[n]["Totals"]["CDtot"]+CD0_components+CD_ind21[n-1])#add drag buildup
        LD.append(CL[n-1]/CD[n-1])

    # fig, ax = plt.subplots()
    # ax.plot(alpha,CD,".")
    # plt.show()
    # fig, ax = plt.subplots()
    # ax.plot(alpha,CL,".")
    # plt.show()
    # fig, ax = plt.subplots()
    # ax.plot(alpha,LD,".")
    # plt.show()

    # print(alpha)
    # print(results[alpha.index(max(alpha))]["StripForces"]["forward_wing"]["Yle"])
    # print(max(LD))

    datum = round(alpha[LD.index(max(LD))]+3) #get CL dists at best LD
    # print(results.keys())
    #obtain lift distributions and sort them
    forward_wing_y_array = results[datum]['StripForces']['forward_wing']['Yle']
    forward_wing_cl_array = results[datum]['StripForces']['forward_wing']['c_cl']
    rear_wing_y_array = results[datum]['StripForces']['rear_wing']['Yle']
    rear_wing_cl_array = results[datum]['StripForces']['rear_wing']['c_cl']

    # Convert to numpy arrays
    arr_p_f = np.array(forward_wing_y_array)#primary array, forward
    arr_a_f = np.array(forward_wing_cl_array)#associate array, forward
    arr_p_r = np.array(rear_wing_y_array)#primary array, rear
    arr_a_r = np.array(rear_wing_cl_array)#associate array, rear

    # Get the sorting index map
    indices_f = np.argsort(arr_p_f)
    indices_r = np.argsort(arr_p_r)

    # Apply index map to both arrays
    forward_wing_y_array = arr_p_f[indices_f]
    forward_wing_cl_array = arr_a_f[indices_f]
    rear_wing_y_array = arr_p_r[indices_r]
    rear_wing_cl_array = arr_a_r[indices_r]

    # print(forward_wing_y_array)
    rwng_ypts = [-(rwng_span0+rwng_span1+rwng_span2),-(rwng_span0+rwng_span1),-(rwng_span0),0,(rwng_span0),(rwng_span0+rwng_span1),(rwng_span0+rwng_span1+rwng_span2)]
    rwng_chrd_pts = [rwng_chrd3,rwng_chrd2,rwng_chrd1,rwng_chrd0,rwng_chrd1,rwng_chrd2,rwng_chrd3]
    rwng_dihe_pts = [rwng_dihe2,rwng_dihe1,rwng_dihe0,rwng_dihe0,rwng_dihe1,rwng_dihe2]
    rwng_swep_pts = [rwng_swep2,rwng_swep1,rwng_swep0,rwng_swep0,rwng_swep1,rwng_swep2]

    fwng_ypts = [-(fwng_span0+fwng_span1+fwng_span2),-(fwng_span0+fwng_span1),-(fwng_span0),0,(fwng_span0),(fwng_span0+fwng_span1),(fwng_span0+fwng_span1+fwng_span2)]
    fwng_chrd_pts = [fwng_chrd3,fwng_chrd2,fwng_chrd1,fwng_chrd0,fwng_chrd1,fwng_chrd2,fwng_chrd3]
    fwng_dihe_pts = [fwng_dihe2,fwng_dihe1,fwng_dihe0,fwng_dihe0,fwng_dihe1,fwng_dihe2]
    fwng_swep_pts = [fwng_swep2,fwng_swep1,fwng_swep0,fwng_swep0,fwng_swep1,fwng_swep2]
    # print(rear_wing_y_array)
    # print(rwng_ypts)

    forward_wing_cellwidth = []
    forward_wing_load = []
    forward_wing_chord = []
    forward_wing_dihe = []
    forward_wing_swep = []
    forward_wing_cellsparlength = []
    forward_wing_shear = [0]
    forward_wing_bend = [0]

    rear_wing_cellwidth = []
    rear_wing_load = []
    rear_wing_chord = []
    rear_wing_dihe = []
    rear_wing_swep = []
    rear_wing_cellsparlength = []
    rear_wing_shear = [0]
    rear_wing_bend = [0]
    x = len(rear_wing_y_array)
    y = 0
    while y < x-1: #linear interpolation to get chord, sweep, dihedral for every y pt
        forward_wing_cellwidth.append(forward_wing_y_array[y+1]-forward_wing_y_array[y])
        rear_wing_cellwidth.append(rear_wing_y_array[y+1]-rear_wing_y_array[y])

        if rear_wing_y_array[y]>rwng_ypts[0] and rear_wing_y_array[y]<rwng_ypts[1]:
            rear_wing_chord.append(rwng_chrd_pts[0]+(rwng_chrd_pts[1]-rwng_chrd_pts[0])*(rear_wing_y_array[y]-rwng_ypts[0])/((rwng_ypts[1]-rwng_ypts[0])))
            rear_wing_dihe.append(rwng_dihe_pts[0])
            rear_wing_swep.append(rwng_swep_pts[0])
            rear_wing_cellsparlength.append(np.sqrt((rear_wing_cellwidth[y]/np.cos(rear_wing_swep[y]))**2+(rear_wing_cellwidth[y]*np.tan(rear_wing_dihe[y]))**2))

        if rear_wing_y_array[y]>rwng_ypts[1] and rear_wing_y_array[y]<rwng_ypts[2]:
            rear_wing_chord.append(rwng_chrd_pts[1]+(rwng_chrd_pts[2]-rwng_chrd_pts[1])*(rear_wing_y_array[y]-rwng_ypts[1])/((rwng_ypts[2]-rwng_ypts[1])))
            rear_wing_dihe.append(rwng_dihe_pts[1])
            rear_wing_swep.append(rwng_swep_pts[1])
            rear_wing_cellsparlength.append(np.sqrt((rear_wing_cellwidth[y]/np.cos(rear_wing_swep[y]))**2+(rear_wing_cellwidth[y]*np.tan(rear_wing_dihe[y]))**2))

        if rear_wing_y_array[y]>rwng_ypts[2] and rear_wing_y_array[y]<rwng_ypts[3]:
            rear_wing_chord.append(rwng_chrd_pts[2]+(rwng_chrd_pts[3]-rwng_chrd_pts[2])*(rear_wing_y_array[y]-rwng_ypts[2])/((rwng_ypts[3]-rwng_ypts[2])))
            rear_wing_dihe.append(rwng_dihe_pts[2])
            rear_wing_swep.append(rwng_swep_pts[2])
            rear_wing_cellsparlength.append(np.sqrt((rear_wing_cellwidth[y]/np.cos(rear_wing_swep[y]))**2+(rear_wing_cellwidth[y]*np.tan(rear_wing_dihe[y]))**2))

        if rear_wing_y_array[y]>rwng_ypts[3] and rear_wing_y_array[y]<rwng_ypts[4]:
            rear_wing_chord.append(rwng_chrd_pts[3]+(rwng_chrd_pts[4]-rwng_chrd_pts[3])*(rear_wing_y_array[y]-rwng_ypts[3])/((rwng_ypts[4]-rwng_ypts[3])))
            rear_wing_dihe.append(rwng_dihe_pts[3])
            rear_wing_swep.append(rwng_swep_pts[3])
            rear_wing_cellsparlength.append(np.sqrt((rear_wing_cellwidth[y]/np.cos(rear_wing_swep[y]))**2+(rear_wing_cellwidth[y]*np.tan(rear_wing_dihe[y]))**2))

        if rear_wing_y_array[y]>rwng_ypts[4] and rear_wing_y_array[y]<rwng_ypts[5]:
            rear_wing_chord.append(rwng_chrd_pts[4]+(rwng_chrd_pts[5]-rwng_chrd_pts[4])*(rear_wing_y_array[y]-rwng_ypts[4])/((rwng_ypts[5]-rwng_ypts[4])))
            rear_wing_dihe.append(rwng_dihe_pts[4])
            rear_wing_swep.append(rwng_swep_pts[4])
            rear_wing_cellsparlength.append(np.sqrt((rear_wing_cellwidth[y]/np.cos(rear_wing_swep[y]))**2+(rear_wing_cellwidth[y]*np.tan(rear_wing_dihe[y]))**2))

        if rear_wing_y_array[y]>rwng_ypts[5] and rear_wing_y_array[y]<rwng_ypts[6]:
            rear_wing_chord.append(rwng_chrd_pts[5]+(rwng_chrd_pts[6]-rwng_chrd_pts[5])*(rear_wing_y_array[y]-rwng_ypts[5])/((rwng_ypts[6]-rwng_ypts[5])))
            rear_wing_dihe.append(rwng_dihe_pts[5])
            rear_wing_swep.append(rwng_swep_pts[5])
            rear_wing_cellsparlength.append(np.sqrt((rear_wing_cellwidth[y]/np.cos(rear_wing_swep[y]))**2+(rear_wing_cellwidth[y]*np.tan(rear_wing_dihe[y]))**2))

        if forward_wing_y_array[y]>fwng_ypts[0] and forward_wing_y_array[y]<fwng_ypts[1]:
            forward_wing_chord.append(fwng_chrd_pts[0]+(fwng_chrd_pts[1]-fwng_chrd_pts[0])*(forward_wing_y_array[y]-fwng_ypts[0])/((fwng_ypts[1]-fwng_ypts[0])))
            forward_wing_dihe.append(fwng_dihe_pts[0])
            forward_wing_swep.append(fwng_swep_pts[0])
            forward_wing_cellsparlength.append(np.sqrt((forward_wing_cellwidth[y]/np.cos(forward_wing_swep[y]))**2+(forward_wing_cellwidth[y]*np.tan(forward_wing_dihe[y]))**2))

        if forward_wing_y_array[y]>fwng_ypts[1] and forward_wing_y_array[y]<fwng_ypts[2]:
            forward_wing_chord.append(fwng_chrd_pts[1]+(fwng_chrd_pts[2]-fwng_chrd_pts[1])*(forward_wing_y_array[y]-fwng_ypts[1])/((fwng_ypts[2]-fwng_ypts[1])))
            forward_wing_dihe.append(fwng_dihe_pts[1])
            forward_wing_swep.append(fwng_swep_pts[1])
            forward_wing_cellsparlength.append(np.sqrt((forward_wing_cellwidth[y]/np.cos(forward_wing_swep[y]))**2+(forward_wing_cellwidth[y]*np.tan(forward_wing_dihe[y]))**2))

        if forward_wing_y_array[y]>fwng_ypts[2] and forward_wing_y_array[y]<fwng_ypts[3]:
            forward_wing_chord.append(fwng_chrd_pts[2]+(fwng_chrd_pts[3]-fwng_chrd_pts[2])*(forward_wing_y_array[y]-fwng_ypts[2])/((fwng_ypts[3]-fwng_ypts[2])))
            forward_wing_dihe.append(fwng_dihe_pts[2])
            forward_wing_swep.append(fwng_swep_pts[2])
            forward_wing_cellsparlength.append(np.sqrt((forward_wing_cellwidth[y]/np.cos(forward_wing_swep[y]))**2+(forward_wing_cellwidth[y]*np.tan(forward_wing_dihe[y]))**2))

        if forward_wing_y_array[y]>fwng_ypts[3] and forward_wing_y_array[y]<fwng_ypts[4]:
            forward_wing_chord.append(fwng_chrd_pts[3]+(fwng_chrd_pts[4]-fwng_chrd_pts[3])*(forward_wing_y_array[y]-fwng_ypts[3])/((fwng_ypts[4]-fwng_ypts[3])))
            forward_wing_dihe.append(fwng_dihe_pts[3])
            forward_wing_swep.append(fwng_swep_pts[3])
            forward_wing_cellsparlength.append(np.sqrt((forward_wing_cellwidth[y]/np.cos(forward_wing_swep[y]))**2+(forward_wing_cellwidth[y]*np.tan(forward_wing_dihe[y]))**2))

        if forward_wing_y_array[y]>fwng_ypts[4] and forward_wing_y_array[y]<fwng_ypts[5]:
            forward_wing_chord.append(fwng_chrd_pts[4]+(fwng_chrd_pts[5]-fwng_chrd_pts[4])*(forward_wing_y_array[y]-fwng_ypts[4])/((fwng_ypts[5]-fwng_ypts[4])))
            forward_wing_dihe.append(fwng_dihe_pts[4])
            forward_wing_swep.append(fwng_swep_pts[4])
            forward_wing_cellsparlength.append(np.sqrt((forward_wing_cellwidth[y]/np.cos(forward_wing_swep[y]))**2+(forward_wing_cellwidth[y]*np.tan(forward_wing_dihe[y]))**2))

        if forward_wing_y_array[y]>fwng_ypts[5] and forward_wing_y_array[y]<fwng_ypts[6]:
            forward_wing_chord.append(fwng_chrd_pts[5]+(fwng_chrd_pts[6]-fwng_chrd_pts[5])*(forward_wing_y_array[y]-fwng_ypts[5])/((fwng_ypts[6]-fwng_ypts[5])))
            forward_wing_dihe.append(fwng_dihe_pts[5])
            forward_wing_swep.append(fwng_swep_pts[5])
            forward_wing_cellsparlength.append(np.sqrt((forward_wing_cellwidth[y]/np.cos(forward_wing_swep[y]))**2+(forward_wing_cellwidth[y]*np.tan(forward_wing_dihe[y]))**2))

        #make spar length array 
        #make chord array
        y = y + 1

    # print(rear_wing_chord)
    rear_wing_cellwidth.append(0)
    rear_wing_chord.append(rwng_chrd3)
    rear_wing_dihe.append(rwng_dihe2)
    rear_wing_swep.append(rwng_swep2)
    rear_wing_cellsparlength.append(0)

    forward_wing_cellwidth.append(0)
    forward_wing_chord.append(fwng_chrd3)
    forward_wing_dihe.append(fwng_dihe2)
    forward_wing_swep.append(fwng_swep2)
    forward_wing_cellsparlength.append(0)
    CD1 = min(CD)#pick two points from drag bucket to calculate k and CD0
    CL1 = CL[CD.index(min(CD))]
    CD2 = CD[LD.index(max(LD))]
    CL2 = CL[LD.index(max(LD))]

    k = abs((CD2-CD1)/((CL2**2)-(CL1**2)))
    CD0 = abs(CD1 - k * (CL1**2))

    # print(CD0)
    # print(k)
    CL_LDmax = sqrt(CD0/k)#use drag polar to find POI
    CD_LDmax = CD0 + k * CL_LDmax**2

    CL_minpower = sqrt(3*CD0/k)
    CD_minpower = CD0 + k * CL_minpower**2

    CL_max = 1.5
    CD_CLmax = CD0 + k * CL_max**2

    forward_wing_mass = 0
    rear_wing_mass = 0
    total_mass = mass_payload * 4 #initial guess for total mass
    total_mass_minus_one = total_mass
    # print(total_mass)
    massplot = [total_mass]

    prop_eff = 0.7
    SFC = 7.7e-7 #N/Ws #8.7e-7 for jab, 7.7e-7 fo the rotax one, 8.5 else
    x = 0
    tol = 0.002
    err = 1
    while (err > tol and x < 1000 and err < 10 and total_mass < 100000000) or x < 5:
        #rearrange Breguet range and endurance eqns to find required mass of fuel to hit targets
        mass_fuel_range = total_mass*(1-np.exp(-aircraft_range*SFC*1*CD_LDmax/(prop_eff*CL_LDmax)))
        mass_fuel_endurance = total_mass - np.exp((np.log(((aircraft_endurance*SFC*CD_minpower)/(prop_eff*sqrt(2*1.225*wing_area)*(CL_minpower**1.5)))+(total_mass*9.81)**(-0.5)))/(-0.5))/9.81

        mass_fuel = max(mass_fuel_range, mass_fuel_endurance)#take the maximum


        V_mindrag = sqrt((total_mass*9.81/(0.5*1.225*wing_area))*sqrt(k/CD0))
        V_maximum = 1.4 * V_mindrag
        # print("Vmax = {}".format(V_maximum))
        CL_Vmax = total_mass*9.81/(0.5*1.225*(V_maximum**2)*wing_area)
        CD_Vmax = CD0 + k * CL_Vmax**2
        # print("CL_Vmax = {}".format(CL_Vmax))
        # print("CD_Vmax = {}".format(CD_Vmax))
        V_stall = sqrt((total_mass*9.81)/(CL_max * 0.5 * 1.225 * wing_area))
        V_minpower = sqrt((total_mass*9.81)/(CL_minpower * 0.5 * 1.225 * wing_area))
        P_Vstall = CD_CLmax * 0.5 * 1.225 * (V_stall**3) * wing_area / prop_eff
        P_Vmax = CD_Vmax * 0.5 * 1.225 * (V_maximum**3) * wing_area / prop_eff
        # print("P_Vstall = {}".format(P_Vstall))
        # print("P_Vmax = {}".format(P_Vmax))
        mass_prop = 20 + 1 * P_Vmax/1000 #simple linear approximation for a piston engine mass using maximum power

        # print(mass_prop)

        # mass_structure = wing_loading_structure * wing_area
        forward_wing_load = []
        forward_wing_mass_dist = []
        forward_wing_shear = [0]
        forward_wing_bend = [0]
        rear_wing_load = []
        rear_wing_mass_dist = []
        rear_wing_shear = [0]
        rear_wing_bend = [0]

        mass_lbs = total_mass * 2.205
        nmax = min(2.1 + (24000/(mass_lbs+10000)),3.8) #find required load factor from FAR 23
        
        y = 0
        array_length = len(rear_wing_y_array)
        while y < array_length:
            forward_wing_load.append(0.5*(forward_wing_cellwidth[y]+forward_wing_cellwidth[y-1])*(forward_wing_cl_array[y]))
            rear_wing_load.append(0.5*(rear_wing_cellwidth[y]+rear_wing_cellwidth[y-1])*(rear_wing_cl_array[y]))
            forward_wing_mass_dist.append(0.5*(forward_wing_cellsparlength[y]+forward_wing_cellsparlength[y-1])*(1-0.0*2*(min(1,abs((sum(forward_wing_cellsparlength[:y])/sum(forward_wing_cellsparlength))-0.5)))))
            rear_wing_mass_dist.append(0.5*(rear_wing_cellsparlength[y]+rear_wing_cellsparlength[y-1])*(1-0.0*2*(min(1,abs((sum(rear_wing_cellsparlength[:y])/sum(rear_wing_cellsparlength))-0.5)))))
            y = y + 1

        # print(nmax)
        # print(sum(forward_wing_mass_dist))
        forward_wing_mass_dist = np.multiply((forward_wing_mass*9.81*nmax)/max(1,(sum(forward_wing_mass_dist))),forward_wing_mass_dist)
        rear_wing_mass_dist = np.multiply((rear_wing_mass*9.81*nmax)/max(1,(sum(rear_wing_mass_dist))),rear_wing_mass_dist)

        # print(sum(forward_wing_mass_dist))

        # fig, ax = plt.subplots()
        # ax.plot(forward_wing_y_array,forward_wing_mass_dist,".")
        # plt.show()

        total_forward_wing_load = sum(forward_wing_load)
        # print(total_forward_wing_load)
        total_rear_wing_load = sum(rear_wing_load)
        # print("total_forward_wing_load = {}".format(total_forward_wing_load))
        # print("total_rear_wing_load = {}".format(total_rear_wing_load))
        total_load = total_forward_wing_load+total_rear_wing_load
        # print("total_load")

        load_dist_adjustment = ((total_mass)*9.81*nmax)/total_load

        forward_wing_maxload = np.multiply(forward_wing_load,load_dist_adjustment)#scale to MTOW and max load factor no more than 3.8G
        # print(sum(forward_wing_maxload))
        rear_wing_maxload = np.multiply(rear_wing_load,load_dist_adjustment)
        forward_wing_maxload_original = forward_wing_maxload
        rear_wing_maxload_original = rear_wing_maxload
        wing_weight_subtraction = 1
        # forward_wing_maxload = np.subtract(forward_wing_maxload,np.multiply(wing_weight_subtraction,forward_wing_mass_dist))#subtract weight of wings at max load factor no more than 3.8 G (assuming constant spanwise distribution)
        # rear_wing_maxload = np.subtract(rear_wing_maxload,np.multiply(wing_weight_subtraction,rear_wing_mass_dist))

        # print("forward_wing_maxload = {}".format(sum(forward_wing_maxload)))
        # print("rear_wing_maxload = {}".format(sum(rear_wing_maxload)))
        total_maxload = sum(forward_wing_maxload)+sum(rear_wing_maxload)

        forward_wing_load_pc = abs(sum(forward_wing_maxload)/total_maxload)
        rear_wing_load_pc = abs(sum(rear_wing_maxload)/total_maxload)
        # print("forward_wing_load_pc = {}".format(forward_wing_load_pc))
        # print("rear_wing_load_pc = {}".format(rear_wing_load_pc))

        forward_wing_reaction = forward_wing_load_pc*total_maxload*np.sign(total_forward_wing_load)
        rear_wing_reaction = rear_wing_load_pc*total_maxload*np.sign(total_rear_wing_load)
        # print("forward_wing_reaction = {}".format(forward_wing_reaction))
        # print("rear_wing_reaction = {}".format(rear_wing_reaction))
        midpoint = round(array_length/2)
        # print(midpoint)

        # forward_wing_maxload[midpoint-1] = forward_wing_maxload[midpoint-1] - forward_wing_reaction/2
        # forward_wing_maxload[midpoint] = forward_wing_maxload[midpoint] - forward_wing_reaction/2
        # rear_wing_maxload[midpoint-1] = rear_wing_maxload[midpoint-1] - rear_wing_reaction/2
        # rear_wing_maxload[midpoint] = rear_wing_maxload[midpoint] - rear_wing_reaction/2
        y = 1

        # print("forward wing max load with reaction = ", sum(forward_wing_maxload))
        forward_wing_spar_mass = [0]
        rear_wing_spar_mass = [0]
        forward_wing_other_mass = [0]
        rear_wing_other_mass = [0]
        while y<array_length:
            forward_wing_shear.append(forward_wing_shear[y-1]+(forward_wing_maxload[y])/np.cos(forward_wing_dihe[y])-(forward_wing_mass_dist[y]*np.cos(forward_wing_dihe[y])))
            rear_wing_shear.append(rear_wing_shear[y-1]+(rear_wing_maxload[y])/np.cos(rear_wing_dihe[y])-(rear_wing_mass_dist[y]*np.cos(rear_wing_dihe[y])))
            y = y + 1
        forward_wing_imaginary_reaction = forward_wing_shear[array_length-1]
        rear_wing_imaginary_reaction = rear_wing_shear[array_length-1]

        # forward_wing_shear[midpoint-1] = forward_wing_shear[midpoint-1] - forward_wing_imaginary_reaction/2
        # forward_wing_shear[midpoint] = forward_wing_shear[midpoint] - forward_wing_imaginary_reaction/2
        # rear_wing_shear[midpoint-1] = rear_wing_shear[midpoint-1] - rear_wing_imaginary_reaction/2
        # rear_wing_shear[midpoint] = rear_wing_shear[midpoint] - rear_wing_imaginary_reaction/2

        forward_wing_shear[midpoint:] = np.subtract(forward_wing_shear[midpoint:],forward_wing_imaginary_reaction)
        rear_wing_shear[midpoint:] = np.subtract(rear_wing_shear[midpoint:],rear_wing_imaginary_reaction)
        y = 1
        # print("array_length = {}".format(array_length))
        while y < array_length:
            forward_wing_bend.append(forward_wing_bend[y-1]+forward_wing_cellsparlength[y-1]*forward_wing_shear[y-1])
            rear_wing_bend.append(rear_wing_bend[y-1]+rear_wing_cellsparlength[y-1]*rear_wing_shear[y-1])
            



            forward_wing_max_thickness = 0.14*forward_wing_chord[y]
            # print(forward_wing_bend)
            forward_wing_max_moment = abs(forward_wing_bend[y])
            
            forward_wing_spar_width = max(0.005,6*forward_wing_max_moment/((forward_wing_max_thickness**2)*250e6))
            
            forward_wing_spar_area = forward_wing_spar_width*forward_wing_max_thickness
            forward_wing_spar_volume = forward_wing_spar_area*forward_wing_cellsparlength[y]
            forward_wing_spar_mass.append(forward_wing_spar_volume*1800)
            forward_wing_other_mass.append(1.15*forward_wing_spar_mass[y]*forward_wing_chord[y])


            rear_wing_max_thickness = 0.14*rear_wing_chord[y]
            rear_wing_max_moment = abs(rear_wing_bend[y])
            
            rear_wing_spar_width = max(0.005,6*rear_wing_max_moment/((rear_wing_max_thickness**2)*250e6))
            rear_wing_spar_area = rear_wing_spar_width*rear_wing_max_thickness
            rear_wing_spar_volume = rear_wing_spar_area*rear_wing_cellsparlength[y]
            rear_wing_spar_mass.append(rear_wing_spar_volume * 1800)
            rear_wing_other_mass.append(1.15*rear_wing_spar_mass[y]*rear_wing_chord[y])

            y = y + 1

        forward_wing_spar_mass = sum(forward_wing_spar_mass)  
        forward_wing_other_mass = sum(forward_wing_other_mass)      
        forward_wing_mass = forward_wing_spar_mass +forward_wing_other_mass#* 2.55#2.55
        rear_wing_spar_mass = sum(rear_wing_spar_mass)
        rear_wing_other_mass = sum(rear_wing_other_mass)
        rear_wing_mass = rear_wing_spar_mass +rear_wing_other_mass#* 2.55#2.55

        # print("fwrd_wing_max_bend = {}".format(forward_wing_max_moment))
        # print("fwrd_wing_max_thck = {}".format(forward_wing_max_thickness))
        # print("forward_wing_spar_width = {}".format(forward_wing_spar_width))
        # print("forward_wing_spar_length = {}".format(sum(forward_wing_cellsparlength)))
        # print("fwrd_wing_mass = {}".format(forward_wing_mass))
        # print("rear_wing_mass = {}".format(rear_wing_mass))
        # print("rear_wing_max_bend = {}".format(rear_wing_max_moment))
        # print("rear_wing_max_thickness = {}".format(rear_wing_max_thickness))
        # print("rear_wing_spar_width = {}".format(rear_wing_spar_width))
        # print("rear_wing_spar_area = {}".format(rear_wing_spar_area))
        # print("forward_wing_mass = {}".format(forward_wing_mass))
        # print("rear_wing_mass = {}".format(rear_wing_mass))
        tail_mass = max(forward_wing_mass/fwng_area_not_projected,rear_wing_mass/rwng_area_not_projected)*(tail_chrd0+tail_chrd1)*0.5*tail_span

        wdiff = rwng_xpos0-fwng_xpos0
        tdiff = tail_xpos-rwng_xpos0
        othermass = total_mass-mass_prop-forward_wing_mass-rear_wing_mass

        fuse_shear_est1= -mass_prop*9.81*nmax
        fuse_bend_est1 = fuse_shear_est1*2.5#estimate at fwng

        fuse_shear_est2= fuse_shear_est1+forward_wing_reaction-(forward_wing_mass*9.81*nmax)
        fuse_bend_est2 = fuse_bend_est1+0.5*wdiff*(fuse_shear_est2)#estimate inbetween wings
                                                                     
        fuse_shear_est3= fuse_shear_est2-(wdiff/(wdiff+tdiff))*othermass*9.81*nmax
        fuse_bend_est3 = fuse_bend_est2+0.5*wdiff*(fuse_shear_est3)#estimate at rwng from front

        # fuse_shear_est4= fuse_shear_est3+rear_wing_reaction-(rear_wing_mass*9.81*nmax)
        fuse_bend_est4 = 0.5*tdiff*(tdiff/(wdiff+tdiff))*othermass*9.81*nmax#estimate at rwng from back


        fuselage_max_moment = max(abs(fuse_bend_est1),abs(fuse_bend_est2),abs(fuse_bend_est3),abs(fuse_bend_est4))
        
        fuselage_spar_area = 2*max(0.02**2,fuselage_max_moment/(fuselage_width*250e6))
        # print(fuselage_spar_area)
        
        # fuselage_spar_area = fuselage_spar_width*(fuselage_width/2)
        fuselage_spar_volume = fuselage_spar_area*fuselage_length
        fuselage_spar_mass = fuselage_spar_volume * 1800
        fuselage_mass = fuselage_spar_mass * 18
        # fuselage_spar_mass = fuselage_spar_mass * 0.7
        # print("fuselage_mass = {}".format(fuselage_mass))
        mass_structure = fuselage_mass + forward_wing_mass + rear_wing_mass + tail_mass
        
        mass_avionics = 0#12
        total_mass_new = mass_payload + mass_fuel + mass_prop + mass_structure + mass_avionics
        mass_landinggear_pc = 0.05
        mass_systems_pc = 0.15
        total_mass_new = total_mass_new/(1-(mass_landinggear_pc+mass_systems_pc))
        mass_empt = total_mass_new-mass_payload-mass_fuel
        err1 = abs(total_mass_new-total_mass)/total_mass
        err2 = abs(total_mass_new-total_mass_minus_one)/total_mass_minus_one
        err = max(err1,err2)
        if x == 0:
            err = tol + 1
        
        # print(err)
        relaxation = 0.1
        total_mass_minus_one = total_mass
        total_mass = (total_mass_new-total_mass)*relaxation + total_mass
        massplot.append(total_mass)

        x = x + 1

        if err > 5 or x == 1000:
            print("mass did not converge!!")

        # fig, ax = plt.subplots()
        # ax.plot(forward_wing_y_array,forward_wing_maxload,".")
        # ax.plot(rear_wing_y_array,rear_wing_maxload,".")
        # plt.show()
        # fig, ax = plt.subplots()
        # ax.plot(forward_wing_y_array,forward_wing_shear,".")
        # ax.plot(rear_wing_y_array,rear_wing_shear,".")
        # plt.show()
        # fig, ax = plt.subplots()
        # ax.plot(forward_wing_y_array,forward_wing_bend,".")
        # ax.plot(rear_wing_y_array,rear_wing_bend,".")
        # plt.show()
        # print(total_mass)
        
        # print(x)
        # uncomment next block for mass breakdown for each mass iteration
        # print("mass_fuel = {}".format(mass_fuel))
        # print("mass_payl = {}".format(mass_payload))
        # print("mass_prop = {}".format(mass_prop))
        # print("fwng_mass = {}".format(forward_wing_mass))
        # print("rwng_mass = {}".format(rear_wing_mass))
        # print("mass_fuse = {}".format(fuselage_mass))
        # print("mass_strc = {}".format(mass_structure))
        # print("mass_totl = {}".format(total_mass))
        # print("=================================")


    # uncomment next block for final mass breakdown and performance
    # print("mass_fuel = {}".format(mass_fuel))
    # print("mass_prop = {}".format(mass_prop))
    # print("mass_fwng = {}".format(forward_wing_mass))
    # print("mass_rwng = {}".format(rear_wing_mass))
    # print("mass_tail = {}".format(tail_mass))
    # print("mass_fuse = {}".format(fuselage_mass))
    # print("mass_sysm = {}".format(mass_systems_pc*total_mass))
    # print("mass_strc = {}".format(mass_structure))
    # print("mass_empt = {}".format(mass_empt))
    # print("mass_totl = {}".format(total_mass))
    # print("=================================")
    # print("LD_m = {}".format(max(LD)))
    # print("V_md = {}".format(V_mindrag))
    # print("V_mp = {}".format(V_minpower))
    # print("V_st = {}".format(V_stall))
    # print("prop = {}".format(mass_prop))
    # print("powr = {}".format(P_Vmax))
    # print("mass = {}".format(total_mass))
    # print("wloa = {}".format(total_mass/wing_area))
    # print("=================================")


    # print("fwng_span = {}".format(fwng_span))
    # print("rwng_span = {}".format(rwng_span))

    # print(mass_fuel_range)
    # print(mass_fuel_endurance)

    # print(P_Vmax)
    # fig, ax = plt.subplots()
    # ax.plot(massplot)
    # plt.show()
    # ax.plot(forward_wing_y_array,forward_wing_maxload)
    # ax.plot(rear_wing_y_array,rear_wing_maxload,".")

    # # # ax.plot(results[alpha.index(max(alpha))]["StripForces"]["forward_wing"]["Yle"],results[alpha.index(max(alpha))]['StripShearMoments']['forward_wing']['Mx/(q*Bref*Sref)'])
    # ax.plot(results[alpha.index(max(alpha))]["StripForces"]["forward_wing"]["Yle"],results[alpha.index(max(alpha))]["StripForces"]["forward_wing"]["c_cl"],".")
    # ax.plot(results[alpha.index(max(alpha))]["StripForces"]["rear_wing"]["Yle"],results[alpha.index(max(alpha))]["StripForces"]["rear_wing"]["c_cl"],".")
    # fig, ax = plt.subplots()
    # ax.plot([fwng_0_le_pnt_y,fwng_1_le_pnt_y,fwng_2_le_pnt_y,fwng_3_le_pnt_y], #wing front view
    #         [fwng_0_le_pnt_z-cbar*fwng_chrd0*np.sin(fwng_inci0),fwng_1_le_pnt_z-cbar*fwng_chrd1*np.sin(fwng_inci1),fwng_2_le_pnt_z-cbar*fwng_chrd2*np.sin(fwng_inci2),fwng_3_le_pnt_z-cbar*fwng_chrd3*np.sin(fwng_inci3)])
    # ax.plot([rwng_0_le_pnt_y,rwng_1_le_pnt_y,rwng_2_le_pnt_y,rwng_3_le_pnt_y],
    #         [rwng_0_le_pnt_z-cbar*rwng_chrd0*np.sin(rwng_inci0),rwng_1_le_pnt_z-cbar*rwng_chrd1*np.sin(rwng_inci1),rwng_2_le_pnt_z-cbar*rwng_chrd2*np.sin(rwng_inci2),rwng_3_le_pnt_z-cbar*rwng_chrd3*np.sin(rwng_inci3)])
    # plt.show()
    # fig, ax = plt.subplots()
    # ax.plot([fwng_0_le_pnt_y,fwng_1_le_pnt_y,fwng_2_le_pnt_y,fwng_3_le_pnt_y], #wing top view
    #         [fwng_0_le_pnt_x+cbar*fwng_chrd0*np.cos(fwng_inci0),fwng_1_le_pnt_x+cbar*fwng_chrd1*np.cos(fwng_inci1),fwng_2_le_pnt_x+cbar*fwng_chrd2*np.cos(fwng_inci2),fwng_3_le_pnt_x+cbar*fwng_chrd3*np.cos(fwng_inci3)])
    # ax.plot([rwng_0_le_pnt_y,rwng_1_le_pnt_y,rwng_2_le_pnt_y,rwng_3_le_pnt_y],
    #         [rwng_0_le_pnt_x+cbar*rwng_chrd0*np.cos(rwng_inci0),rwng_1_le_pnt_x+cbar*rwng_chrd1*np.cos(rwng_inci1),rwng_2_le_pnt_x+cbar*rwng_chrd2*np.cos(rwng_inci2),rwng_3_le_pnt_x+cbar*rwng_chrd3*np.cos(rwng_inci3)])
    # plt.show()
    

    # fig, ax = plt.subplots()
    # ax.plot(forward_wing_y_array,forward_wing_maxload_original,".")
    # ax.plot(rear_wing_y_array,rear_wing_maxload_original,".")
    # plt.show()       
    # fig, ax = plt.subplots()
    # ax.plot(forward_wing_y_array,forward_wing_mass_dist,".")
    # ax.plot(rear_wing_y_array,rear_wing_mass_dist,".")
    # plt.show() 
    # fig, ax = plt.subplots()
    # ax.plot(forward_wing_y_array,forward_wing_maxload,".")
    # ax.plot(rear_wing_y_array,rear_wing_maxload,".")
    # plt.show()
    # fig, ax = plt.subplots()
    # ax.plot(forward_wing_y_array,forward_wing_shear,".")
    # ax.plot(rear_wing_y_array,rear_wing_shear,".")
    # plt.show()
    # fig, ax = plt.subplots()
    # ax.plot(forward_wing_y_array,forward_wing_bend,".")
    # ax.plot(rear_wing_y_array,rear_wing_bend,".")
    # plt.show()

    # print(max(rear_wing_bend))
    # print(total_mass)

    X1 = results[3]['Totals']['Xref']
    Cma1 = results[3]['StabilityDerivatives']['Cma']
    #change ref pt
    ref_pnt = avl.Point(x = 0.5*(fwng_xpos0+rwng_xpos0)+3,
                        y=0, z=0)

    ## put it into avl format

    aircraft = avl.Aircraft(name='aircraft',
                            reference_area=wing_area,
                            reference_chord=wing_mac,
                            reference_span=wing_span,
                            reference_point=ref_pnt,
                            mach=mach,
                            surfaces=[forward_wing, rear_wing, tail_wing])

    ##start session with specified geometry

    simple_case = avl.Case(name='zero_aoa',
        alpha=0,beta=0)
                
    session = avl.Session(geometry=aircraft, cases=[simple_case])


    result = session.run_all_cases()

    X2 = result[1]['Totals']['Xref']
    Cma2 = result[1]['StabilityDerivatives']['Cma']

    #use X and Cma pts to find N0, then take away SM*MAC

    N0 = (Cma1*X2-Cma2*X1)/(Cma1-Cma2)
    SM = 0.05
    # print(N0)
    # print(wing_mac)
    Xcg = N0-SM*wing_mac
    # print(Xcg)
    Ycg = 0
    Zcg = le_pnt.z

    #change ref pt
    ref_pnt = avl.Point(x = Xcg,
                            y=0, z=0)
    
    ## put it into avl format
    
    aircraft = avl.Aircraft(name='aircraft',
                                reference_area=wing_area,
                                reference_chord=wing_mac,
                                reference_span=wing_span,
                                reference_point=ref_pnt,
                                mach=mach,
                                surfaces=[forward_wing, rear_wing, tail_wing])
    
    ##start session with specified geometry
    
    simple_case = avl.Case(name='stability_test',
            alpha=datum-3,beta=0)
                    
    session = avl.Session(geometry=aircraft, cases=[simple_case])

    
    result = session.run_all_cases()

    # print(result[1]['StabilityDerivatives']['Cma'])

    #estimate moments of inertia
    #intended for dynamic stability, now redundant

    Ixx = (1/12)*forward_wing_mass*fwng_span**2 + forward_wing_mass*(fwng_zpos0-Zcg)**2 + (1/12)*rear_wing_mass*rwng_span**2 + rear_wing_mass*(rwng_zpos0-Zcg)**2
    Iyy = (forward_wing_mass*(fwng_xpos0-Xcg)**2 + rear_wing_mass*(rwng_xpos0-Xcg)**2 + mass_prop*((fwng_xpos0 - 2)-Xcg)**2)
    Izz = (1/12)*forward_wing_mass*fwng_span**2 + forward_wing_mass*(fwng_xpos0-Xcg)**2 + (1/12)*rear_wing_mass*rwng_span**2 + rear_wing_mass*(rwng_xpos0-Xcg)**2 + mass_prop*((fwng_xpos0 - 2)-Xcg)**2
    Ixy = 0 #left-right symmetry
    Ixz = forward_wing_mass*fwng_xpos0*fwng_zpos0 + rear_wing_mass*rwng_xpos0*rwng_zpos0
    Iyz = 0 #left-right symmetry

    C0 = Ixx*Izz - Ixz**2
    C1 = Izz/C0
    C2 = Ixz/C0
    C3 = C2*(Ixx - Iyy + Izz)
    C4 = C1*(Iyy - Izz)-C2*Ixz
    C5 = 1/Iyy
    C6 = C5*Ixz
    C7 = C5*(Izz - Ixx)
    C8 = Ixx/C0
    C9 = C8*(Ixx-Iyy) + C2*Ixz

    def Cz(psi):
        return np.array([[np.cos(psi),np.sin(psi),0],[-np.sin(psi),np.cos(psi),0],[0,0,1]])

    def Cy(theta):
        return np.array([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])

    def Cx(phi):
        return np.array([[1,0,0],[0,np.cos(phi),np.sin(phi)],[0,-np.sin(phi),np.cos(phi)]])

    # initial state
    V = 50
    X0 = [V,0,1,0,0,0,0,0,0]
    alpha = X0[2]/X0[0]
    beta = X0[1]/X0[0]
    pb2V = X0[3]*wing_span/(2*V)
    qc2V = X0[4]*wing_mac/(2*V)
    rb2V = X0[5]*wing_span/(2*V)

    Cbe = Cx(-X0[6])@(Cy(-X0[7])@Cz(-X0[8]))
    grav = np.array([0,0,-9.8])
    bodygrav = Cbe @ grav

    CL = result[1]['Totals']['CLtot']
    CD = result[1]["Totals"]["CDtot"]+CD0_components
    CY = result[1]['Totals']['CYtot']
    Cl = result[1]['Totals']['Cltot']
    Cm = result[1]['Totals']['Cmtot']
    Cn = result[1]['Totals']['Cntot']
    Cma= result[1]['StabilityDerivatives']['Cma']
    Clb= result[1]['StabilityDerivatives']['Clb']
    Cnb= result[1]['StabilityDerivatives']['Cnb']
    # print("Cma = ",Cma)
    # print("Clb = ",Clb)
    # print("Cnb = ",Cnb)
    stability = max(Cma, Clb, -Cnb)
    # print(stability)
    return total_mass, stability

#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================