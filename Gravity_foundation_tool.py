'''
This Python code was designed to provide a preliminary sizing and costing for gravity-based fixed & gravity anchors offshore foundations.'
The input parameters are user-defined. The module uses assumptions and formulae from different guidelines, primarily '
the handbook for marine geotechnical engineering. The user can define the location data, device data, and some design parameters: '
whether the foundation is embedded or not, and whether a perimeter skirt is provided or not.
Design checks are performed for sliding, bearing (drained and undrained), and overturning resistance in both x and y directions. '
The outputs is a CSV file including a design space of foundation sizes and design checks as well as the minimum viable foundation size '
a preliminary CAPEX estimate based on the unit cost specified for the foundation material.

'''



import math
import pandas as pd
from enum import Enum

class soil_types_granularOrcohesive(Enum): 
    granular_soil = 1
    cohesive_soil  = 2
#----------------------------------------------------------------------
# INPUTS
#----------------------------------------------------------------------
# Location data
beta_x_seabedSlope_xdir = 0 # unit: degree, # it is assumed that beta_y_seabedSlope_ydir = 0

rho_w_SeaWaterUnitWeight = 10 # unit: kN/m3
soilType = soil_types_granularOrcohesive.granular_soil

phi_drainedFrictionAngle = 38 # unit: degree
c_drainedCohesion = 0 # unit: kPa
phi_undrainedFrictionAngle = 38 # unit: degree
c_undrainedCohesion = 0 # unit: kPa
s_skirtTip_undrainedShearStrength = 0 # unit: kPa, undrained shear strength at the skirt tip (assuming that a perimeter skirt is specified)
s_avg_undrainedShearStrengthAveragedAlongContactSoil = 0 # unit: kPa, undrained shear strength averaged over the side soil contact zone.

rho_soil_sat_soilSaturatedUnitWeight = 20 # unit: kN/m3 
soil_sensitivity = 1 

# mu_frictionCoefficient = 0.65 # coefficient of friction between the foundation base and soil or between the soil and soil when shear keys cause this type of sliding failure.
mu_frictionCoefficient = 0.78 # with shear key

# Device data
Lmin_minimumFoundationLength = 5 # unit: m, it should be at least 1 m.
Lmax_maximumFoundationLength = 15 # unit: m
dL_foundationLengthIncrement = 0.25 # unit: m
t_foundationThickness = 1.5 # unit:m

rho_conc_concreteUnitWeight = 24 # unit: kN/m3
ConcreteUnitCost = 120 # unit: Euro/m3

zs_skirtDepth = 0.1*Lmin_minimumFoundationLength # if it is 0, means there is not any perimeter skirt/ shear key
# zs_skirtDepth = 0 # if it is 0, means there is not any perimeter skirt/ shear key
de_foundationEmbeddedDepth = 0.1*Lmin_minimumFoundationLength

V_verticalLoad = 0 # unit: kN, 'positive' value for vertical load in gravity direction and 'negative' value for opposite direction of gravity (for anchor gravity foundation)
Hx_horizontalLoad_xdir = 3430.9 # unit: kN
Hy_horizontalLoad_ydir = 0 # unit: kN
Mx_moment_xdir = 0 # unit: kN.m
My_moment_ydir = 0 # unit: kN.m

# Design parameters
gamma_p_soilParametersSafetyFactor = 1.15
gamma_d_soilDesignCapacitySafetyFactor = 1.3 

#----------------------------------------------------------------------
# Defining required functions
#----------------------------------------------------------------------
def degree_to_radian(angle_degree):
    angle_radian = angle_degree*(math.pi/180)
    return angle_radian
def calculate_designSoilParameters(soilPar, safetyFac):
    soilPar_d = soilPar / safetyFac
    return soilPar_d
def calculate_slidingLoad(H,V,W,theta):
    P = H * math.cos(theta) + (V + W) * math.sin(theta)
    return P
def calculate_normalLoad(H,V,W,theta):
    N = (V + W)*math.cos(theta) - H * math.sin(theta)
    return N
def check_designCriteria(f , F, gamma):
    if F / gamma > f: design_check = True
    else: design_check = False
    return design_check
def calculate_soilSlidingCapacity(N,zs,hs,B,A,soilType,c_d,phi_d,rho_buoyant,mu,s_skirtTip,s_avg):
    if soilType == soil_types_granularOrcohesive.granular_soil:
        Rp_passiveSoilResistance = math.tan(math.pi/4 + phi_d/2)\
                **2 * rho_buoyant *zs**2 * B/2 # unit:kN
        Q = mu * N + Rp_passiveSoilResistance
    elif soilType == soil_types_granularOrcohesive.cohesive_soil:
        if zs == 0: # there is no perimeter skirt/ shear key
            Q_ST = min(s_skirtTip * A , mu * N)
        else: 
            Q_ST = s_skirtTip * A + 2 * s_avg * hs * B
        Q_LT = c_d * A + mu * N
        Q = min(Q_ST,Q_LT)
    return Q
def calculate_resultantMomentatBaseCenter(H,V,M,Wf,Ws,t,zs,theta):
    M_rb = M + H *(t + zs) * math.cos(theta) + Wf * (t/2 + zs) * math.sin(theta) + \
        Ws * zs/2 * math.sin(theta) + V * (t + zs) * math.sin(theta)
    return M_rb
def calculate_loadEccentricity(M,V):
    e = abs(M / V)
    return e
def check_foundationMinimumLength_usingEccentricity(e,L):
    if 6 * e >= L:
        minLength_check = False
    else:
        minLength_check = True
    return minLength_check

def calculate_soilBearingCapacity(H,V,ds,L,B,Lprime,Bprime,c,phi,rho_buoyant,s_uz,S_t,theta):
    N_phi = (math.tan(math.pi / 4 + phi / 2)) ** 2
    N_q = math.exp(math.pi*math.tan(phi)) * N_phi
    if phi == 0: N_c = 2 + math.pi
    else: N_c = (N_q-1)/math.tan(phi)
    m = (2+Lprime/Bprime)/(1+Lprime/Bprime) * math.cos(theta)**2 + \
        (2+Bprime/Lprime) /(1+Bprime/Lprime)* math.sin(theta)**2
    if (1-H/(V + Bprime*Lprime *c*1/math.tan(phi))) < 0: # this is added to prevent a complex value for iq
        iq = 0
    else: iq = (1-H/(V + Bprime*Lprime *c*1/math.tan(phi))) ** m
    if phi == 0: ic = 1-m*H/(Bprime*Lprime*s_uz*N_c) 
    else: ic = iq - (1-iq)/(N_c*math.tan(phi))
    if ic < 0: ic = 0
    sc = 1 + Bprime/Lprime*N_q/N_c
    bc = 1
    gc = 1
    dc = 1 + 2*(1-math.sin(phi))**2*math.atan(ds/Bprime)*N_q/N_c
    Kc = ic * sc * dc * bc * gc
    q_c = s_uz * N_c * Kc 

    sq = 1 + Bprime/Lprime*math.tan(phi)
    bq = 1
    gq = 1
    dq = 1 + 2*(1-math.sin(phi))**2*math.atan(ds/Bprime)*math.tan(phi)
    Kq = iq * sq * dq * bq * gq
    fz = 1
    q_q = rho_buoyant*ds*(1+(N_q*Kq-1)*fz)

    if (1-H/(V + Bprime*Lprime *c*1/math.tan(phi))) < 0: # this is added to prevent a complex value for iq
        i_gamma = 0
    else: i_gamma = (1-H/(V + Bprime*Lprime*c*1/math.tan(phi)))**(m + 1)
    s_gamma = 1-0.4*Bprime/Lprime    
    b_gamma = 1
    g_gamma = 1
    d_gamma = 1 
    K_gamma = i_gamma * s_gamma * d_gamma * b_gamma * g_gamma
    N_gamma = 2*(1 + N_q)*math.tan(phi)*math.tan(math.pi/4 + phi/5)
    q_gamma = rho_buoyant*(Bprime/2)*N_gamma*K_gamma*fz

    p = 2*(L + B) # unit: m
    zavg = 0.5 * ds
    delta = phi-5*math.pi/180
    Q_bearingCapacity = Lprime * Bprime*(q_c + q_q + q_gamma) +\
          p*ds*(c/S_t + rho_buoyant*zavg*math.tan(delta)) # unit: kN
    
    return Q_bearingCapacity 
def calculate_stabilityMomentatFoundationCorner(V,Wf,L,theta):
    M_sf = (V + Wf) * math.cos(theta) * L/2
    return M_sf
def calculate_overturningMomentatFoundationCorner(H,V,M,Wf,t,L,theta):
    M_of = M + (H * math.cos(theta) + V * math.sin(theta)+ Wf * math.sin(theta)/2) \
                    * t + H * math.sin(theta) * L/2
    return M_of

#----------------------------------------------------------------------
# CALCULATION LOOP
#----------------------------------------------------------------------
# Required calculation parameters
beta_x_seabedSlope_xdir_radian = degree_to_radian(beta_x_seabedSlope_xdir)
phi_drainedFrictionAngle_radian = degree_to_radian(phi_drainedFrictionAngle)
phi_undrainedFrictionAngle_radian = degree_to_radian(phi_undrainedFrictionAngle)

cdd_designDrainedCohesion = calculate_designSoilParameters(c_drainedCohesion,\
                                                  gamma_p_soilParametersSafetyFactor)
phi_dd_designDrainedFrictionAngle_radian = math.atan(calculate_designSoilParameters(\
                            math.tan(phi_drainedFrictionAngle_radian)\
                                          ,gamma_p_soilParametersSafetyFactor))
cdu_designUndrainedCohesion = calculate_designSoilParameters(c_undrainedCohesion,\
                                                  gamma_p_soilParametersSafetyFactor)
phi_du_designUndrainedFrictionAngle_radian = math.atan(calculate_designSoilParameters(\
                            math.tan(phi_undrainedFrictionAngle_radian)\
                                          ,gamma_p_soilParametersSafetyFactor))

dsoil_sideSoilContactDepth = zs_skirtDepth + de_foundationEmbeddedDepth

# Resizing the foundation, in both x and y directions and check 1. sliding in both x & y direstions, 2. bearing capacity, and 3. overturning moment in both x & y directions
n = int(round((Lmax_maximumFoundationLength-Lmin_minimumFoundationLength)\
    /dL_foundationLengthIncrement,0) + 1) # Total iterations for foundation resizing would be n * n
data = list()
for i in range(n):
    B_foundationWidth = Lmin_minimumFoundationLength + i*dL_foundationLengthIncrement
    
    for j in range(n):
        L_foundationLength = Lmin_minimumFoundationLength + j*dL_foundationLengthIncrement

        # Foundation properties
        A_foundationArea = L_foundationLength*B_foundationWidth
        V_foundationVolume = A_foundationArea*t_foundationThickness
        Wf_foundationWeight = V_foundationVolume*rho_conc_concreteUnitWeight # unit:kN
        cost_foundation = V_foundationVolume*ConcreteUnitCost # unit:Euro
        
        # Calculating submerged weights
        rho_soil_buoyant_soilBuoyantUnitWeight = rho_soil_sat_soilSaturatedUnitWeight - \
                                                rho_w_SeaWaterUnitWeight
        Wf_submergedFoundationWeight = V_foundationVolume * \
            (rho_conc_concreteUnitWeight - rho_w_SeaWaterUnitWeight) # unit:kN
        Ws_submergedSoilPlugWeight = A_foundationArea * zs_skirtDepth \
            * rho_soil_buoyant_soilBuoyantUnitWeight # unit:kN
        Wtot_totalSubmergedWeight = Wf_submergedFoundationWeight + Ws_submergedSoilPlugWeight

        #----------------------------------------------------------------------
        # SLIDING CHECK
        #----------------------------------------------------------------------
        # Calculating critical loads for sliding check, only one case, downslope moment, is critical.
        Px_slidingLoad_xdir = calculate_slidingLoad(Hx_horizontalLoad_xdir,V_verticalLoad,\
                                                    Wtot_totalSubmergedWeight, beta_x_seabedSlope_xdir_radian)
        N_normalLoad = calculate_normalLoad(Hx_horizontalLoad_xdir,V_verticalLoad,\
                                            Wtot_totalSubmergedWeight, beta_x_seabedSlope_xdir_radian)
        if N_normalLoad <= 0: 
            sliding_check_xdir = False 
            sliding_check_ydir = False 
            bearing_check = False
            overturning_check_aroundYdir = False
            overturning_check_aroundXdir = False
            data.append([B_foundationWidth,L_foundationLength,t_foundationThickness,\
                     Wf_foundationWeight,cost_foundation,bearing_check,\
                        sliding_check_xdir,sliding_check_ydir,\
                            overturning_check_aroundYdir,overturning_check_aroundXdir])
            continue

        Py_slidingLoad_ydir = Hy_horizontalLoad_ydir

        # Calculating sliding capacity in both x and y directions
        Q_soilSlidingCapacity_xdir = calculate_soilSlidingCapacity(N_normalLoad,zs_skirtDepth,\
                                    dsoil_sideSoilContactDepth,B_foundationWidth,A_foundationArea,\
                                       soilType,cdd_designDrainedCohesion,phi_dd_designDrainedFrictionAngle_radian,\
                                           rho_soil_buoyant_soilBuoyantUnitWeight, mu_frictionCoefficient\
                                            ,s_skirtTip_undrainedShearStrength,\
                                            s_avg_undrainedShearStrengthAveragedAlongContactSoil)
        Q_soilSlidingCapacity_ydir = calculate_soilSlidingCapacity(N_normalLoad,zs_skirtDepth,\
                                    dsoil_sideSoilContactDepth,L_foundationLength,A_foundationArea,\
                                       soilType,cdd_designDrainedCohesion,phi_dd_designDrainedFrictionAngle_radian,\
                                           rho_soil_buoyant_soilBuoyantUnitWeight, mu_frictionCoefficient\
                                            ,s_skirtTip_undrainedShearStrength,\
                                            s_avg_undrainedShearStrengthAveragedAlongContactSoil)
        
        # Checking sliding in both x and y directions
        sliding_check_xdir = check_designCriteria(Px_slidingLoad_xdir, Q_soilSlidingCapacity_xdir, \
                                          gamma_d_soilDesignCapacitySafetyFactor)
        sliding_check_ydir = check_designCriteria(Py_slidingLoad_ydir, Q_soilSlidingCapacity_ydir, \
                                          gamma_d_soilDesignCapacitySafetyFactor)
        
        #----------------------------------------------------------------------
        # BEARING CHECK
        #----------------------------------------------------------------------
        # Checking bearing capacity for both case 1 and case 2
        for k in range(2):
            if k == 1: # Case 2
                Hx_horizontalLoad_xdir = - Hx_horizontalLoad_xdir 
                My_moment_ydir = - My_moment_ydir
            # Calculating critical loads for bearing check, one extra case, upslope moment, might be critical for bearing check.
            Px_slidingLoad_xdir = calculate_slidingLoad(Hx_horizontalLoad_xdir,V_verticalLoad,\
                                                    Wtot_totalSubmergedWeight, beta_x_seabedSlope_xdir_radian)
            N_normalLoad = calculate_normalLoad(Hx_horizontalLoad_xdir,V_verticalLoad,\
                                            Wtot_totalSubmergedWeight, beta_x_seabedSlope_xdir_radian)
            # Calculating resultant moment at base to be used for finding loads eccentricities
            M_ry_resultantMomentatBaseCenter_aroundYdir = calculate_resultantMomentatBaseCenter(\
            Hx_horizontalLoad_xdir, V_verticalLoad, My_moment_ydir,Wf_submergedFoundationWeight,\
                Ws_submergedSoilPlugWeight,t_foundationThickness,zs_skirtDepth,beta_x_seabedSlope_xdir)
            M_rx_resultantMomentatBaseCenter_aroundXdir = Mx_moment_xdir + Hy_horizontalLoad_ydir\
                *(t_foundationThickness + zs_skirtDepth)
            
            # Calculating eccentricities in x and y direction
            ex_eccentricity_xdir = calculate_loadEccentricity(\
                    M_ry_resultantMomentatBaseCenter_aroundYdir, N_normalLoad)
            ey_eccentricity_ydir = calculate_loadEccentricity(\
                    M_rx_resultantMomentatBaseCenter_aroundXdir, N_normalLoad)
            
            # Check minimum length
            minLength_check = [check_foundationMinimumLength_usingEccentricity\
                            (ex_eccentricity_xdir, L_foundationLength),\
                                check_foundationMinimumLength_usingEccentricity\
                                (ey_eccentricity_ydir,B_foundationWidth)]
            if not all(minLength_check):
                bearing_check = False
            else:
                # Calculating effective lengths
                if (L_foundationLength - 2 * ex_eccentricity_xdir)>=\
                    (B_foundationWidth - 2 * ey_eccentricity_ydir):
                    Lprime_effectiveFoundationLength = L_foundationLength - 2 * ex_eccentricity_xdir
                    Bprime_effectiveFoundationWidth = B_foundationWidth - 2 * ey_eccentricity_ydir
                    if Px_slidingLoad_xdir == 0:
                        theta = math.pi/2
                    else:
                        theta = math.atan(Py_slidingLoad_ydir/Px_slidingLoad_xdir)
                else:
                    Lprime_effectiveFoundationLength = B_foundationWidth - 2 * ey_eccentricity_ydir
                    Bprime_effectiveFoundationWidth = L_foundationLength - 2 * ex_eccentricity_xdir
                    if Py_slidingLoad_ydir == 0:
                        theta = math.pi/2
                    else:
                        theta = math.atan(Px_slidingLoad_xdir/Py_slidingLoad_ydir)     
                # Calculating bearing capacity for both drained and undrained condition
                Q_bearingCapacity_drained = calculate_soilBearingCapacity(math.sqrt(Px_slidingLoad_xdir**2\
                                                + Py_slidingLoad_ydir**2), N_normalLoad\
                                            ,dsoil_sideSoilContactDepth,L_foundationLength,\
                                                B_foundationWidth\
                                                ,Lprime_effectiveFoundationLength,Bprime_effectiveFoundationWidth\
                                                    ,cdd_designDrainedCohesion\
                                    ,phi_dd_designDrainedFrictionAngle_radian,\
                                        rho_soil_buoyant_soilBuoyantUnitWeight\
                                            ,s_skirtTip_undrainedShearStrength,soil_sensitivity,theta)
                Q_bearingCapacity_undrained = calculate_soilBearingCapacity(math.sqrt(Px_slidingLoad_xdir**2\
                                                + Py_slidingLoad_ydir**2), N_normalLoad\
                                            ,dsoil_sideSoilContactDepth,L_foundationLength,\
                                                B_foundationWidth\
                                                ,Lprime_effectiveFoundationLength,Bprime_effectiveFoundationWidth\
                                                    ,cdu_designUndrainedCohesion\
                                    ,phi_du_designUndrainedFrictionAngle_radian,\
                                        rho_soil_buoyant_soilBuoyantUnitWeight\
                                            ,s_skirtTip_undrainedShearStrength,soil_sensitivity,theta)
                Q_bearingCapacity = min(Q_bearingCapacity_drained,Q_bearingCapacity_undrained)
                
                # Checking bearing capacity
                bearing_check = check_designCriteria(N_normalLoad, Q_bearingCapacity\
                                                     , gamma_d_soilDesignCapacitySafetyFactor)
            if bearing_check == False: break
       
        #----------------------------------------------------------------------
        # OVERTURNING MOMENT CHECK
        #----------------------------------------------------------------------
        M_sy_stabilityMoment_aroundYdir = calculate_stabilityMomentatFoundationCorner(\
                             V_verticalLoad, Wf_submergedFoundationWeight,L_foundationLength\
                                ,beta_x_seabedSlope_xdir_radian)
        # M_sx_stabilityMoment_aroundXdir = calculate_stabilityMomentatFoundationCorner(\
        #                      V_verticalLoad, Wf_submergedFoundationWeight,B_foundationWidth ,0)
        M_sx_stabilityMoment_aroundXdir = (V_verticalLoad + Wf_submergedFoundationWeight)\
                                         * B_foundationWidth / 2

        M_oy_overturningMoment_aroundYdir = calculate_overturningMomentatFoundationCorner(\
            Hx_horizontalLoad_xdir, V_verticalLoad,My_moment_ydir,Wf_submergedFoundationWeight\
                ,t_foundationThickness,L_foundationLength, beta_x_seabedSlope_xdir_radian)
        # M_ox_overturningMoment_aroundXdir = calculate_overturningMomentatFoundationCorner(\
        #     Hy_horizontalLoad_ydir, V_verticalLoad,Mx_moment_xdir,Wf_submergedFoundationWeight\
        #         ,t_foundationThickness,B_foundationWidth, 0)
        M_ox_overturningMoment_aroundXdir = Mx_moment_xdir + \
                                           Hy_horizontalLoad_ydir * t_foundationThickness 

        overturning_check_aroundYdir = check_designCriteria(M_oy_overturningMoment_aroundYdir, \
                                    M_sy_stabilityMoment_aroundYdir, gamma_d_soilDesignCapacitySafetyFactor)
        overturning_check_aroundXdir = check_designCriteria(M_ox_overturningMoment_aroundXdir, \
                                    M_sx_stabilityMoment_aroundXdir, gamma_d_soilDesignCapacitySafetyFactor)

        # RESULTS
        data.append([B_foundationWidth,L_foundationLength,t_foundationThickness,\
                     Wf_foundationWeight,cost_foundation,bearing_check,\
                        sliding_check_xdir,sliding_check_ydir,\
                            overturning_check_aroundYdir,overturning_check_aroundXdir])

#----------------------------------------------------------------------
# OUTPUT
#----------------------------------------------------------------------
df = pd.DataFrame(data,columns=['Length (m)','Width (m)','Thickness (m)','Weight (kN)','Cost (Euro)',\
                                'Bearing','Sliding_x','Sliding_y','Overturning_y','Overturning_x'])
print(df)
df.to_csv('SELKIE Gravity Foundation Design.csv')
# print(df.to_string())
