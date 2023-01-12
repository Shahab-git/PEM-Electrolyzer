from cmath import asinh, sqrt
import numpy as np
from numpy import log as ln
import math
            #MODEL DESCRIPTION
#The overall reaction of water splitting is given by
"H2O(l) + electrical energy => H2(g) + 1/2 O2(g)"
#The electrolyser consists of two half-cells of PEM and alkaline electrolysis.
"Anode: H2O(l)  =>1/2 O2(g) + 2H+(aq) + 2e-"
V0 = 1.23                          #voltage(V)
#For the reduction, the following reaction occurs in the cathode by a mercury pool electrode:
"Cathode: 2 H2O(l) + 2 e- => H2(g) + 2 OH-(aq)"
V0 = 0.85                          #voltage(V)
#Also, in the mixing chamber the following reaction occurs:
"2H+(aq) + 2 OH-(aq) => 2 H2O(l)"
V0 = 0.85                          #voltage(V)

#Input Parameters
m_inlet_water=0.1                  #Rate of incoming water(kg/s)
T_inlet_water =353                 #Temperature of incoming water(K)
P_inlet_water=100                  #Pressure of incoming water(kPa)
P_outlet_H2=100                    #Enthalpy of outlet hydrogen(kpa)

#Inlet Water
enthalpy= 285.9                    #kJ/mol
"H_inlet_water=enthalpy(Water,T=T_inlet_water,P=P_inlet_water) "      #Enthalpy of incoming water
"s_inlet_water=entropy(Water,T=T_inlet_water,P=P_inlet_water)"        #Entropy of incoming water

#Outlet Oxygen
T_outlet_O2=T_inlet_water                                             #Temperature of outlet oxygen
P_outlet_O2=P_outlet_H2                                               #Pressure of outlet oxygen
"H_outlet_O2=enthalpy(Oxygen,T=T_outlet_O2,P=P_outlet_O2)"            #Enthalpy of outlet oxygen
"s_outlet_O2=entropy(Oxygen,T=T_outlet_O2,P=P_outlet_O2)"             #Entropy of outlet oxygen

#Outlet Hydrogen
T_outlet_H2=T_inlet_water                                             #Temperature of outlet hydrogen
"H_outlet_H2=Enthalpy(Hydrogen,T=T_outlet_H2,P=P_outlet_H2)"          #Enthalpy of outlet hydrogen
"s_outlet_H2=Entropy(Hydrogen,T=T_outlet_H2,P=P_outlet_H2)"           #Entropy of outlet hydrogen

#Electrochemical relationships
J = 6000
F = 96485.3365                                                        #Faraday's constant (C/mol)
N_dot_H2O_in = 100                                                    #Molar flow rate of entering water
N_dot_H2O_reacted=J / (2 * F)                                         #Molar rate of H2O consumed in the reaction
N_dot_H2O_out = N_dot_H2O_in - N_dot_H2O_reacted                      #Molar flow rate of exiting water
N_dot_H2_out = N_dot_H2O_reacted                                      #Molar outlet flow rate of H2
N_dot_O2_out = 0.5 * N_dot_H2_out                                     #Molar outlet flow rate of O2

#Mass Balance around electrolyzer
NH2O = 1.5                                                            #Flow water input to Electrolyzer (mol/s)
MwH2O = 18                                                            #Molecular Weight of H2O (g/mol)
MwO2 = 32                                                             #Molecular weight of O2 (g/mol)
MwH2 = 2.016                                                          #Molecular weight of H2 (g/mol)
NH2 = NH2O                                                            #All water turned into H2 and O2
MH2 = NH2*MwH2                                                        #Mass flow of H2 (g/s)
NO2 = NH2O/2                                                          #All water turned into H2 and O2
MO2 = NO2*MwO2                                                        #Mass flow of O2 (g/s)
Istack = NH2O*2*F
print("Istack: ",Istack)

#Calculating Voltage reversible & Overpotential
Ncell = 1000                                                                        #Number of cell per stack
Nstack = 100                                                                        #Number of stack needed
Icell = Istack/Ncell/Nstack                                                         #Current per cell
T0 = 273                                                                            #Zero Degree Celsius to Kelvin (K)
Tcell = T0+80                                                                       #Cell Operating Temperature
Vrev = 1.5298-1.5421*1e-03*Tcell+9.523*1e-05*Tcell*ln(Tcell)+9.84*10e-08*Tcell**2   #Voltage reversible calculation
print("Vrev: ",Vrev)
Vact = 0.0514*Icell+0.2098                                                          #Calculating activation overpotential
print("Vact: ",Vact)
Vohm = 0.08*Icell                                                                   #Calculating ohmic overpotential
print("Vohm: ",Vohm)
Vovr = Vrev+Vact+Vohm                                                               #Overall voltage
print("Vovr: ",Vovr)
dV = Vact+Vohm                                                                      #Overpotential Difference
CellEff = Vrev/Vovr*100                                                             #Cell Efficiency Calculation
print("Cell Efficiency: ",CellEff,"%")

#Reversible potential
R = 8.3144621                                                          #Gas Constant (J/molK)
T_PEME_K = 353
V_0 = 1.229 - (8.5e-4) * (T_PEME_K - 298)

#Activation overpotential
E_act_a = 76000
E_act_c = 18000
J_ref_a = 1.7e5
J_ref_c = 4.6e3
k=273
J_0_a = J_ref_a * math.exp(-E_act_a / (R * T_PEME_K))                                   #Anode exchange current density(A/m2)
J_0_c = J_ref_c * math.exp(-E_act_c / (R * T_PEME_K))                                   #Cathode exchange current density(A/m2)
V_act_a=(R*T_inlet_water/F)*(J/(2*J_0_a))                                               #Anode activation overpotentia
V_act_c=(R*T_inlet_water/F)*(J/(2*J_0_c))                                               #Cathode activation overpotential
eta_ohm= J * R
eta_act_c = (R * T_PEME_K / F) * ln(J / (2 * J_0_c) + sqrt(J / (2 * J_0_c ** 2 + 1)))
eta_act_a = (R * T_PEME_K / F) * ln(J / (2 * J_0_a) + sqrt(J / (2 * J_0_a ** 2 + 1)))
V = V_0 + eta_act_a + eta_act_c + eta_ohm
V_0=1.229-(8.5*0.0001//(1/k)*(T_inlet_water-298))+(R*T_inlet_water/F)*ln(sqrt(P_outlet_H2*P_outlet_O2)/P_inlet_water) #Reversible potential voltage

#Ohmic overpotential
lambda_a = 14
lambda_c = 10
L = 50e-6
Deltam = 5e-05                                                                     #Thickness of membrane (m)
CD = 2000
lambda_x=(((lambda_a-lambda_c)/L))+lambda_c                                        #water content at the photon exchange membrane edges
sigma_lambda_x=(0.5139*lambda_x-0.326)*math.exp(1268*(1/303)-(1/T_inlet_water))    #Local ionic PEM conductivity of the membrane
dydx=1/sigma_lambda_x
R_ohm="integral(dydx,x,0,L)"		                                               #The overall ohmic resistance
V_ohm = CD*Deltam/sigma_lambda_x		                                           #Overall ohmic overpotential
print("V_ohm: ",V_ohm)

#Total required potential voltage
N_cells = 1000                                                  #Number of cell per stack
V_cell=V_0+V_act_a+V_act_c+V_ohm		                        #cell required potential voltage
V_stack =V_cell*N_cells		                                    #Stack required potential voltage
Icell = Istack/Ncell/Nstack                                     #Current per cell
print("Icell: ",Icell)

#Electrolyzer performance Calculation
LHV_H2 = 119950                                                 #kJ/kg
W_dot_elec=J*V_stack*1000 		                                #Rate of electric energy input
eta_th_H2=(LHV_H2*N_dot_H2_out)/(W_dot_elec)		            #Hydrogen production efficiency

#Mass flow Calculation"
m_H2_out=2*N_dot_H2_out*0.001		                            #Outlet flow rate of H2"
m_dot_H2O_out=18*N_dot_H2O_out*0.001		                    #Outlet flow rate of H2O"
m_O2_out=32*N_dot_O2_out*0.001		                            #Outlet flow rate of O2"
m_O2_H2O_out=m_O2_out+m_dot_H2O_out		                        #Outlet flow rate of O2 and H2O"
m_dot_H2O_reacted=18*N_dot_H2O_reacted*0.001		            #Reacted flow rate of H2O"

#Alternative Ohmic overpotential caculation
R_ohm2=-10**(-7)*T_inlet_water+4.5*10**(-5)
V_ohm2=J*R_ohm2
print("V_ohm2: ",V_ohm2)