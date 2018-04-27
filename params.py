import numpy as np




def mesh_1D(num_point):
    """
    Creating Mesh Point according
    to number of grid slices
    
    returns respectively latitude in radian, degree 
    and difference of two neighbour latitude
    """
    delta_latt   = np.pi/num_point
    latt         = -np.pi/2 + np.arange(num_point+1)*delta_latt
    latt_degree  = np.rad2deg(latt)
    
    return latt, latt_degree, delta_latt

def orbit_single_star(semi_major_axis ,star_mass ,earth_yr ,e ,phi ,phi_peri ,spin_obliq ,azim_obliq ):
    """
    phi_peri ,spin_obliq ,azim_obliq is not needed!
    """
    # is this needed?
#    orbital_period = np.sqrt(semi_major_axis**3/star_mass)
#    orb_freq = 2.0*np.pi/(orbital_period*earth_yr)
#    h_ang = np.sqrt(star_mass*semi_major_axis*(1.0-e**2))  
    #you may put thsi another way
    [phi_peri, spin_obliq,azim_obliq,phi] = np.rad2deg([phi_peri, spin_obliq, azim_obliq, phi])
    r = semi_major_axis*(1.0-e**2)/(1.0+e*np.cos(phi-phi_peri))
    if phi>2.0*np.pi :
        phi = np.mod(phi,2.0*np.pi) 
    phidot = 2.0*np.pi*(np.sqrt(star_mass*semi_major_axis*(1.0-e**2))/(r*r))/earth_yr
    return phi_peri ,spin_obliq ,azim_obliq ,phi ,phidot ,r



 
def INSOLATION(num_point,q0,Lstar,r,lat, spin_obliq, phi, phi_peri, azim_obliq):
    cos_H = np.ones(num_point+1)
    insol = np.ones(num_point+1)*(q0*Lstar)/(np.pi*r*r)
    sind = -np.sin(spin_obliq)*np.cos(phi-phi_peri - azim_obliq)     
    delta = np.arcsin(sind)
    cosd = np.sqrt(1.0-sind*sind)
    tand = np.tan(delta)
    
    for i in range(num_point+1):
        cos_H[i] = -np.tan(lat[i])*tand
        if abs(cos_H[i])>1.0 :
            cos_H[i] = cos_H[i]/abs(cos_H[i])
        H = np.arccos(cos_H[i])
        insol[i] = insol[i]*(H*np.sin(lat[i])*sind + np.cos(lat[i])*cosd*np.sin(H))
    return insol
    
    
def ALBEDO(num_point,spin_obliq,phi,phi_peri,azim_obliq,lat,f_ocean,T):
    albedo = np.ones(num_point+1)
    cos_H = np.ones(num_point+1)
    sind = -np.sin(spin_obliq)*np.cos(phi-phi_peri - azim_obliq)     
    delta = np.arcsin(sind)
    cosd = np.sqrt(1.0-sind*sind)
    tand = np.tan(delta)
    for i in range(num_point+1):
        cos_H[i] = -np.tan(lat[i])*tand
        if abs(cos_H[i])>1.0 :
            cos_H[i] = cos_H[i]/abs(cos_H[i])
        H = np.arccos(cos_H[i]) 
        if H>0 :
            mu = np.sin(lat[i])*sind + (np.cos(lat[i])*cosd*np.sin(H))/H
        elif H == 0.0 :
            mu = 0.5
        mu *= 0.87
        ao = (0.026/(1.1*mu**1.7+0.065))+0.15*(mu-0.1)*(mu-0.5)*(mu-1)
        fi = max(0.0,(1-np.exp((T[i]-273.0)/10.0)))
        ac = max(0.19,-0.07+8.0E-3*np.arccos(mu)*57.2957799513)
        part1 = (1-fi)*(ao*(0.33)+ac*0.67)
        part2 = fi*(0.62*(0.5)+ac*0.5)
        part3 = (1-fi)*(0.2*(0.5)+ac*0.5)
        part4 = fi*(0.85*(0.5)+ac*0.5)
        albedo[i] = (f_ocean)*(part1+part2)+(1-f_ocean)*(part3 + part4)
    return albedo
        
        
def DIFFUSION(D, rot_period, p, p0):
    return D*rot_period*rot_period*(p/p0)
    
def water_boil(p):
    return -0.055*p**6 + 1.2502*p**5 - 11.111*p**4 + 48.87*p**3-112.46*p**2+144.07*p+31.096 + 273
    

def C_tot(num_point,p, p0, C_atm0,f_ocean,freeze,T):
    C_atm = (p/p0)*C_atm0
    C_land = 1.0E9+C_atm
    C_ocean = 2.1E11+C_atm
    f_land = 1 - f_ocean
    f_ice = np.ones(num_point+1)
    
    for i in range(num_point+1):
        f_ice[i] = 1.0 - np.exp(-(freeze-T[i])/10.0)
        if f_ice[i] < 0 :
            f_ice[i] = 0.0
        if T[i] >= freeze :
            C_ice = 0.0
        elif T[i]>freeze and T>freeze-10.0 :
            C_ice = 53.1E9
        elif T[i] <= freeze - 10.0:
            C_ice = 11.1E9
    return f_land*C_land + f_ocean*(f_ice*C_ice + (1.0-f_ice)*C_ocean)           
            


    
def tau_ir(p, p0, T, freeze):
    return 0.79*(p/p0)*(T/freeze)**3

def OLR(sigma_SB,T,tau_ir):
    return sigma_SB*(T**4)/(1.0 + 0.75*tau_ir)




def Q(insol, albedo, infrared):
    return insol*(1.0-albedo) - infrared
    
def hab(num_point, T,freeze, boil):
    hab = np.zeros(T.size, dtype = np.float64)
    for i in range(num_point+1):
        if  T[i]>=freeze and T[i]<=boil :
            hab[i]=1.0
    return hab