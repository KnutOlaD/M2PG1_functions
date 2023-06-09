'''
Package containing functions that does various calculations for M2PG1 output and/or input data that
does not have the main purpose of producing plots/input/output files.  

author: Knut Ola Dølven

Content:

volume_to_mass(volume,density,depth,gas='methane'): 
Function that calculates the molar mass of a volume of gas
'''
#Imports
import numpy as np
import matplotlib.pyplot as plt

#####################################
############ FUNCTIONS ##############
#####################################

def volume_to_mass(volume,temperature,pressure,gas='methane',method='ideal'):
    '''
    Function that calculates the number of moles of a volume of gas using either the ideal gas law or
    van der Waals equation of state. The function can be used to calculate the mass of a gas
    volume at a given depth using a single value or a list of values for volume, density and depth.

    van der Waals equation
    (P+a*n^2/V^2)(V-n*b)=n*R*T  --> n^3(-a*b/V^2)+n^2(a/V)-n(P*b+R*T)+P*V = 0
    where a and b are constants for the gas type, V is the volume, n is the number of moles, R is
    the gas constant, T is the temperature and P is the pressure.

    Inputs: 
    volume: (N,)array_like or scalar
       Volume of gas [m3]
    temperature: (N,)array_like or scalar
         temperature of gas (degC)
    pressure: (N,)array_like or scalar
        pressure [Pa]
    gas: string
        Gas type, default is methane, supported types are 'methane', 'carbondioxide', 'nitrogen',
        'argon' and 'oxygen'
    method: string
        Method for calculating molar mass, options are
          'ideal' - Ideal gas law
          'vdw' - van der Waals equation of state
          'Mario' - Mario's equation of state
        Default is 'ideal'
    
    Outputs:
    moles: (N,)array_like or scalar
        Number of moles of gas [mol]
    '''
    #Check if inputs are arrays or scalars
    if isinstance(volume,(list,np.ndarray)):
        volume = np.array(volume)
        temperature = np.array(temperature)
        pressure = np.array(pressure)
    else:
        volume = np.array([volume])
        temperature = np.array([temperature])
        pressure = np.array([pressure])

    #Check if inputs are of same length
    if not (len(volume) == len(temperature) == len(pressure)):
        raise ValueError('Inputs must be of same length')
    
    ### Set input parameters ###

    #Get temperature to Kelvin
    temperature = 273.15 + temperature

    #Check what gas it is and get the correct gas constants
    if gas == 'methane':
        molar_mass = 16.04 #g/mol
        a = 2.283*10**(-1) #Pa*m^6/mol^2
        b = 4.301*10**(-5) #m^3/mol
        R = 8.314 #J/mol*K
    elif gas == 'carbondioxide':
        molar_mass = 44.01
        raise ValueError('Gas not supported yet')
    elif gas == 'nitrogen':
        molar_mass = 28.01
        raise ValueError('Gas not supported yet')
    elif gas == 'argon':
        molar_mass = 39.95
        raise ValueError('Gas not supported yet')
    elif gas == 'oxygen':
        molar_mass = 32.00
        raise ValueError('Gas not supported yet')
    
    #Set up equation
    if method == 'ideal':
        #Ideal gas law
        moles = pressure*volume/(R*temperature)
        coeff = []
    elif method == 'vdw':
        #van der Waals equation of state
        #Solve the cubic equation of state for n
        #Coefficients for cubic equation. 
        coeff = np.array([-a*b/volume**2,a/volume,-(pressure*b+R*temperature),pressure*volume])
        #Find out if there's more than 1 equation
        if len(coeff.shape) > 1:
            #Then loop over and find all the roots
            moles = np.zeros(len(coeff[0,:]))
            for i in range(len(coeff[0,:])):
                #Solve the cubic equation and extract real root of the third root
                moles[i] = np.real(np.roots(coeff[:,i])[2].astype(float))
        else:
            #Solve the cubic equation and extract real root of the third root
            moles = np.real(np.roots(coeff[:,0])[2].astype(float))
            #Calculate mass from n
            #mass = n*molar_mass
    elif method == 'Mario':
        #Calculate the density of the gas at the given temperature and pressure
        density = pressure*molar_mass/(R*temperature)
        #Dont understand this.. 

    return(moles)
    

###############################################################################################

def PresfromDepth(depth,rho):
    '''
    Function that calculates the hydrostatic pressure (not including atmospheric pressure) from 
    density and depth assuming hydrostatic equilibrium and incompressible fluid.

    Inputs:
    depth: (N,)array_like or scalar
        Depth [m]
    rho: (N,)array_like or scalar  
        Density [kg/m3]

    Outputs:
    pressure: (N,)array_like or scalar
        Pressure [Pa]
                '''
    #Check if inputs are arrays or scalars
    if isinstance(depth,(list,np.ndarray)):
        depth = np.array(depth)
        rho = np.array(rho)
    else:
        depth = np.array([depth])
        rho = np.array([rho])
    
    pressure = rho*9.81*depth

    return(pressure)

######## testrun ########
temperature = 10
pressure = 10132500
volume = 1
gas = 'methane'
method = 'ideal'
mass = volume_to_mass(volume,temperature,pressure,gas,method)

print(mass)

####### testrun with gas at different pressure #########
#taking a gas volume of 1 m3 at 10 degC and 1 atm pressure and
#sinking it gradually down to 1500 m depth
#the gas is assumed to be methane
#the mass of the gas is calculated using the ideal gas law and van der Waals equation of state
#the mass of the gas is then plotted as a function of depth
#Define depth array
depth = np.arange(0,1500,1)
#Define pressure array
pressure = PresfromDepth(depth,1027)
#Define temperature array
temperature = np.ones(len(depth))*10
#Define volume array
volume = np.ones(len(depth))*1
#Calculate mass using ideal gas law
mass_ideal = volume_to_mass(volume,temperature,pressure,gas='methane',method='ideal')
#Calculate mass using van der Waals equation of state
mass_vdw = volume_to_mass(volume,temperature,pressure,gas='methane',method='vdw')

#Plot the results
plt.plot(depth,mass_ideal,label='Ideal gas law')
plt.plot(depth,mass_vdw,label='van der Waals equation of state')
plt.xlabel('Depth [m]')
plt.ylabel('[mol]')
plt.legend()
plt.show()

