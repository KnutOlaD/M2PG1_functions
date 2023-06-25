'''
Package containing functions that does various calculations for M2PG1 output and/or input data that
does not have the main purpose of producing plots/input/output files.  

Author: Knut Ola Dølven

Content:

volume_to_mass(volume,density,depth,gas='methane'): 
Function that calculates the molar mass of a volume of gas

PresfromDepth(depth,rho)
Function that calculates the hydrostatic pressure from depth and density

get_density(temperature,salinity)
Function that calculates the density of water using the UNESCO 1981 equation of state

get_brsm(bubble_size_range=[0.1,10],
                temperature=5,
                salinity=35,
                model=None,
                tuning=None,
                resolution=100)
Function that uses different bubble rising speed models to return the bubble rising
speed for all bubbles in the bubble size range with resolution resolution.

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

##########################################################################################

def get_density(temperature,salinity,pressure=None):
    '''
    Function that calculates the density of water using the UNESCO 1981 equation of state 
    neglecting the effect of pressure, i.e. simplified equation. Might add in pressure
    dependence later. 

    Inputs:
    temperature: scalar or 1d array of temperature values in degrees Celsius
    salinity : scalar or 1d array of salinity values in g/kg

    Outputs:
    density: scalar or 1d array of density values in kg/m^3
    '''
    
    #check if the inputs are arrays or scalars and make them arrays if they are scalars
    if np.isscalar(temperature):
        temperature = np.array([temperature])
    if np.isscalar(salinity):
        salinity = np.array([salinity])
    if np.isscalar(pressure):
        pressure = np.array([pressure])

    #Define coefficients
    a0 = 999.842594
    a1 = 6.793952*10**-2
    a2 = -9.095290*10**-3
    a3 = 1.001685*10**-4
    a4 = -1.120083*10**-6
    a5 = 6.536332*10**-9

    b0 = 8.24493*10**-1
    b1 = -4.0899*10**-3
    b2 = 7.6438*10**-5
    b3 = -8.2467*10**-7
    b4 = 5.3875*10**-9

    c0 = -5.72466*10**-3
    c1 = 1.0227*10**-4
    c2 = -1.6546*10**-6

    d0 = 4.8314*10**-4

    #Calculate density
    density = a0 + a1*temperature + a2*temperature**2 + \
        a3*temperature**3 + a4*temperature**4 + a5*temperature**5 + \
            (b0 + b1*temperature + b2*temperature**2 + b3*temperature**3 \
             + b4*temperature**4)*salinity + (c0 + c1*temperature + \
                c2*temperature**2)*salinity**1.5 + d0*salinity**2

    return density

##########################################################################################

def get_brsm(bubble_size_range=[0.1,10],
             temperature=5,
             salinity=35,
             model=None,
             tuning=None,
             resolution=100):
    '''
    Function that uses different bubble rising speed models to return the bubble rising 
    speed for all bubbles in the bubble size range with resolution resolution. 
    Supported models are listed below in the Inputs section.

    Inputs: #----------------------------------------------
    bubble_size_range: 1d-array 
        Array containing the minimum and maximum bubble size in mm.
        Default is [0.1,2.5] mm
    temperature: scalar 
        Gives temperature values in degrees Celsius.
        Default is 5 degrees Celsius
    salinity: scalar 
        Salinity in g/kg.
        Default is 35 g/kg
    model: string
        bubble rising speed model.
        Default is the model from Woolf (1993). 
    
        Alternatives are
    
        'Woolf1993' (Woolf, 1993). This is the default model. Linear increase in rising speed with 
        bubble size, then constant at 0.25m s^-1 for all bubbles larger a certain threshold size. 

        'Fan1990' (Fan and Tsuchiyua, 1990), has tuning input [c,d] which determines if the 
        model is for freshwater/saltwater and clean/dirty bubbles. The c parameter is set
        to 1.4 for seawater and 1.2 for pure freshwater. The d parameter is set to 1.6 for completely
        clean bubbles and 0.8 for dirty bubbles. I'm guessing there are valid intermediate
        values for these parameters as well, but see Fan & Tsuchiyua, 1990 and Leifer & Patro, 2002
        for details.

        'Woolf1991' (Woolf and Thorpe, 1991)

        'Leifer2002' (Leifer & Patro, 2002)

        'Veloso2015' (Leifer, 2000, clean bubbles)

    tuning: Some models comes with different parameter settings related to e.g. 
        clean/surfactant covered. See description of bubble rising speed models
        for details. This parameter can be various things.  
    
    resolution: scalar
        Number of points in the bubble size range to calculate the bubble rising speed for.    
        Default is 100 points.
        
    
    Outputs: #--------------------------------------------------
    brsm: 2d-array
        Array containing the bubble sizes in column 1 and the bubble rising speed for all bubbles 
        in the bubble size range

    #end
    '''

    ### Calculate input parameters
    viscosity_sw= get_viscosity(temperature,salinity) #Viscosity of the seawater [Pa s]
    density = get_density(temperature,salinity) #kg/m^3 #density of seawater [kg/m^3]
    #Would be better to use gsw, but using this function to avoid dependencies. 
    g_acc = 9.81 #m/s^2 #Gravitational acceleration [m/s^2]
    surf_tens = 0.0728 #N/m #Surface tension of water is 72 dynes/cm, but some suggest 0.0728 N/m - 
    #guessing this is the surface tension of seawater then. [N/m]
    visc_kinematic = viscosity_sw/density #m^2/s #Kinematic viscosity of seawater [m^2/s]
    #Recalculate the bubble size range to meters
    bubble_size_range = np.array(bubble_size_range)/1000 #[m]

    #Create the bubble array
    b_size_array = np.linspace(bubble_size_range[0],bubble_size_range[1],resolution) #[mm]
    BRSM = np.zeros((resolution,2)) #Bubble rising speed matrix
    BRSM[:,0] = b_size_array #Bubble size array

    #Check what model to use
    if model == None:
        model = 'Woolf1993'
        print('No model chosen, using default model Woolf, 1993')

    if model == 'Woolf1993':
        BRSM[:,1] = 0.172*b_size_array**1.28*g_acc**0.76*visc_kinematic**-0.56
        #Set all rs higher than 0.25 m/s to 0.25 m/s
        BRSM[BRSM[:,1]>0.25,1] = 0.25
    elif model == 'Fan1990':
        if tuning == None:
            print('This model needs tuning parameters, specify tuning parameters [c,d],'
                  'see documentation for details. Using parameters for clean ' 
                  'bubbles in seawater as default.')
            c = 1.4
            d = 1.6
        else:
            c = tuning[0]
            d = tuning[1]
        #Calculate the Morton number
        Morton = (g_acc*viscosity_sw**4)/(density*surf_tens**3)
        BRSM[:,1] = (((density*g_acc*b_size_array**2)/
                        (3.68*(Morton**-0.038)*viscosity_sw))**-d+
                    ((c*surf_tens)/(density*b_size_array)+
                        g_acc*b_size_array)**(-d/2))**(-1/d)

    return BRSM
