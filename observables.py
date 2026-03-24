"""
TU Delft Computational Physics - Project 1: Molecular Dynamics of Argon
Authors: Kyproula Mitsidi, Konstantinos Pourgourides
February - March 2026

~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

This module hosts all the functions that are necessary for the calculation and plotting of
observables, as well as their errors (when applicable).
"""

import molecular_dynamics_argon as mda
import numpy as np
import matplotlib.pyplot as plt
import random as rand
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

plt.rcParams.update({
    'axes.labelsize': 14,   # all axis titles
    'xtick.labelsize': 14,  # tick labels
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'axes.labelpad': 15
})

def get_after_lrs(obs,lrs):
    """
    Returns any observable with shape(num_tsteps) post rescaling.

    Parameters
    ----------
    obs : np.ndarray
        Observable array with shape(num_tsteps)
    lrs : int
        Last rescaling step

    Returns
    -------
    equil_obs : np.ndarray
        Observation post rescaling
    """

    # CALCULATIONS
    #=============
    equil_obs = obs[lrs:]

    # RETURN
    #=======
    return equil_obs

def autocorrelation_function(A, fit_until, show_plot=False, save_plot=False):
    """
    Applies the autocorrelation function method to calculate errors of time-averages.

    Parameters
    ----------
    A : np.ndarray
        Any time-series we want to apply the autocorrelation function to
    fit_until : int
        The time step until calculations are carried out and fit is performed (well-behaved part)
    show_plot : bool
        Show autocorrelation function fit to find chi or not (default False)
    save_plot : str,bool
        Name of plot in case of save (default False)

    Returns
    -------
    A_error : float
        Error of the time average of observable A
    tau : float
        Exponential fit parameter of the autocorrelation function
    """
    
    def model(t,tau):
        """
        Returns the fit equation  exp(-t/tau) for autocorrelation function.

        Parameters
        ----------
        t : np.ndarray
            Time grid
        tau : float
            Exponential fit parameter

        Returns
        -------
        equation : np.ndarray
            Exponential fit function
        """

        # CALCULATIONS
        #=============
        equation = np.exp(-t/tau)

        # RETURN
        #=======
        return equation
        
    def fit(chi):
        """
        Performs the exponential fit on the autocorrelation function.

        Parameters
        ----------
        chi : np.ndarray
            The autocorrelation function with shape(num_tsteps)

        Returns
        -------
        tau : float
            The exponential fit parameter
        """

        # CALCULATIONS
        #=============
        params, covariance = curve_fit(model, mda.time_grid[:fit_until], chi, p0=[70], bounds=(0.9, np.inf))
        if show_plot:
            fig,ax = plt.subplots(figsize=(10,5),dpi=300)
            #plt.title('Autocorrelation Function')
            plt.xlabel('timestep')
            plt.ylabel(r'$\chi_A(x)$')
            plt.plot(mda.time_grid[:fit_until], chi, label=r'$\chi_A(x)$')
            plt.plot(mda.time_grid[:fit_until], model(mda.time_grid[:fit_until],*params),linestyle='--',label=rf'$exp(-t/ \tau) $, $\tau$ ={params[-1]:.2f}')
            plt.legend()
            plt.tight_layout()
            if save_plot!=False:
                fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
            plt.show(close=True)
        tau = params[-1]

        # RETURN
        #========
        return tau

    # CALCULATIONS
    #=============
    chi = np.zeros((fit_until))
    N = len(A)
    for i,t in enumerate(mda.time_grid[:fit_until]):
        a = (N-i)*sum(A[0:N-i]*A[i:N])
        b = np.sum(A[0:N-i])*np.sum(A[i:N])
        c = (N-i)*np.sum(A[0:N-i]**2) - (np.sum(A[0:N-i]))**2
        d = (N-i)*np.sum(A[i:N]**2) - (np.sum(A[i:N]))**2
        if np.sqrt(c)*np.sqrt(d) != 0:
            chi[i] = (a - b)/(np.sqrt(c)*np.sqrt(d))
        else:
            chi[i] = 0
    tau = int(np.round(fit(chi),decimals=0))
    A_error = np.sqrt(2*tau/N)*np.sqrt(np.mean(A**2) - (np.mean(A))**2)

    # RETURN
    #=======
    return A_error, tau
    

def block_bootstrap(A, n, tau = False, resample = True):
    """
    Implements the block bootsrtap method.
    
    Parameters
    ----------
    A : np.array
        Values of observable over time
    n : int
        Repetitions of resampling method
    tau : float
        Exponential fit parameter of the autocorrelation function, if known (default False)
    resample : bool
         Determines whether the method is performed for 1 or n repetitions
         
    Returns
    --------
    if resample is False:
    
    avgA_per_n[i] : float
        Average of observable for 1 resampling
    avgAsqr_per_n[i] : float
        Average of squared observable for 1 resampling
        
    ---------
    
    if resample is True:
    
    avgA : float
        Average of observable over n resamplings
    sigma : float
        Estimated error of the observable
    """

    # CALCULATIONS
    #=============
    if tau == False:
        tau = int(autocorrelation_function(A,100)[1])
    A = A[::tau]
    avgA_per_n = np.zeros((n))
    avgAsqr_per_n = np.zeros((n))
    
    for i in range(n):
        sample = np.array(rand.choices(A, k=len(A)))
        avgA_per_n[i]  = np.mean(sample)
        avgAsqr_per_n[i] = np.mean(sample**2)
        if not resample: 

            # RETURN
            #=======
            return avgA_per_n[i], avgAsqr_per_n[i]
                
    avgA = np.mean(avgA_per_n)
    avgAsqr = np.mean(avgAsqr_per_n)
    sigma = np.sqrt(avgAsqr - avgA**2)

    # RETURN
    #=======
    return avgA, sigma
    
def plot_temperature_evolution(temp_evolution, lrs, save_plot = False):
    """
    Plots the temperature evolution with relation to the target temperature over time.
    
    Parameters
    ----------
    temperature_evolution : np.ndarray
        Temperature over time with shape(num_tsteps)
    lrs : int
        Last rescaling step
    save_plot : str,bool
        Name of plot in case of save (default False)
    """

    # SETUP CANVAS
    #=============
    fig,ax = plt.subplots(figsize=(10,5),dpi=300)
    ax.set_title(rf'Temperature over time for {mda.num_atoms} atoms and T = {(mda.k_b*mda.temperature/mda.epsilon):.3f}$k_B/\epsilon$')
    ax.set_xlabel('Time Step')
    ax.set_xticks(np.arange(0,mda.num_tsteps+3000,3000))
    ax.set_ylabel(rf'Temperature $(\epsilon / k_B)$')
    ax.set_ylim(bottom=0, top=3*(mda.k_b*mda.temperature/mda.epsilon))

    # PLOT
    #=====
    ax.plot(mda.time_grid,temp_evolution,color='deepskyblue',label=r'$T(t)$',alpha=1)
    ax.plot(mda.time_grid,[mda.k_b*mda.temperature/mda.epsilon]*len(mda.time_grid),'k',label=rf'$T_{{target}}$ = {mda.temperature*mda.k_b/mda.epsilon} $\epsilon / k_B$')
    plt.axvline(x=lrs, color='gray', linestyle='--', label=rf'$t_{{lrs}}$ = {lrs}')
    

    # SAVE (OPTIONAL) AND DISPLAY
    #============================
    plt.tight_layout()
    plt.legend()
    if save_plot != False:
        fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
    plt.show(close=True)

    # RETURN
    #=======
    return None



def plot_energy_evolution(ke, pe, lrs, save_plot = False):
    """
    Plots the kinetic, potential and total energy over time with rescaling
    
    Parameters
    ----------
    ke : np.ndarray
        kinetic energy with shape (num_tsteps)
    pe : np.ndarray
        potential energy with shape (num_tsteps)
    lrs : int
        Last rescaling step
    save_plot : str,bool
        Name of plot in case of save (default False)
    """
    
    # SETUP CANVAS
    #=============
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,5),dpi=300)
    ax.set_title(rf'Energy per atom vs time step for {mda.num_atoms} atoms and T = {(mda.k_b*mda.temperature/mda.epsilon):.3f}$k_B/\epsilon$')
    ax.set_ylabel(r'Energy ($\varepsilon$)')
    ax.set_xlabel('time step')
    ax.set_xticks(np.arange(0,len(mda.time_grid),int(len(mda.time_grid)/10)))

    # PLOT
    #=====
    ax.plot(mda.time_grid,pe/mda.num_atoms,color='orange',label='Potential')
    ax.plot(mda.time_grid,ke/mda.num_atoms,color='blue',label='Kinetic')
    ax.plot(mda.time_grid,(ke + pe)/mda.num_atoms,color='green',label='Total')
    plt.axvline(x=lrs, color='gray', linestyle='--', label=rf'$t_{{lrs}}$ = {lrs}')

    # SAVE (OPTIONAL) AND DISPLAY
    #============================
    plt.legend()
    plt.tight_layout()
    if save_plot != False:
        fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
    plt.show(close=True)

    # RETURN
    #=======
    return None


def plot_initial_velocity_distribution(init_vel, save_plot = False):
    """
    Plots the initial velocity distribution.

    Parameters
    ----------
    init_vel : np.ndarray
        Initial velocities with shape(num_atoms,dim)
    save_plot : str,bool
        Name of plot in case of save (default False)
    """

    # CALCULATIONS
    #=============
    speeds = (mda.beta/mda.alpha)*np.linalg.norm(init_vel, axis=1)
    v = np.linspace(0, speeds.max(), 200)

    # SETUP CANVAS
    #=============
    fig,ax = plt.subplots(1,1,figsize=(10,5),dpi=300)
    ax.set_title(rf'Speed distribution vs Maxwell-Boltzmann for {mda.num_atoms} atoms and T = {(mda.k_b*mda.temperature/mda.epsilon):.3f}$k_B/\epsilon$')
    ax.set_ylim(0,1.1)
    ax.set_xlabel(r'Speed ($\sigma / \beta$)')
    ax.set_ylabel('Probability density')

    # PLOT
    #=====
    ax.hist(speeds, bins=25, density=True, color='darkblue', label='Initial Speed Distribution', alpha=0.5)
    ax.plot(v, v**2 * np.sqrt(2/np.pi)*np.exp(-0.5*v**2),color='black', label=r'Maxwell-Boltzmann $(\frac{2}{\pi})v^2 e^{-\frac{v^2}{2}}$')

    # SAVE (OPTIONAL) AND DISPLAY
    #============================
    plt.legend()
    plt.tight_layout()
    if save_plot != False:
        fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
    plt.show(close=True)

    # RETURN
    #=======
    return None
    
def plot_pair_correlation_function(counts, lrs, save_plot = False):
    """
    Calculates and plots the pair correlation function g(r) with errors.
    
    Parameters
    ----------
    counts : np.ndarray
        Number of pairwise distances r_ij between [r,r+dr] per time step 
        to be used in the pair correlation calculations with shape(num_bins_g,num_tsteps)
    lrs : int
        Last rescaling step
    save_plot : str,bool
        Name of plot in case of save (default False)

    Returns
    --------
    g : np.ndarray
        The pair correlation function g(r)
    """
    
    # CALCULATIONS
    #=============
    L_dimless = mda.L/mda.alpha
    bin_edges = np.linspace(0, L_dimless, mda.num_bins_g+1)
    braket_n = np.mean(counts, axis = 1)
    num_timesteps = mda.num_tsteps - lrs
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dr = bin_edges[1] - bin_edges[0]
    g = L_dimless**3 * braket_n / ((mda.num_atoms - 1) * mda.num_atoms * 4 * np.pi * bin_centers**2 * dr)
    errorbars = np.zeros((mda.num_bins_g))
    for i in range(mda.num_bins_g):
        if g[i]>=0.4:
            errorbars[i] = autocorrelation_function(counts[i, :], 100)[0]
    errorbars = L_dimless**3 * errorbars / ((mda.num_atoms - 1) * mda.num_atoms * 4 * np.pi * bin_centers**2 * dr)
    upper = g + errorbars
    lower = g - errorbars
    
    
    # SETUP CANVAS
    #=============
    fig,ax = plt.subplots(figsize=(10,5),dpi=300)
    ax.set_title(rf'Pair correlation function g(r) for {mda.num_atoms} atoms and T = {(mda.k_b*mda.temperature/mda.epsilon):.3f}$k_B/\epsilon$')
    ax.set_xlabel(rf'r ($\sigma$)')
    ax.set_ylabel('g(r)')

    # PLOT
    #=====
    plt.plot(bin_centers,g,color='k',markersize=5, linewidth=1)
    plt.fill_between(bin_centers, lower, upper, color='red', alpha=0.7)
    plt.axhline(y=1.0, color='r', linestyle='--', label='Ideal Gas', alpha=0.3)

    # SAVE (OPTIONAL) AND DISPLAY
    #============================
    plt.legend()
    plt.tight_layout()
    if save_plot != False:
        fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
    plt.show(close=True)
    
    # RETURN
    #=======
    return g

def get_pressure(virial, lrs):
    """
    Calculates pressure of system with error.
    
    Parameters
    ----------
    virial : np.ndarray
        The virial term per time step to be used in the pressure calculations with shape(num_tsteps)
    lrs : int
        Last rescaling step

    Returns
    ------
    pressure : float
        The pressure of the system
    error_pressure : float
        The error in the pressure of the system
    """
    
    # CALCULATIONS
    #=============
    pressure = 1 - (mda.thermodynamic_beta/(3*mda.num_atoms))*np.mean(mda.epsilon*get_after_lrs(virial,lrs))
    error_virial = autocorrelation_function(mda.epsilon*get_after_lrs(virial,lrs), fit_until=300, show_plot=False)[0]
    error_pressure = mda.thermodynamic_beta/(3*mda.num_atoms)*error_virial
    
    # RETURN
    #=======
    return pressure, error_pressure

def get_specific_heat_capacity(ke, lrs):
    """
    Calculates specific heat capacity of the system by using block bootstrap resampling.
    
    Parameters
    ----------
    ke : np.ndarray
        kinetic energy of the system with shape(num_tsteps)
    lrs : int
        Last rescaling step

    Returns
    ------
    cv : float
        Specific heat capacity of the system
    error_cv : float
        The error in the specific heat capacity of the system
    """

    # CALCULATIONS
    #=============
    n = 1500
    cv_samples = np.zeros((n))
    tau = int(autocorrelation_function(get_after_lrs(ke,lrs),100)[1])
    
    for i in range(n):
        avgA_per_n,  avgAsqr_per_n = block_bootstrap(get_after_lrs(ke,lrs), n, tau, resample = False)
        deltaK2 =  avgAsqr_per_n - avgA_per_n**2  
        cv_samples[i] = (3*avgA_per_n**2) / (2*avgA_per_n**2 - 3*mda.num_atoms*deltaK2)
    
    cv_mean = np.mean(cv_samples)
    cv_error = np.std(cv_samples)
    
    # RETURN
    #=======
    return cv_mean, cv_error
    
def plot_msd(pos, lrs, save_plot=False):
    """
    Calculates mean squared displacement
    
    Parameters
    ----------
    pos : np.ndarray
        positions of atoms with shape(num_atoms,num_tsteps,dim)
    lrs : int
        Last rescaling step
    save_plot : str,bool
        Name of plot in case of save (default False)
        
    Returns
    -------
    D_dimless : float
        The diffusion coefficient
    D_err_dimless : float
        The error in diffusion coefficient
    """   

    # CALCULATIONS
    #=============
    pos_equil = pos[:,lrs:,:]
    lag_grid = np.arange(1,mda.num_tsteps-lrs,10)
    msd_per_lag = np.zeros((len(lag_grid)))
    std_msd_per_lag = np.zeros((len(lag_grid)))
    
    for i, lag in enumerate(lag_grid):
        displacement_per_atom = (pos_equil - np.roll(pos_equil,-lag, axis = 1)) [:,:-lag,:]
        squared_distance_per_atom = np.sum(displacement_per_atom**2, axis = 2)
        msd_per_lag[i] = np.mean(squared_distance_per_atom, axis = (0,1))
        std_msd_per_lag[i] = np.std(squared_distance_per_atom, axis = (0,1))

    # ERROR PROPAGATION
    #==================
    errors_of_D = std_msd_per_lag[-10::]/(6*lag_grid[-10::])
    weights = 1/errors_of_D**2
    sum_weights = np.sum(weights)
    D_dimless_arr = msd_per_lag[-10::]/(6*lag_grid[-10::])
    D_dimless = np.sum(D_dimless_arr*weights)/sum_weights
    D_err_dimless = np.sqrt(1/sum_weights)

    # SETUP CANVAS
    #=============
    fig,ax = plt.subplots(figsize=(10,5),dpi=300)
    ax.set_title(rf'MSD over lag for {mda.num_atoms} atoms and T = {(mda.k_b*mda.temperature/mda.epsilon):.3f}$k_B/\epsilon$')
    ax.set_xlabel('Lag')
    ax.set_ylabel(rf'MSD ($\sigma^2$)')
    ax.set_xticks(np.arange(0,mda.num_tsteps+1e3,1e3))

    # PLOT
    #=====
    plt.plot(lag_grid, msd_per_lag, color='blue', label='MSD for different lag values')

    # SAVE (OPTIONAL) AND DISPLAY
    #============================
    plt.tight_layout()
    plt.legend()
    if save_plot != False:
        fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
    plt.show(close=True)

    return D_dimless, D_err_dimless
    

def plot_position_evolution_animation(pos, save_plot=False):
    """
    Makes an animation of the positions as they evolve in time.

    Parameters
    ----------
    pos : np.ndarray
        positions of atoms with shape(num_atoms,num_tsteps,dim)
    save_plot : str,bool
        Name of plot in case of save (default False)

    Returns
    ------
    ani : gif
        Gif of the animation (save_plot must be True)
    """

    from matplotlib.animation import FuncAnimation

    frame_step = 10
    num_atoms, num_tsteps, dim = pos.shape

    # SETUP CANVAS
    #=============
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, mda.L/mda.alpha)
    ax.set_ylim(0, mda.L/mda.alpha)
    ax.set_zlim(0, mda.L/mda.alpha)

    # ANIMATION
    #==========
    scat = ax.scatter(pos[:,0,0], pos[:,0,1], pos[:,0,2], s=20)
    frames = range(0, num_tsteps, frame_step)
    def update(t):
        scat._offsets3d = (pos[:,t,0], pos[:,t,1], pos[:,t,2])
        ax.set_title(f"timestep = {t}")
        return scat,
    ani = FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=50,
        cache_frame_data=False
    )

    # SAVE (OPTIONAL) AND DISPLAY
    #============================
    if save_plot:
        ani.save(f'journal_plots/{save_plot}.gif', writer='pillow', fps=20)

    # RETURN
    #=======
    return ani