"""
TU Delft Computational Physics - Project 1: Molecular Dynamics of Argon
Authors: Kyproula Mitsidi, Konstantinos Pourgourides
February - March 2026

~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

This module hosts all the functions that are necessary for the simulation to run.
The functions are used to initialize the atom positions and velocities, and to make
the time evolution of the system according to the velocity-Verlet algorithm while 
taking all the Lennard-Jones interactions into account.
"""

import numpy as np


def set_physical_constants(density_input, temperature_input, n_input):
    """
    Declares the global physical constants that are used throughout the simulation. Units are in SI.
    
    Parameters
    ----------
    density_input : float
        Atomic density of the system
    temperature_input : float
        Absolute temperature of system
    n_input : int
        Length of box in lattice constant units (L=n*a)
    """
    
    global density, lattice_constant, n, num_atoms, temperature, k_b, mass, epsilon, sigma, alpha, beta, gamma, thermodynamic_beta

    # PHYSICAL CONSTANTS
    #===================
    density = density_input
    n = n_input
    k_b = 1.38e-23
    mass = 6.63e-26
    epsilon = k_b*119.8
    temperature = temperature_input*epsilon/k_b
    sigma = 3.41e-10
    alpha = sigma
    beta = np.sqrt(mass*sigma**2/(epsilon))
    gamma = mass
    thermodynamic_beta = 1/(k_b*temperature)
    num_atoms = 4*n**3
    lattice_constant = (num_atoms/density)**(1/3)*alpha/n


def set_simulation_constants(num_tsteps_input):
    """
    Declares the global simulation constants that are used throughout the simulation.

    Parameters
    ----------
    num_tsteps_input : int
        Total number of time steps for the simulation
    """    
    
    global dim, num_tsteps, time_grid, h, L, L_dimless, L_lattice, num_bins_g

    #USER INPUTS
    #===========
    num_tsteps = num_tsteps_input
    time_grid = np.arange(0,num_tsteps,1)

    #NON-CUSTOMIZABLE
    #===========
    h = 2e-15
    num_bins_g = 500
    dim = 3
    L = n*lattice_constant    
    
def simulate(rescale_period, rescale_force):
    """
    Molecular dynamics simulation using Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    rescale_period : int
        Number of time steps between each rescaling check
    rescale_force : int
        Number of steps for which the temperature needs to converge in order to break rescaling process, else rescaling is forced

    Returns
    -------
    pos : np.ndarray
        Atom positions with shape(num_atoms,num_tsteps,dim)
    ke : np.ndarray
        Kinetic energy of the system with shape(num_tsteps)
    pe : np.ndarray
        Kinetic energy of the system with shape(num_tsteps)
    temperature_evolution : np.ndarray
        Temperature of the system over time with shape(num_tsteps)
    counts : np.ndarray
        Number of pairwise distances r_ij between [r,r+dr] per time step 
        to be used in the pair correlation calculations with shape(num_bins_g,num_tsteps)
    virial : np.ndarray
        The virial term per time step to be used in the pressure calculations with shape(num_tsteps)
    lrs : int
        The last rescaling time step
    """   

    # TIME EVOLUTION
    #===============
    pos, counts, ke, pe, temperature_evolution, virial, lrs = verlet(rescale_period, rescale_force)

    # RETURN
    #=======
    return pos/alpha, ke/epsilon, pe/epsilon, k_b*temperature_evolution/epsilon, counts, virial/epsilon, lrs


def verlet(rescale_period, rescale_force):
    """
    Time evolution of the system using Verlet's algorithm. 

    Parameters
    ----------
    rescale_period : int
        Number of time steps between each rescaling check
    rescale_force : int
        Number of steps for which the temperature needs to converge in order to break rescaling process

    Returns
    -------
    pos : np.ndarray
        Atom positions with shape(num_atoms,num_tsteps,dim)
    counts : np.ndarray
        Number of pairwise distances r_ij between [r,r+dr] per time step 
        to be used in the pair correlation calculations with shape(num_bins_g,num_tsteps)
    ke : np.ndarray
        Kinetic energy of the system with shape(num_tsteps)
    pe : np.ndarray
        Kinetic energy of the system with shape(num_tsteps)
    temperature_evolution : np.array
        Temperature of the system over time with shape(num_tsteps)
    virial : np.ndarray
        The virial term per time step to be used in the pressure calculations with shape(num_tsteps)
    lrs : int
        The last rescaling time step
    """   

    # VARIABLES
    #==========
    rescaling = True
    lrs = num_tsteps-1

    # MEMORY ALLOCATION
    #==================
    pos = np.zeros((num_atoms, num_tsteps, dim))
    ke = np.zeros(num_tsteps)
    pe = np.zeros(num_tsteps)
    temperature_evolution = np.zeros(num_tsteps)
    forces_i = np.zeros((num_atoms, dim))
    forces_iplus1 = np.zeros((num_atoms, dim))
    rel_dist = np.zeros((num_atoms,num_atoms-1))
    du_dr = np.zeros((num_atoms, num_atoms-1))
    virial = np.zeros((num_tsteps))
    counts = np.zeros((num_bins_g, num_tsteps))
    bin_edges = np.linspace(0, L, num_bins_g + 1)
    bin_width = bin_edges[1]
    
    # INITIALIZATION
    #===============
    pos[:,0,:] = init_position()
    vel = init_velocity()
    rel_pos_current, rel_dist = atomic_distances(pos[:,0,:])
    forces_i, du_dr = lj_force(rel_pos_current, rel_dist)
    virial[0] = 0.5*np.sum(rel_dist*du_dr, axis=(0,1))
    ke[0] = kinetic_energy(vel)
    pe[0] = potential_energy(rel_dist)
    temperature_evolution[0] = 2*ke[0] / (3*(num_atoms-1)*k_b)
    
    # VELOCITY-VERLET TIME EVOLUTION
    #===============================
    for i in range(num_tsteps-1):

        # UPDATE OF PHYSICAL QUANTITIES
        #=============================
        pos[:,i+1,:] = (pos[:,i,:] + vel*h + 0.5*(h**2)*(forces_i/mass)) % L
        rel_pos_current, rel_dist = atomic_distances(pos[:,i+1,:])
        forces_iplus1, du_dr = lj_force(rel_pos_current, rel_dist)
        virial[i+1] = 0.5*np.sum(rel_dist*du_dr, axis=(0,1))
        vel += (0.5*h/mass)*(forces_i + forces_iplus1)
        forces_i = forces_iplus1
        bin_idx = ((rel_dist) // bin_width).astype(int)
        counts[:, i] = np.bincount(bin_idx.ravel(), minlength=num_bins_g)

        # VELOCITY RESCALING CHECK
        #=========================
        if rescaling:
            if i>(rescale_period-1) and i%rescale_period == 0:
                vel, rescaling =  rescale(vel, rescale_force, temperature_evolution, i+1)
                if not rescaling: 
                    lrs = i
        
        # ENERGIES AND TEMPERATURE POST VELOCITY RESCALING
        #=================================================
        ke[i+1] = kinetic_energy(vel)
        pe[i+1] = potential_energy(rel_dist)
        temperature_evolution[i+1] = 2*ke[i+1] / (3*(num_atoms-1)*k_b)

    # RETURN
    #=======                                            
    return pos, counts, ke, pe, temperature_evolution, virial, lrs


def rescale(vel, rescale_force, temperature_evolution, current_moment, threshold=5):
    """
    Rescales velocity and checks whether desired condition is reached to terminate the iterations.

    Parameters
    ----------
    vel : np.ndarray
        The velocities of the atoms at a specific time step with shape(num_atoms,dim)
    rescale_force : int
        Number of steps for which the temperature needs to converge in order to break rescaling process
    temperature_evolution : np.ndarray
        The temperature over time with shape (num_tsteps)
    current_moment : int
        The time step for which rescaling is applied

    Returns
    -------
    vel : np.ndarray
        Rescaled velocities at a specific time step with shape(num_atoms, dim)
    rescaling: bool
        Indicates the termination of rescaling when False
    """   

    # CALCULATIONS
    #=============
    rescaling = True
    vel -= np.mean(vel, axis=0)
    factor = np.sqrt((num_atoms-1)*3*k_b*temperature / np.sum(mass*vel**2, axis=(0,1)))
    vel *= factor
    temperature_current = 2*0.5*np.sum(mass*vel**2, axis=(0,1)) / (3*(num_atoms-1)*k_b)
    
    # RESCALING CHECK
    #================
    if current_moment >= rescale_force:
        recent_temp = temperature_evolution[current_moment-rescale_force:current_moment]
        if np.all(np.abs(recent_temp - temperature) < threshold) and abs(temperature_current - temperature) < threshold:
            rescaling = False
        
    # RETURN
    #=======
    return vel, rescaling

    
def atomic_distances(pos):
    """
    Calculates relative positions and distances between particles.

    Parameters
    ----------
    pos : np.ndarray
        The positions of the atoms at a specific time step with shape(num_atoms,dim)
        
    Returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles at a specific time step with shape(num_atoms,num_atoms-1,dim)
    rel_dist : np.ndarray
        The distance between particles at a specific time step with shape(num_atoms,num_atoms-1)
    """
    
    # CALCULATIONS (INCLUDING PERIODIC BCs)
    #=====================================
    rel_pos = pos[:, None, :] - pos[None, :, :]
    bigger_distance_pos = rel_pos > L/2
    bigger_distance_neg = rel_pos < -L/2
    rel_pos[bigger_distance_pos] = -(L - rel_pos[bigger_distance_pos])
    rel_pos[bigger_distance_neg] = -(-L - rel_pos[bigger_distance_neg])

    # REMOVAL OF DIAGONAL ELEMENTS
    #=============================
    i, j = np.where(~np.eye(num_atoms, dtype=bool))
    rel_dist = np.linalg.norm(rel_pos, axis=2)[i, j].reshape(num_atoms, num_atoms-1)
    rel_pos = rel_pos[i, j].reshape(num_atoms, num_atoms-1, 3)

    # RETURN
    #=======
    return rel_pos, rel_dist


def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    net_force : np.ndarray
        The net force acting on particle i due to all other particles at a specific time step with shape(num_atoms,dim)
    du_dr : np.ndarray
        Derivative of LJ potential at specific time step to be used in the virial term with shape(num_atoms,num_atoms-1)
    """

    # CALCULATIONS
    #=============
    du_dr = 4*epsilon*((6*sigma**6/rel_dist**7) - (12*sigma**12/rel_dist**13))  
    net_force = np.sum(-du_dr[:, :, np.newaxis]* rel_pos / rel_dist[:, :, np.newaxis], axis=1) 

    # RETURN
    #=======
    return net_force, du_dr

def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system at a specific time step.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particles for a specific timestep with shape(num_atoms,dim)

    Returns
    -------
    ke : float
        The total kinetic energy of the system for a specific time step
    """

    # CALCULATIONS
    #============
    ke = np.sum(0.5*mass*vel**2, axis=(0,1))

    # RETURN
    #=======
    return ke


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system at a specific time step.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    u : float
        The total potential energy of the system for a specific time step
    """

    # CALCULATIONS
    #=============
    u = 0.5*np.sum(4*epsilon*(sigma**12/rel_dist**12 - sigma**6/rel_dist**6), axis=(0,1))

    # RETURN
    #=======
    return u


def init_velocity():
    """
    Initializes the system with Gaussian distributed velocities.
    
    Returns
    -------
    initial_vel : np.ndarray
        Array of particle velocities with shape(num_atoms,dim)
    """

    # CALCULATIONS
    #=============
    sigma_normal = np.sqrt(k_b*temperature/mass)
    init_vel = np.random.normal(loc = 0.0, scale = sigma_normal, size = (num_atoms,dim)) 

    # RETURN
    #=======
    return init_vel 


def init_position():
    """
    Initializes the system with atoms on an fcc lattice.
    
    Returns
    -------
    initial_pos : np.ndarray
        Array of particle positions with shape(num_atoms,dim)
    """

    # CALCULATIONS
    #=============
    basis = np.array([[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
    init_pos = []

    for i in range(n):
        for j in range(n):
            for k in range(n):
                cell_origin = np.array([i,j,k])*lattice_constant 
                for b in basis:
                    init_pos.append(cell_origin + b*lattice_constant)
    init_pos = np.array(init_pos)

    # RETURN
    #=======
    return init_pos
