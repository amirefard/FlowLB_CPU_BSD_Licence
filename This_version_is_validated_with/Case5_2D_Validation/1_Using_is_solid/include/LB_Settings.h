#ifndef LB_SETTINGS_H_
#define LB_SETTINGS_H_
/// *********************
/// SIMULATION PARAMETERS
/// *********************
// These are the relevant simulation parameters.
// They can be changed by the user.

///===================================================================================================================================== ///
///                                               GENERAL (FLOW FIELD, PARTICLE, THERMAL)                                                ///
///===================================================================================================================================== ///
///Thermal Stencil
#define D2Q5_T
//#define D2Q9_T

///Particle or Not?
//#define Flow_With_Point_Particle
//#define Flow_With_Resolved_Particle
#define Flow_Without_Particle

///Thermal or Not?
#define Flow_With_Thermal_Effect
//#define Flow_Without_Thermal_Effect

///Flow Force Scheme?
#define Shan_Chen
//#define Guo
//#define Luo

///SRT or MRT?
#define SRT
//#define MRT  /// Only with Guo model

///===================================================================================================================================== ///
///                                                       FLOW FIELD BC                                                                  ///
///===================================================================================================================================== ///
/// X BC
//#define X_Periodic

#if !defined X_Periodic
    #define X_Wall
    //#define X_velIn_ZeroGradOut
#endif

/// Y BC
//#define Y_Periodic

#if !defined Y_Periodic
    #define Y_Wall
#endif

/// Pressure BC
//#define Ladd_Pressure      /// NOT active in this version
//#define Zou_He_Pressure  /// NOT active in this version

/// Velocity BC
#define Ladd_Velocity
//#define Zou_He_Velocity  /// NOT active in this version

/// Stationary Wall BC
//#define Ladd_Stationary_Wall      /// NOT active in this version
//#define Zou_He_Stationary_Wall  /// NOT active in this version

///===================================================================================================================================== ///
///                                                 FLOW WITHOUT PARTICLE                                                                ///
///===================================================================================================================================== ///
#if defined Flow_Without_Particle
//#define Read_Geometry_From_File

#define No_curved_boundary_treatment //Curved boundary treatment
//#define Curved_boundary_treatment //Curved boundary treatment

#if defined Curved_boundary_treatment //Important: with stationary solid only.
    #define Bouzidi
    //#define Filippova
    //#define Mei_Luo
#endif

//Solid BC Type (Select only one):
#if defined Flow_With_Thermal_Effect
    #define Solid_Constant_Temperature
    //#define Solid_Adiabatic
    //#define Solid_Conjugate
#endif

#endif //End of Flow_Without_Particle Conditions.

///===================================================================================================================================== ///
///                                                           THERMAL BC, ETC                                                            ///
///===================================================================================================================================== ///
///Thermal Conditions:
#if defined Flow_With_Thermal_Effect

#if !defined X_Periodic
    #define Left_Fixed_Temp
    //#define Left_Adiabatic
    #define Right_Fixed_Temp
    //#define Right_Adiabatic
#endif

#if !defined Y_Periodic
    //#define Bottom_Fixed_Temp
    #define Bottom_Adiabatic
    //#define Top_Fixed_Temp
    #define Top_Adiabatic
#endif

#define Force_Thermal_1
//#define Force_Thermal_2

#define Equilibrium_Thermal_1  //The approach that equilibrium distribution function is defined: like paper Mohamad, A. A., M. El-Ganaoui, and R. Bennacer. "Lattice Boltzmann simulation of natural convection in an open ended cavity." International Journal of Thermal Sciences 48, no. 10 (2009): 1870-1875.
//#define Equilibrium_Thermal_2    //The approach that equilibrium distribution function is defined: Lai, Feng-Hsiang, and Yue-Tzu Yang. "Lattice Boltzmann simulation of natural convection heat transfer of Al 2 O 3/water nanofluids in a square enclosure." International Journal of Thermal Sciences 50, no. 10 (2011): 1930-1941.
#endif///End of Thermal Conditions

///===================================================================================================================================== ///
///                                                     RESOLVED PARTICLE                                                                ///
///===================================================================================================================================== ///
#if defined Flow_With_Resolved_Particle

//#define STATIONARY_CYLINDER
//#define STATIONARY_ELLIPSE
#define MOVING_CYLINDER
//#define MOVING_ELLIPSE

#define SDF
//#define MDF        // See. Kang & Hassan (2011):https://doi.org/10.1002/fld.2304 *** ZL Wang et al. IJMF (2008): https://doi.org/10.1016/j.ijmultiphaseflow.2007.10.004 *** ZL Wang et al. IJHMT (2009) https://doi.org/10.1016/j.ijheatmasstransfer.2009.03.048 *** Kun Luo, Zeli Wang, et al. PRE (2007) https://doi.org/10.1103/PhysRevE.76.066709
const int mdfstep = 5;//Number of MDF steps. Adjust this line.

//#define Without_added_mass
#define With_added_mass

//#define Without_added_heat
#define With_added_heat

//#define Two_Point_delta
#define Four_Point_Delta
//#define Kernel_Delta

//#define Circle_Collision_1     // R. Glowinski et al, J. Comput. Phys., 169 (2001) 363�426.///Z.-G. Feng, E.E. Michaelides,  806(Journal of Computational Physics, 195, (2004), pp.602�628)
//#define Circle_Collision_2   // Wan & Turek, Int. J. Numer. Meth. Fluids, 51 (2006), 531-566./// Liang Wang et al., 814(Commun. Comput. Phys., 13 (2013), pp.1151-1172)
#define Circle_Collision_3   //Feng ZG, Michaelides EE. Journal of Computational Physics,202 (2005), pp.20-51.

//#define Ellipse_Collision_1     // R. Glowinski et al, J. Comput. Phys., 169 (2001) 363�426.///Z.-G. Feng, E.E. Michaelides,  806(Journal of Computational Physics, 195, (2004), pp.602�628)
//#define Ellipse_Collision_2
//#define Ellipse_Collision_3
#define Ellipse_Collision_4

#if defined Flow_With_Thermal_Effect
#define Particle_Fixed_Temperature
//#define Particle_Varying_Temperature
#endif

#define Random_Arrangement_ON
//#define Random_Arrangement_OFF

#endif///End of Resolved Particle Conditions:

///===================================================================================================================================== ///
///                                                        POINT PARTICLE                                                                ///
///===================================================================================================================================== ///
///Point Particle Conditions:
#if defined Flow_With_Point_Particle
#define Without_added_mass
//#define With_added_mass

//#define Two_Point_delta
#define Four_Point_Delta
#endif///End of Point Particle Conditions

///===================================================================================================================================== ///
///                                                    OTHER SETTINGS                                                                    ///
///===================================================================================================================================== ///
///Recovery Conditions
#define Recovery_OFF
//#define Recovery_ON

#define ASCII_Recovery_File
//#define BINARY_Recovery_File

///VTK Type
//#define VTK_ASCII //Paraview file format in ASCII format
#define VTK_BINARY //Paraview file format in BINARY format
//#define CASE //Ensight file format

///Use Residual or Not?
//#define Stop_by_Residual

///=====================================================================================================================================
//////////////////////// EXPLANATION /////////////////////////////////////
////////////////////// DO NOT TOUCH //////////////////////////////////////

// Read_Geometry_From_File::
// Condition of Read_Geometry_From_File is just for Flow_Without_Particle
// In this case you should:
// 1: convert the geometry to a matrix of 0 and 1 by the MATLAB program: image_to_matrix.m .This MATLAB rogram creates geometry.txt
// 2: put the geometry.txt in the same folder as where parameters.dat or fluid_data.dat are created.
// 3: Enter manually Nx and Ny values in ALBORZ program.
// Below condition can be used by different boundary conditions.
// If gravity has a non-zero value then Permeability is calculated and reported and saved to fluid_data.dat
// If you deactive below line you can still define some solid nodes in the geometry by putting is_solid=TRUE (See beginning of Flow_Without_Particle section)

// Force Calculation method.
// Guo: Z. Guo, C. Zheng, B. Shi, Discrete lattice effects on the forcing term in the lattice Boltzmann method, Phys. Rev. E 65 (046308) (2002) 1�6.
// Shan-Chen: Shan, X., Chen, H. (1993), Lattice Boltzmann model for simulating flows with multiple phase and components, Physical Review E, 47, 1815-1819
// Luo: L.-S. Luo, Lattice-Gas Automata and Lattice Boltzmann Equations for Two-Dimensional Hydrodynamics, Ph.D. thesis, Georgia Institute of Technology,1993.
// Guo model is a second order model but has less stability than Shan-Chen in low relaxation times. It is also somewhat slower.

// Select Single Relaxation Time (SRT) or Multi Relaxation Time (MRT)
// -MRT model is JUST adjusted for "Guo" force scheme.

// Boundary types
// Exactly one of the boundary options has to be defined.
// X_Periodic
// - Boundary condition of X-Direction: Periodic inlet and outlet
// X_NOT_Periodic
// - Boundary condition of X-Direction: Not periodic

// Output visualisation file format: VTK or CASE(Ensight)

/// Particle types
// Exactly one of the following options has to be defined.
// STATIONARY_CYLINDER
// - for simulation of, e.g., Karman vortex street
// - the cylinder is kept stationary in space, it does not move
// STATIONARY_ELLIPSE
// DEFORMABLE_CYLINDER
// - for simulation of a moving cylinder.
// - If you choose "Direct_Forcing" then the cylinder CANNOT deform but moves in the fluid.
// DEFORMABLE_ELLIPSE
// - for simulation of a moving ellipse


// SDF: Single-Direct-Forcing.
// MDF: Multi-Direct-Forcing. Do NOT use it when there is NO particle inside the flow.
// use MDF in case of Flow_Without_Thermal_Effect

//For description of added mass refer to page 88 and 89 of Kang PhD Thesis (2010)

// Delta Function choice
// - You can choose Original delta function which was first developed by Kr�ger.NOT available in this version
// - You can choose 2 point function: Eq. 67 Kang thesis (804)
// - You can choose 4 point function: Eq 68 Kang thesis (804)
// - You can choose Kernel function: Eq. 17 Feng (806) Z.-G. Feng, E.E. Michaelides / Journal of Computational Physics 202 (2005) pp.20�51

// Collision model for CIRCLE
// - select collision model for more than 1 moving CIRCULAR particles or particle-wall collision. (Just for Direct forcing).

// Activate "Recovery_ON" if you want to continue the program from a recovery step (for example if the program has stopped due to electercity failure but you do NOT want to start from 1st iteration)
// In case of "Recovery_ON", initial values for fluid and particles are read from "data_fluid_recovery.dat" (See Initialize subroutine)
// Activate "Recovery_OFF" if you want to start the program from the 1st iteration.
// Only ONE of these options must be active at each run.

//Stop_by_Residual: Stop the program when residual falls below a certain value

//Curved boundary treatment
// Bouzidi: M. Bouzidi, M. Firdaouss, P. Lallemand, Momentum transfer of a Boltzmann-lattice fluid with boundaries, Physics of Fluids 13, 3452 (2001);
// Filippova
// Mei_Luo: R. Mei, L. Luo, Wl Shy, An Accurate Curved Boundary Treatment in the Lattice Boltzmann Method, Journal of Computational Physics 155, 307�330 (1999)

#endif
