

def write_pflotran_steady_flow(filename, vol_flow_rate, diffusion_coef, porosity):

    pflotran_file = f""" 
#================================================

SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE RICHARDS
    /
  /
  CHECKPOINT
    FORMAT HDF5
  /

END
SUBSURFACE


# #=========================== numerical methods ================================
# NUMERICAL_METHODS FLOW
#   LINEAR_SOLVER
#    SOLVER DIRECT
#   /
# END

#=========================== discretization ===================================
GRID
  TYPE unstructured_explicit full_mesh_vol_area.uge 
  GRAVITY 0.d0 0.d0 0.d0
END

#=========================== fluid properties =================================
FLUID_PROPERTY
  DIFFUSION_COEFFICIENT {diffusion_coef:0.2e}
END

DATASET Permeability
  FILENAME dfn_properties.h5
END


#=========================== material properties ==============================
MATERIAL_PROPERTY rock
  ID 1
  POROSITY {porosity:0.2f}
  TORTUOSITY 1.d0
  CHARACTERISTIC_CURVES curve_1   
  ROCK_DENSITY 2305.15   # Gypsum Density in kg/m3
  PERMEABILITY
    DATASET Permeability
  /
/
#================================================= relative permeability-saturation function curve properties =============================
CHARACTERISTIC_CURVES curve_1
  SATURATION_FUNCTION VAN_GENUCHTEN    
  LIQUID_RESIDUAL_SATURATION 0.1d0
  M 0.8d0
  ALPHA 1.d-4
  /
  PERMEABILITY_FUNCTION MUALEM_VG_LIQ  
  LIQUID_RESIDUAL_SATURATION 0.d0     
  M 0.8d0
  /
/


#=========================== output options ===================================
OUTPUT
  PRINT_PRIMAL_GRID
  ACKNOWLEDGE_VTK_FLAW
  FORMAT VTK
  MASS_FLOWRATE
  MASS_BALANCE
  VARIABLES
    LIQUID_PRESSURE
    PERMEABILITY
  /
END


#=========================== times ============================================
TIME
  #INITIAL_TIMESTEP_SIZE  1.d-8 d
  FINAL_TIME 1.d2 y
  MAXIMUM_TIMESTEP_SIZE 10.d0 y
END

#=========================== regions ==========================================
REGION All
  COORDINATES
    -1.d20 -1.d20 -1.d20
    1.d20 1.d20 1.d20
  /
END 

REGION pinned 
  FILE pinned.ex
END

REGION injection_side
  FILE boundary_left_w.ex
END

REGION production_side
  FILE boundary_right_e.ex
END


#=========================================== flow conditions (flow of water) ============================================================
FLOW_CONDITION Initial_condition
  TYPE
  LIQUID_PRESSURE DIRICHLET    
  TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE 1.01325d6
  TEMPERATURE 25 C             
END

FLOW_CONDITION pinned
  TYPE
    LIQUID_PRESSURE dirichlet
  /
  LIQUID_PRESSURE 1.01325d6
END

FLOW_CONDITION injection
    TYPE
        RATE SCALED_VOLUMETRIC_RATE VOLUME
    /
    RATE {vol_flow_rate} m^3/d
END

FLOW_CONDITION extraction 
    TYPE
        RATE SCALED_VOLUMETRIC_RATE VOLUME 
    /
    RATE -{vol_flow_rate} m^3/d
END


#=========================== condition couplers ===============================
INITIAL_CONDITION initial
  FLOW_CONDITION Initial_condition
  REGION All
END

SOURCE_SINK inflow
    FLOW_CONDITION injection 
    TRANSPORT_CONDITION inject
    REGION injection_side
END

SOURCE_SINK outflow
    FLOW_CONDITION extraction 
    TRANSPORT_CONDITION inject  
    REGION production_side 
END


BOUNDARY_CONDITION PINNED 
  FLOW_CONDITION pinned
  TRANSPORT_CONDITION initial  
  REGION pinned
END


#=========================== stratigraphy couplers ============================
STRATA
  REGION All
  MATERIAL rock
END


END_SUBSURFACE

"""

    with open(filename, "w") as fp:
        fp.write(pflotran_file)



def write_pflotran_rxn(filename, ijob, vol_flow_rate, diffusion_coef, porosity, gypsum_rate_const, calcite_rate_constant, gypsum_surface_area, calcite_surface_area):

    pflotran_file = f""" 

SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE RICHARDS
    /
    SUBSURFACE_TRANSPORT transport
      MODE GIRT
      OPTIONS
        SKIP_RESTART
      /
    /
  /
  RESTART
    FILENAME dfn_flow_{ijob:02d}-restart.h5
    RESET_TO_TIME_ZERO
  END
END

SUBSURFACE

#=========================== numerical methods ================================
#NUMERICAL_METHODS FLOW
#  NEWTON_SOLVER
#    ATOL 1.d-8
#    RTOL 1.d-4
#  END
#END
#
#NUMERICAL_METHODS TRANSPORT
#  NEWTON_SOLVER
#    ATOL 1.d-8
#    RTOL 1.d-4
#  END
#END

#=========================== discretization ===================================
GRID
  TYPE unstructured_explicit full_mesh_vol_area.uge 
  GRAVITY 0.d0 0.d0 0.d0
END

#=========================== fluid properties =================================
FLUID_PROPERTY
  DIFFUSION_COEFFICIENT {diffusion_coef:0.2e}
END

DATASET Permeability
  FILENAME dfn_properties.h5
END


#=========================== material properties ==============================
MATERIAL_PROPERTY rock
  ID 1
  POROSITY {porosity:0.2f}
  TORTUOSITY 1.d0
  CHARACTERISTIC_CURVES curve_1   
  ROCK_DENSITY 2305.15   # Gypsum Density in kg/m3
  PERMEABILITY
    DATASET Permeability
  /
/
#================================================= relative permeability-saturation function curve properties =============================
CHARACTERISTIC_CURVES curve_1
SATURATION_FUNCTION VAN_GENUCHTEN    
LIQUID_RESIDUAL_SATURATION 0.1d0
M 0.8d0
ALPHA 1.d-4
/
PERMEABILITY_FUNCTION MUALEM_VG_LIQ  
LIQUID_RESIDUAL_SATURATION 0.d0     
M 0.8d0
/
/
#=========================== chemistry ========================================
CHEMISTRY
  PRIMARY_SPECIES
    H+
    CO3--
    Ca++
    Na+ 
    SO4--
  /
  SECONDARY_SPECIES
    OH-
    HCO3-
    CO2(aq)
    CaCO3(aq)
    CaHCO3+
    CaOH+
    NaCO3-
    CaSO4(aq)
    HSO4-
    NaSO4-
    H2SO4(aq)
    NaHCO3(aq)
  /
  PASSIVE_GAS_SPECIES
    CO2(g)
  /

 MINERALS     # List all primary minerals first, before listing secondary minerals
     Gypsum
     Calcite
  /
  MINERAL_KINETICS
    Gypsum
      # RATE_CONSTANT 1.2d-8 mol/m^2-sec 
      RATE_CONSTANT {gypsum_rate_const:0.2e} mol/m^2-sec 
    /
    Calcite
      # RATE_CONSTANT 2.43d-5 mol/m^2-sec
      RATE_CONSTANT {calcite_rate_constant:0.2e} mol/m^2-sec
    /
  /
  
  DATABASE /Users/jhyman/src/pflotran/database/hanford.dat 
  LOG_FORMULATION
  ACTIVITY_COEFFICIENTS TIMESTEP
  #UPDATE_MINERAL_SURFACE_AREA
  UPDATE_POROSITY
  UPDATE_PERMEABILITY
  #UPDATE_TORTUOSITY
  MOLAL
  OUTPUT
    MINERALS
    PH
    TOTAL
    ALL
    MINERALS
    MINERAL_VOLUME_FRACTION
    MINERAL_SURFACE_AREA
    MINERAL_SATURATION_INDEX
  /
/

#=========================== output options ===================================
OUTPUT
  MASS_BALANCE
  MASS_BALANCE_FILE
    TOTAL_MASS_REGIONS
      All
    /
  /
   FORMAT VTK
   PERIODIC TIME 1.d1 y
   ACKNOWLEDGE_VTK_FLAW
  VARIABLES
    LIQUID_PRESSURE
    POROSITY
    PERMEABILITY
    LIQUID_SATURATION
    NO_FLOW_VARIABLES
    NO_ENERGY_VARIABLES
    VOLUME
/
END


#=========================== times ============================================
TIME
  # INITIAL_TIMESTEP_SIZE  1.d-2 d
  FINAL_TIME 2.d2 y
  MAXIMUM_TIMESTEP_SIZE 1.d1 y
END



#=========================== regions ==========================================
REGION All
  COORDINATES
    -1.d20 -1.d20 -1.d20
    1.d20 1.d20 1.d20
  /
END 

REGION pinned 
  FILE pinned.ex
END


REGION injection_side
  FILE boundary_left_w.ex
END

REGION production_side
  FILE boundary_right_e.ex
END

#=========================== constraints ======================================
CONSTRAINT initial_constraint
  CONCENTRATIONS
    H+     7.07d0    pH   
    Na+    2.d-10    T
    CO3--  1.d-10    T 
    Ca++   15.648d-3 T 
    SO4--  15.648d-3 M Gypsum 
  /
  MINERALS
    Gypsum  {1 - porosity:0.2f}  {gypsum_surface_area}   m^2/kg     # SSA = 0.75 m2/g 
    Calcite 0.d0    {calcite_surface_area} m^2/m^3   # source: Xu & Jiang
  /
# EQUILIBRATE_AT_EACH_CELL
END


CONSTRAINT injection_constraint
  CONCENTRATIONS
    H+     11.39      pH   
    Na+    1.d0       T
    CO3--  0.5d0      T 
    Ca++   1.d-10     T 
    SO4--  1.d-10     T
  /
  END

#=========================== transport conditions =============================
TRANSPORT_CONDITION initial_conc
  TYPE ZERO_GRADIENT
  CONSTRAINT_LIST
    0.d0 initial_constraint
  /
END


TRANSPORT_CONDITION inject
  TYPE dirichlet_zero_gradient
    CONSTRAINT_LIST
    0.d0 injection_constraint
  /
END



#=========================================== flow conditions (flow of water) ============================================================
FLOW_CONDITION Initial_condition
  TYPE
  LIQUID_PRESSURE DIRICHLET    
  TEMPERATURE DIRICHLET
  /
  LIQUID_PRESSURE 1.01325d6
  TEMPERATURE 25 C             
END

FLOW_CONDITION pinned
  TYPE
    LIQUID_PRESSURE dirichlet
  /
  LIQUID_PRESSURE 1.01325d6
END

FLOW_CONDITION injection
    TYPE
        RATE SCALED_VOLUMETRIC_RATE VOLUME
    /
    RATE {vol_flow_rate} m^3/d
END

FLOW_CONDITION extraction 
    TYPE
        RATE SCALED_VOLUMETRIC_RATE VOLUME 
    /
    RATE -{vol_flow_rate} m^3/d
END


#=========================== condition couplers ===============================
INITIAL_CONDITION initial
  FLOW_CONDITION Initial_condition
  TRANSPORT_CONDITION initial_conc
  REGION All
END


SOURCE_SINK inflow
    FLOW_CONDITION injection 
    TRANSPORT_CONDITION inject
    REGION injection_side
END

SOURCE_SINK outflow
    FLOW_CONDITION extraction 
    TRANSPORT_CONDITION inject  
    REGION production_side 
END


INITIAL_CONDITION All
  FLOW_CONDITION Initial_condition
  TRANSPORT_CONDITION initial_conc
  REGION All
END

BOUNDARY_CONDITION PINNED 
  FLOW_CONDITION pinned
  TRANSPORT_CONDITION initial_conc  
  REGION pinned
END

#=========================== stratigraphy couplers ============================
STRATA
  REGION All
  MATERIAL rock
END

END_SUBSURFACE

"""

    with open(filename, "w") as fp:
        fp.write(pflotran_file)
