#!/usr/local/bin/ Rscript
library(BioCro)
args <- commandArgs(trailingOnly = TRUE)
biocro_initial_state=as.list(read.delim(args[1]))
biocro_parameters=as.list(read.delim(args[2]))
biocro_weather=read.delim(args[3])
# steady_state_module_names<- list(
#     "soil_type_selector",
#     stomata_water_stress = "stomata_water_stress_linear",
#     "parameter_calculator",
#     "soybean_development_rate_calculator",
#     partitioning_coefficients = "partitioning_coefficient_logistic",
#     "soil_evaporation",
#     "solar_zenith_angle",
#     "shortwave_atmospheric_scattering",
#     "incident_shortwave_from_ground_par",
#     "ten_layer_canopy_properties",
#     canopy_photosynthesis = "ten_layer_c3_canopy",
#     "ten_layer_canopy_integrator",
#     partitioning_growth_calculator = "no_leaf_resp_neg_assim_partitioning_growth_calculator",
#     "senescence_coefficient_logistic"
# )
steady_state_module_names=as.list(t(read.delim(args[4])))
derivative_module_names=as.list(t(read.delim(args[5])))

# Gro_wrapper <- function(biocro_initial_state, biocro_parameters, biocro_weather, steady_state_module_names, derivative_module_names){
  
solver_params <- list(
type = 'boost_rkck54',#'Gro_euler',
output_step_size = 1.0,
adaptive_rel_error_tol = 1e-4,
adaptive_abs_error_tol = 1e-4,
adaptive_max_steps = 200)

# biocro_output <- Gro_solver(biocro_initial_state, biocro_parameters, biocro_weather, steady_state_module_names, derivative_module_names, solver_params)

# biocro_output <- run_biocro(biocro_initial_state, biocro_parameters, biocro_weather, steady_state_module_names, derivative_module_names, solver_params,TRUE)

biocro_output <- with(soybean, {run_biocro(
  biocro_initial_state,
  biocro_parameters,
  biocro_weather,
  direct_modules,
  differential_modules,
  solver_params
)})

# biocro_simulation gro(biocro_initial_state, biocro_parameters, biocro_weather, steady_state_module_names, derivative_module_names, solver_params);
# state_vector_map result = gro.run_simulation()


# reorder result columns to be in alphabetical order
biocro_output <- biocro_output[-1,order(names(biocro_output))]
message(sprintf("Hello %s",args[3]))

#   return(biocro_output)
write.table(biocro_output,"biocro_output2.txt",sep="\t",row.names=FALSE,quote=FALSE)
  
# }