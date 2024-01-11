#' biocro_fba
#' 
#' A function to run biocro while taking metabolism into account
#' 
#' @return a dataframe
#' 
#' @examples
#' biocro_fba()
# library(R.utils)
# func_timer <- function(FUN){
#     results <- tryCatch({
#         withTimeout({FUN(args)},timeout=100)
#     }, error = function(e){
#         if(grepl("reached elapsed time limit",e$message))
#         FUN(args) else
#         paste(e$message, "EXTRACTERROR")
#     })
#     if(grepl("EXTRACTERROR",results)){
#         print(gsub("EXTRACTERROR","",results))
#         results <- NULL
#     }
#     return(results)
    }
    
biocro_fba <- function(
    initial_values = list(),
    parameters = list(),
    drivers,
    direct_module_names = list(),
    differential_module_names = list(),
    ode_solver = default_ode_solver,
    verbose = FALSE
){
    # Save inputs where they will be expected
    write.table(initial_values,'./../yamls/tempInitialState.txt',sep="\t",row.names=FALSE,quote=FALSE)
    write.table(parameters,'../yamls/parameter_files/soybean_parameters_biocro3_fba.txt',sep="\t",row.names=FALSE,quote=FALSE)
    write.table(direct_module_names,'../yamls/parameter_files/soybean_ss_modules_fba.txt',sep="\t",row.names=FALSE,quote=FALSE)
    write.table(differential_module_names,'../yamls/parameter_files/soybean_deriv_modules_fba.txt',sep="\t",row.names=FALSE,quote=FALSE)
    
    # Run python
    command = "python3"
    path2python = '"inst/python_biocro/fba_wrapper.py"'
    # sysout <- func_timer(FUN=system2,args=paste(command," ",path2python),stdout=TRUE)
    sysout = system2(command,path2python,stdout=TRUE)
    
    if(sysout){
        # Load results from csv
        path2results = 'inst/python_biocro/yamls/BioCro_output_complete.txt'
        results <- read.table(path2results)
        return(results)
    }
}
