# Do the calculations inside an empty list so that temporary variables are not created in .Global.
manihot_esculenta_initial_state = with(list(), {
    datalines =
    "symbol                  value
    Rhizome                  0.99
    Leaf                     0.02
    Stem                     0.99
    Root                     1
    Grain                    0.0001
    soil_water_content       0.32
    LeafLitter               0
    RootLitter               0
    RhizomeLitter            0
    StemLitter               0
    substrate_pool_leaf      0.03       # No substrate pool can be zero with the utilization growth module
    substrate_pool_grain     0.00005    # No substrate pool can be zero with the utilization growth module
    substrate_pool_stem      0.005      # No substrate pool can be zero with the utilization growth module
    substrate_pool_root      0.005      # No substrate pool can be zero with the utilization growth module
    substrate_pool_rhizome   0.00005    # No substrate pool can be zero with the utilization growth module
    TTc                      0"

    data_frame = utils::read.table(textConnection(datalines), header=TRUE)
    values = as.list(data_frame$value)
    names(values) = data_frame$symbol
    values
})
