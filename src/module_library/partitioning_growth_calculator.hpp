#ifndef PARTITIONING_GROWTH_CALCULATOR_H
#define PARTITIONING_GROWTH_CALCULATOR_H

#include "../modules.h"
#include "../state_map.h"
#include "BioCro.h"  // for resp

/**
 *  @class partitioning_growth_calculator
 *
 *  @brief Uses a set of partitioning coefficients to determine mass
 *  assimilation rates for several plant organs.
 *
 *  ### Partitioning overview
 *
 *  BioCro includes several partitioning growth calculators that determine these
 *  rates using slightly different methods. The different modules can be
 *  distinguished by the sets of tissues they use, the ways they apply
 *  respiration and water stress, and their responses to negative canopy
 *  assimilation rates. (A negative canopy assimilation rate indicates that the
 *  leaves are respiring more than they are photosynthesizing.)
 *
 *  In all partitioning growth calculators, the base growth rate for an organ is
 *  determined from the net canopy assimilation rate and a coefficient that
 *  determines the fraction of the net assimilation that is "partitioned" to
 *  that organ. Then, further modifications may take place to account for water
 *  stress, maintenance respiration, or other processes that affect the amount
 *  of carbon available to the organ for growth. Note that losses due to
 *  senescence and gains due to remobilized carbon from other organs are handled
 *  elsewhere and are not included here.
 *
 *  Respiration is included via the `resp()` function, which implements an
 *  empirical rule for determining the fraction of energy spent on respiration
 *  at a particular temperature. See the following paper for a general
 *  discussion of the importance of respiration in understanding plant growth:
 *  [Amthor, J. S. "The role of maintenance respiration in plant growth" Plant,
 *  Cell & Environment 7, 561–569 (1984)]
 *  (https://doi.org/10.1111/1365-3040.ep11591833).
 *
 *  ### Specifics of this module
 *
 *  In this module, no distinction is made between positive and negative canopy
 *  assimilation rates. Thus, respiratory losses in the leaf that result in a
 *  negative canopy assimilation rate are spread out to the other organs.
 *
 *  This module includes four organs:
 *  - `Leaf`: The leaf growth rate is modified by water stress and then
 *     respiration. Note that this effectively double-counts leaf respiration
 *     because the net canopy assimilation rate already includes it.
 *  - `Stem`: The stem growth rate is modified by respiration.
 *  - `Root`: The root growth rate is modified by respiration.
 *  - `Rhizome`: The rhizome growth rate is modified by respiration.
 *
 *  Here it is assumed that the major effect of water stress on mass
 *  accumulation is a reduction in the leaf growth rate, following
 *  [Boyer, J. S. "Leaf Enlargement and Metabolic Rates in Corn, Soybean, and
 *  Sunflower at Various Leaf Water Potentials" Plant Physiology 46, 233–235 (1970)]
 *  (https://doi.org/10.1104/pp.46.2.233).
 */
class partitioning_growth_calculator : public SteadyModule
{
   public:
    partitioning_growth_calculator(
        const state_map* input_parameters,
        state_map* output_parameters)
        :  // Define basic module properties by passing its name to its parent class
          SteadyModule("partitioning_growth_calculator"),

          // Get pointers to input parameters
          kLeaf_ip(get_ip(input_parameters, "kLeaf")),
          kStem_ip(get_ip(input_parameters, "kStem")),
          kRoot_ip(get_ip(input_parameters, "kRoot")),
          kRhizome_ip(get_ip(input_parameters, "kRhizome")),
          canopy_assimilation_rate_ip(get_ip(input_parameters, "canopy_assimilation_rate")),
          LeafWS_ip(get_ip(input_parameters, "LeafWS")),
          mrc1_ip(get_ip(input_parameters, "mrc1")),
          mrc2_ip(get_ip(input_parameters, "mrc2")),
          temp_ip(get_ip(input_parameters, "temp")),

          // Get pointers to output parameters
          newLeafcol_op(get_op(output_parameters, "newLeafcol")),
          newStemcol_op(get_op(output_parameters, "newStemcol")),
          newRootcol_op(get_op(output_parameters, "newRootcol")),
          newRhizomecol_op(get_op(output_parameters, "newRhizomecol"))
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();

   private:
    // Pointers to input parameters
    const double* kLeaf_ip;
    const double* kStem_ip;
    const double* kRoot_ip;
    const double* kRhizome_ip;
    const double* canopy_assimilation_rate_ip;
    const double* LeafWS_ip;
    const double* mrc1_ip;
    const double* mrc2_ip;
    const double* temp_ip;

    // Pointers to output parameters
    double* newLeafcol_op;
    double* newStemcol_op;
    double* newRootcol_op;
    double* newRhizomecol_op;

    // Main operation
    void do_operation() const;
};

string_vector partitioning_growth_calculator::get_inputs()
{
    return {
        "kLeaf",                     // dimensionless
        "kStem",                     // dimensionless
        "kRoot",                     // dimensionless
        "kRhizome",                  // dimensionless
        "canopy_assimilation_rate",  // Mg / ha / hour
        "LeafWS",                    // dimensionless
        "mrc1",                      // dimensionless
        "mrc2",                      // dimensionless
        "temp"                       // degrees C
    };
}

string_vector partitioning_growth_calculator::get_outputs()
{
    return {
        "newLeafcol",    // Mg / ha / hour
        "newStemcol",    // Mg / ha / hour
        "newRootcol",    // Mg / ha / hour
        "newRhizomecol"  // Mg / ha / hour
    };
}

void partitioning_growth_calculator::do_operation() const
{
    // Collect inputs and make calculations

    double kLeaf = *kLeaf_ip;
    double kStem = *kStem_ip;
    double kRoot = *kRoot_ip;
    double kRhizome = *kRhizome_ip;
    double canopy_assimilation_rate = *canopy_assimilation_rate_ip;
    double LeafWS = *LeafWS_ip;
    double mrc1 = *mrc1_ip;
    double mrc2 = *mrc2_ip;
    double temp = *temp_ip;

    double newLeafcol, newStemcol, newRootcol, newRhizomecol;

    // Calculate the rate of new leaf production
    if (kLeaf > 0) {
        newLeafcol = canopy_assimilation_rate * kLeaf * LeafWS;
        newLeafcol = resp(newLeafcol, mrc1, temp);
    } else {
        newLeafcol = 0.0;
    }

    // Calculate the rate of new stem production
    if (kStem >= 0) {
        newStemcol = canopy_assimilation_rate * kStem;
        newStemcol = resp(newStemcol, mrc1, temp);
    } else {
        newStemcol = 0.0;
    }

    // Calculate the rate of new root production
    if (kRoot > 0) {
        newRootcol = canopy_assimilation_rate * kRoot;
        newRootcol = resp(newRootcol, mrc2, temp);
    } else {
        newRootcol = 0.0;
    }

    // Calculate the rate of new rhizome production
    if (kRhizome > 0) {
        newRhizomecol = canopy_assimilation_rate * kRhizome;
        newRhizomecol = resp(newRhizomecol, mrc2, temp);
    } else {
        newRhizomecol = 0.0;
    }

    // Update the output parameter list
    update(newLeafcol_op, newLeafcol);
    update(newStemcol_op, newStemcol);
    update(newRootcol_op, newRootcol);
    update(newRhizomecol_op, newRhizomecol);
}

#endif
