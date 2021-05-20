#ifndef GRIMM_SOYBEAN_FLOWERING_H
#define GRIMM_SOYBEAN_FLOWERING_H

#include "../modules.h"
#include "../state_map.h"

/**
 * \brief Model for soybean development and flowering based on Grimm et al. (1993). See grimm_soybean_flowering_calculator for details.
 */
class grimm_soybean_flowering : public DerivModule
{
   public:
    grimm_soybean_flowering(
        const state_map* input_quantities,
        state_map* output_quantities) :  // Define basic module properties by passing its name to its parent class
                                                                      DerivModule("grimm_soybean_flowering"),
                                                                      // Get pointers to input quantities
                                                                      grimm_rate_ip(get_ip(input_quantities, "grimm_rate")),
                                                                      // Get pointers to output quantities
                                                                      grimm_physiological_age_op(get_op(output_quantities, "grimm_physiological_age"))
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();

   private:
    // Pointers to input quantities
    const double* grimm_rate_ip;
    // Pointers to output quantities
    double* grimm_physiological_age_op;
    // Main operation
    void do_operation() const;
};

string_vector grimm_soybean_flowering::get_inputs()
{
    return {
        "grimm_rate"};
}

string_vector grimm_soybean_flowering::get_outputs()
{
    return {
        "grimm_physiological_age"};
}

void grimm_soybean_flowering::do_operation() const
{
    // Collect inputs
    const double grimm_rate = *grimm_rate_ip;  // physiological days per hour

    // Update the output quantity list
    update(grimm_physiological_age_op, grimm_rate);
}

#endif
