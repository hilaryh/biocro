#ifndef TWO_LAYER_SOIL_PROFILE_H
#define TWO_LAYER_SOIL_PROFILE_H

#include "../modules.h"
#include "../state_map.h"
#include "AuxBioCro.h"  // For soilML_str
#include "BioCro.h"     // For soilML

class two_layer_soil_profile : public DerivModule {
	public:
		two_layer_soil_profile(state_map const& input_quantities, state_map* output_quantities) :
			// Define basic module properties by passing its name to its parent class
			DerivModule("two_layer_soil_profile"),
			// Get pointers to input quantities
			cws1_ip(get_ip(input_quantities, "cws1")),
			cws2_ip(get_ip(input_quantities, "cws2")),
			soil_depth1_ip(get_ip(input_quantities, "soil_depth1")),
			soil_depth2_ip(get_ip(input_quantities, "soil_depth2")),
			soil_depth3_ip(get_ip(input_quantities, "soil_depth3")),
			precip_ip(get_ip(input_quantities, "precip")),
			canopy_transpiration_rate_ip(get_ip(input_quantities, "canopy_transpiration_rate")),
			soil_field_capacity_ip(get_ip(input_quantities, "soil_field_capacity")),
			soil_wilting_point_ip(get_ip(input_quantities, "soil_wilting_point")),
			soil_saturation_capacity_ip(get_ip(input_quantities, "soil_saturation_capacity")),
			soil_air_entry_ip(get_ip(input_quantities, "soil_air_entry")),
			soil_saturated_conductivity_ip(get_ip(input_quantities, "soil_saturated_conductivity")),
			soil_b_coefficient_ip(get_ip(input_quantities, "soil_b_coefficient")),
			soil_sand_content_ip(get_ip(input_quantities, "soil_sand_content")),
			phi1_ip(get_ip(input_quantities, "phi1")),
			phi2_ip(get_ip(input_quantities, "phi2")),
			wsFun_ip(get_ip(input_quantities, "wsFun")),
			Root_ip(get_ip(input_quantities, "Root")),
			lai_ip(get_ip(input_quantities, "lai")),
			temp_ip(get_ip(input_quantities, "temp")),
			solar_ip(get_ip(input_quantities, "solar")),
			windspeed_ip(get_ip(input_quantities, "windspeed")),
			rh_ip(get_ip(input_quantities, "rh")),
			hydrDist_ip(get_ip(input_quantities, "hydrDist")),
			rfl_ip(get_ip(input_quantities, "rfl")),
			rsec_ip(get_ip(input_quantities, "rsec")),
			rsdf_ip(get_ip(input_quantities, "rsdf")),
			soil_clod_size_ip(get_ip(input_quantities, "soil_clod_size")),
			soil_reflectance_ip(get_ip(input_quantities, "soil_reflectance")),
			soil_transmission_ip(get_ip(input_quantities, "soil_transmission")),
			specific_heat_of_air_ip(get_ip(input_quantities, "specific_heat_of_air")),
			stefan_boltzman_ip(get_ip(input_quantities, "stefan_boltzman")),
			soil_water_content_ip(get_ip(input_quantities, "soil_water_content")),
			// Get pointers to output quantities
			cws1_op(get_op(output_quantities, "cws1")),
			cws2_op(get_op(output_quantities, "cws2")),
			soil_water_content_op(get_op(output_quantities, "soil_water_content"))
		{}
		static string_vector get_inputs();
		static string_vector get_outputs();
	private:
		// Pointers to input quantities
		const double* cws1_ip;
		const double* cws2_ip;
		const double* soil_depth1_ip;
		const double* soil_depth2_ip;
		const double* soil_depth3_ip;
		const double* precip_ip;
		const double* canopy_transpiration_rate_ip;
		const double* soil_field_capacity_ip;
		const double* soil_wilting_point_ip;
		const double* soil_saturation_capacity_ip;
		const double* soil_air_entry_ip;
		const double* soil_saturated_conductivity_ip;
		const double* soil_b_coefficient_ip;
		const double* soil_sand_content_ip;
		const double* phi1_ip;
		const double* phi2_ip;
		const double* wsFun_ip;
		const double* Root_ip;
		const double* lai_ip;
		const double* temp_ip;
		const double* solar_ip;
		const double* windspeed_ip;
		const double* rh_ip;
		const double* hydrDist_ip;
		const double* rfl_ip;
		const double* rsec_ip;
		const double* rsdf_ip;
		const double* soil_clod_size_ip;
		const double* soil_reflectance_ip;
		const double* soil_transmission_ip;
		const double* specific_heat_of_air_ip;
		const double* stefan_boltzman_ip;
		const double* soil_water_content_ip;
		// Pointers to output quantities
		double* cws1_op;
		double* cws2_op;
		double* soil_water_content_op;
		// Main operation
		void do_operation() const;
};

string_vector two_layer_soil_profile::get_inputs() {
	return {
		"cws1",
		"cws2",
		"soil_depth1",
		"soil_depth2",
		"soil_depth3",
		"precip",
		"canopy_transpiration_rate",
		"soil_field_capacity",
		"soil_wilting_point",
		"soil_saturation_capacity",
		"soil_air_entry",
		"soil_saturated_conductivity",
		"soil_b_coefficient",
		"soil_sand_content",
		"phi1",
		"phi2",
		"wsFun",
		"Root",
		"lai",
		"temp",
		"solar",
		"windspeed",
		"rh",
		"hydrDist",
		"rfl",
		"rsec",
		"rsdf",
		"soil_clod_size",
		"soil_reflectance",
		"soil_transmission",
		"specific_heat_of_air",
		"stefan_boltzman",
		"soil_water_content"
	};
}

string_vector two_layer_soil_profile::get_outputs() {
	return {
		"cws1",
		"cws2",
		"soil_water_content"
	};
}

void two_layer_soil_profile::do_operation() const {
	// Collect inputs and make calculations
	double cws1 = *cws1_ip;
	double cws2 = *cws2_ip;
	double soil_depth1 = *soil_depth1_ip;
	double soil_depth2 = *soil_depth2_ip;
	double soil_depth3 = *soil_depth3_ip;
	double soil_water_content = *soil_water_content_ip;

	double cws[] = {cws1, cws2};
	double soil_depths[] = {soil_depth1, soil_depth2, soil_depth3};

	struct soilML_str soilMLS = soilML(*precip_ip, *canopy_transpiration_rate_ip, cws, soil_depth3, soil_depths,
			*soil_field_capacity_ip, *soil_wilting_point_ip, *soil_saturation_capacity_ip, *soil_air_entry_ip, *soil_saturated_conductivity_ip,
			*soil_b_coefficient_ip, *soil_sand_content_ip, *phi1_ip, *phi2_ip, *wsFun_ip,
			2 /* Always uses 2 layers. */, *Root_ip, *lai_ip, 0.68, *temp_ip,
			*solar_ip, *windspeed_ip, *rh_ip, *hydrDist_ip, *rfl_ip,
			*rsec_ip, *rsdf_ip, *soil_clod_size_ip, *soil_reflectance_ip, *soil_transmission_ip,
			*specific_heat_of_air_ip, *stefan_boltzman_ip);

	double layer_one_depth = soil_depth2 - soil_depth1;
	double layer_two_depth = soil_depth3 - soil_depth2;
	double cws_mean = (soilMLS.cws[0] * layer_one_depth + soilMLS.cws[1] * layer_two_depth) / (layer_one_depth + layer_two_depth);

	// Update the output quantity list
	update(cws1_op, soilMLS.cws[0] - cws1);
	update(cws2_op, soilMLS.cws[1] - cws2);
	update(soil_water_content_op, cws_mean - soil_water_content);
}

#endif
