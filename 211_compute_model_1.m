load "genus_6_find_and_test_model.m";
load "get_more_gonal_maps.m";
load "qc_modular.m";
load "model_equation_finder_new.m";

//1st model data
image_curve_1, image_curve_eq_1, model_map_1, image_ratpts_1,good_primes_for_single_model_1, \
infinite_point_primes_data_1,bad_point_primes_data_1, defining_eq_bad_disks_1, image_curve_non_monic_leading_coefficient_1, image_curve_non_monic_eq_1 \
:= find_and_test_model(X - Y - Z, V - W - 2*Y, U + 2*Y + Z, V - W - 2*Y, X0_plus_211, x, y, ratpts_211 : printlevel := 1, max_prime := 32);

infinite_point_primes_data_1;
//rec<recformat<p, infinite_points> | 
//p := 31,
//infinite_points := [ (0 : 0 : 0 : 0 : 0 : 1), (12 : 1 : 30 : 0 : 1 : 0), (18 : 1 : 30 
//    : 0 : 1 : 0), (0 : 0 : 0 : 0 : 0 : 1), (12 : 1 : 30 : 0 : 1 : 0), (18 : 1 : 30 : 0
//    : 1 : 0) ]>

bad_point_primes_data_1;
//no bad points for 31

printf "Q(x,y) = %o\np = 31\n", image_curve_eq_1;

preimage_schemes := {pt@model_map_1@@model_map_1 : pt in ratpts_211};
preimage_additional_rat_pts := &join{RationalPoints(sch) : sch in preimage_schemes} diff {X0_plus_211!pt : pt in ratpts_211};
printf "additional points: %o\n", preimage_additional_rat_pts;
// {@ @}



// 2nd model data
image_curve_2, image_curve_eq_2, model_map_2, image_ratpts_2,good_primes_for_single_model_2, \
infinite_point_primes_data_2,bad_point_primes_data_2, defining_eq_bad_disks_2, image_curve_non_monic_leading_coefficient_2, image_curve_non_monic_eq_2 \
:= find_and_test_model_new(U + 2*Y + Z, X - Y - Z, V - W - 2*Y, X - Y - Z, X0_plus_211, x, y, ratpts_211 : printlevel := 1, primes_to_test := [31]);

infinite_point_primes_data_2;
//rec<recformat<p, infinite_points> | 
//p := 31,
//infinite_points := [ (1 : 1 : 0 : 0 : 0 : 1), (0 : 0 : 0 : 0 : 1 : 0), 
//    (0 : 13 : 18 : 1 : 1 : 0), (0 : 19 : 12 : 1 : 1 : 0), (1 : 0 : 1 : 1
//    : 1 : 0) ]>
bad_point_primes_data_2;
// [], i.e. no bad points for 31

// clear denominators:
printf "Q(x,y) = %o\np = 31\n", image_curve_eq_2 * 16;

preimage_schemes := {pt@model_map_2@@model_map_2 : pt in ratpts_211};
preimage_additional_rat_pts := &join{RationalPoints(sch) : sch in preimage_schemes} diff {X0_plus_211!pt : pt in ratpts_211};
printf "additional points: %o\n", preimage_additional_rat_pts;
// {@ @}

/*image_curve_eq_2 := -x^4 + 4*x^3*y^2 - 6*x^3*y - 2*x^3 - 8*x^2*y^4 - 8*x^2*y^3 - 16*x^2*y^2 + 
10*x^2*y + 8*x^2 + 8*x*y^5 - 48*x*y^4 - 56*x*y^3 + 36*x*y^2 + 22*x*y - 6*x +
16*y^6 - 24*y^5 - 56*y^4 + 16*y^3 + 24*y^2 - 10*y + 1; // manually coerce to univariate polynomial ring
print image_curve_eq_2;
*/