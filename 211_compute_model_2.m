load "genus_6_find_and_test_model.m";
load "get_more_gonal_maps.m";
load "qc_modular.m";
load "model_equation_finder_211.m";

// 2nd model data
image_curve_2, image_curve_eq_2, model_map_2, image_ratpts_2,good_primes_for_single_model_2, \
infinite_point_primes_data_2,bad_point_primes_data_2, defining_eq_bad_disks_2, image_curve_non_monic_leading_coefficient_2, image_curve_non_monic_eq_2 \
:= find_and_test_model_new_211(U + 2*Y + Z, X - Y - Z, V - W - 2*Y, X - Y - Z, X0_plus_211, x, y, ratpts_211 : printlevel := 1, primes_to_test := [31]);

infinite_point_primes_data_2;
/*[
    rec<recformat<p, infinite_points> | 
        p := 31,
        infinite_points := [ (1 : 1 : 0 : 0 : 0 : 1), (0 : 0 : 0 : 0 : 1 : 0), 
            (0 : 13 : 18 : 1 : 1 : 0), (0 : 19 : 12 : 1 : 1 : 0), (1 : 0 : 1 : 1
            : 1 : 0) ]>
]*/
bad_point_primes_data_2;
// [], i.e. no bad points for 31

// clear denominators:
printf "Q(x,y) = %o\np = 31\n", image_curve_eq_2;

Q := x^5*y + 2*x^5 + 13/3*x^4*y^2 + 28/3*x^4*y - 17/3*x^4 + 22/3*x^3*y^3 + 
16*x^3*y^2 - 12*x^3*y - 22/3*x^3 + 6*x^2*y^4 + 38/3*x^2*y^3 - 28/3*x^2*y^2 - 
13*x^2*y - 2/3*x^2 + 7/3*x*y^5 + 14/3*x*y^4 - 3*x*y^3 - 22/3*x*y^2 - 2*x*y + 
4/3*x + 1/3*y^6 + 2/3*y^5 - 1/3*y^4 - 4/3*y^3 - 2/3*y^2 + 1/3*y + 1/3;

Q *:= 3; // make Q monic in y

print Q;

preimage_schemes := {pt@model_map_2@@model_map_2 : pt in ratpts_211};
preimage_additional_rat_pts := &join{RationalPoints(sch) : sch in preimage_schemes} diff {X0_plus_211!pt : pt in ratpts_211};
printf "additional points: %o\n", preimage_additional_rat_pts;
// {@ @}

