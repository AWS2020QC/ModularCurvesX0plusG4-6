load "genus_6_find_and_test_model.m";
load "get_more_gonal_maps.m";
load "qc_modular.m";
load "model_equation_finder.m";


//1st model data
image_curve_1, image_curve_eq_1, model_map_1, image_ratpts_1,good_primes_for_single_model_1, \
infinite_point_primes_data_1,bad_point_primes_data_1, defining_eq_bad_disks_1, image_curve_non_monic_leading_coefficient_1, image_curve_non_monic_eq_1
:= find_and_test_model(W+Y, X, U+V, X, X0_plus_359, x, y, ratpts_359);

infinite_point_primes_data_1;
// referent part of the output:
//    rec<recformat<p, infinite_points> |
//        p := 7,
//        infinite_points := [ (0 : 1 : 0 : 0 : 6 : 1), (0 : 6 : 0 : 0 : 6 : 1),
//            (0 : 1 : 0 : 0 : 6 : 1), (0 : 6 : 0 : 0 : 6 : 1) ]>,
// each point is listed twice because denominators (X) are the same

bad_point_primes_data_1;
// there are no bad points for p := 7

good_pts_1, bool_1, bad_pts_1, data_1, fake_pts_1, bad_disks_1 := QCModAffine(image_curve_eq_1, 7 : printlevel:=2);
// The curve has only the known rational points outside the bad residue disks and
// disks at infinity

 good_pts_1;
//[
//    [ -2, 1 ],
//    [ 0, 1 ],
//    [ 0, 0 ],
//    [ -1, 1 ],
//    [ -1/3, 1/3 ]
//]
bool_1;
//true
fake_pts_1;
// [* *]
bad_disks_1;
// []

// This confirms that, going back to our canonical model, we need to show that the F_7-points(/disks of) (0 : 1 : 0 : 0 : 6 : 1), (0 : 6 : 0 : 0 : 6 : 1)
// do not contain rational points. We use 2nd model for that.


//2nd model data
image_curve_2, image_curve_eq_2, model_map_2, image_ratpts_2,good_primes_for_single_model_2, \
infinite_point_primes_data_2,bad_point_primes_data_2, defining_eq_bad_disks_2, image_curve_non_monic_leading_coefficient_2, image_curve_non_monic_eq_2
:= find_and_test_model(X, U+V, W+Y, U+V, X0_plus_359, x, y, ratpts_359);

infinite_point_primes_data_2;
// referent part:
// rec<recformat<p, infinite_points> |
//        p := 7,
//        infinite_points := [ (0 : 0 : 0 : 1 : 0 : 0), (0 : 0 : 0 : 1 : 0 : 0) ]>,
//

bad_point_primes_data_2;
// again, there are no bad points for p := 7 (the smallest prime for which there exist bad points is p := 17)

 good_pts_2, bool_2, bad_pts_2, data_2, fake_pts_2, bad_disks_2 := QCModAffine(image_curve_eq_2, 7 : printlevel:=2);
// The curve has only the known rational points outside the bad residue disks and
// disks at infinity

 good_pts_2;
//[
//    [ 1, 0 ],
//    [ 1, -1 ],
//    [ 1, -2 ],
//    [ 0, 1 ],
//    [ 0, -1 ],
//    [ 3, -1 ]
//]


bool_2;
//true
fake_pts_2;
// [* *]
bad_disks_2;
// []

// This confirms that the only rational points are (images of) our known rational points.


load "preimage_checker.m";
 additional_preimage(W+Y, X, U+V, X, X0_plus_359, x, y, ratpts_359);
//{@ @}
//there are no additional pre-images (additional rational points in the image of the known rational points on our 1st plane affine patch)

additional_preimage(X, U+V, W+Y, U+V, X0_plus_359, x, y, ratpts_359);
//{@ @}
//there are no additional pre-images (additional rational points in the image of the known rational points on our 2nd plane affine patch)
