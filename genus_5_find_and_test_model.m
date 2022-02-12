load "model_equation_finder.m";
load "genus_5_canonical_models.m";

/*
 * selectively (un)comment as needed
 */

/*
X0_plus_157_1_of_1, Patch_157_1_of_1_eq, model_map_157_1_of_1, ratpts_157_1_of_1_eq, good_primes_for_single_model_157_1_of_1, infinite_point_primes_data_model_157_1_of_1, bad_point_primes_data_model_157_1_of_1 := find_and_test_model(3*(V + 2*Y + 4*Z), V - 4*W - 3*X - 2*Y - 5*Z, V - W + Y, 2*V - 2*W - Y, X0_plus_157, x, y, ratpts_157);

*/


load "qc_modular.m";
p := 5;
Patch_157_1_of_1_eq := 
y^4 + (26/9*x^5 + 683/18*x^4 + 1772/9*x^3 + 1059/2*x^2 + 721*x + 765/2)*y^3 +
    (13/6*x^10 + 6461/108*x^9 + 80107/108*x^8 + 147358/27*x^7 + 710678/27*x^6 +
    260159/3*x^5 + 1185691/6*x^4 + 307704*x^3 + 312891*x^2 + 186921*x +
    49572)*y^2 + (-1/6*x^15 - 40/27*x^14 + 158267/1944*x^13 + 13058977/5832*x^12
    + 27072977/972*x^11 + 139928825/648*x^10 + 3383819957/2916*x^9 +
    327590123/72*x^8 + 360896509/27*x^7 + 178539895/6*x^6 + 300889759/6*x^5 +
    126224607/2*x^4 + 57559225*x^3 + 71894655/2*x^2 + 13743837*x + 4845663/2)*y
    - 7/18*x^20 - 2015/108*x^19 - 795385/1944*x^18 - 571046195/104976*x^17 -
    5137341277/104976*x^16 - 32912831723/104976*x^15 - 152493625145/104976*x^14
    - 31027531192/6561*x^13 - 482769464965/52488*x^12 + 261014179/216*x^11 +
    501455462005/5832*x^10 + 120320640259/324*x^9 + 27110528372/27*x^8 +
    71482892281/36*x^7 + 107782491011/36*x^6 + 3467242158*x^5 + 6094596371/2*x^4
    + 3943016613/2*x^3 + 1772756685/2*x^2 + 247339494*x + 32234193;
good_pts_157_1_of_1, bool_157_1_of_1, bad_pts_157_1_of_1, data_157_1_of_1, fake_rat_pts_157_1_of_1, bad_disks_157_1_of_1 := QCModAffine(Patch_157_1_of_1_eq, p : N:=20);
printf "good_pts_157_1_of_1: %o\n", good_pts_157_1_of_1;
printf "bool_157_1_of_1: %o\n", bool_157_1_of_1;
printf "bad_pts_157_1_of_1: %o\n", bad_pts_157_1_of_1;
printf "fake_rat_pts_157_1_of_1: %o\n", fake_rat_pts_157_1_of_1;
printf "bad_disks_157_1_of_1: %o\n\n", bad_disks_157_1_of_1;
/*
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_157_1_of_1: [
[ 0, -153/2 ],
[ 0, -153 ],
[ -3/2, 27/16 ],
[ -3/2, -27/32 ],
[ -3, 0 ],
[ -3, -27/2 ],
[ -6, -306 ],
[ 3, 0 ],
[ -1, 0 ]
]
bool_157_1_of_1: true
bad_pts_157_1_of_1: []
fake_rat_pts_157_1_of_1: [* *]
bad_disks_157_1_of_1: []

*/

/*
There isn't much freedom as to what we can pick for x if we want it to be a
3-to-1 map because there's a unique g_3^1 for this curve (the most we can do
is apply a linear fractional transformation to x). On the other hand, we can
basically pick y to be anything we want. However, we want to be careful if we
don't want the images of the rational points to stay distinct.
*/

/*
X0_plus_181_1_of_1, Patch_181_1_of_1_eq, model_map_181_1_of_1, ratpts_181_1_of_1_eq, good_primes_for_single_model_181_1_of_1, infinite_point_primes_data_model_181_1_of_1, bad_point_primes_data_model_181_1_of_1 := find_and_test_model(Z, (X + Y) + 8*Z, (2*17*127)*(V + W + X + 2*Y + 2*Z), 2*V - W + 3*X - 6*Y + 5*Z, X0_plus_181_galbraith, x, y, ratpts_181_galbraith);



load "qc_modular.m";
*/
p := 7;
Patch_181_1_of_1_eq := 
y^3 + (-11088*x^5 + 16124*x^4 - 9178*x^3 + 2548*x^2 - 344*x + 18)*y^2 +
    (37937948*x^10 - 109109424*x^9 + 139717096*x^8 - 104840620*x^7 +
    51026716*x^6 - 16825292*x^5 + 3805736*x^4 - 583156*x^3 + 57968*x^2 - 3380*x
    + 88)*y - 38931018912*x^15 + 165273816264*x^14 - 324830859688*x^13 +
    391968536592*x^12 - 324676043144*x^11 + 195506522152*x^10 - 88396329808*x^9
    + 30555659584*x^8 - 8141070384*x^7 + 1672000152*x^6 - 262592968*x^5 +
    30981584*x^4 - 2659496*x^3 + 156920*x^2 - 5696*x + 96;
good_pts_181_1_of_1, bool_181_1_of_1, bad_pts_181_1_of_1, data_181_1_of_1, fake_rat_pts_181_1_of_1, bad_disks_181_1_of_1 := QCModAffine(Patch_181_1_of_1_eq, p : N:=20);
printf "good_pts_181_1_of_1: %o\n", good_pts_181_1_of_1;
printf "bool_181_1_of_1: %o\n", bool_181_1_of_1;
printf "bad_pts_181_1_of_1: %o\n", bad_pts_181_1_of_1;
printf "fake_rat_pts_181_1_of_1: %o\n", fake_rat_pts_181_1_of_1;
printf "bad_disks_181_1_of_1: %o\n\n", bad_disks_181_1_of_1;
/*
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_181_1_of_1: [
[ 1/4, 3/512 ],
[ 1/4, 0 ],
[ 1/4, -3/512 ],
[ 1/3, 8/243 ],
[ 1/3, 4/243 ],
[ 1/3, 0 ],
[ 1, 648 ],
[ 0, -6 ],
[ 1/2, 9/8 ]
]
/*
bool_181_1_of_1: true
bad_pts_181_1_of_1: []
fake_rat_pts_181_1_of_1: [* *]
bad_disks_181_1_of_1: []

*/

/*
X0_plus_227_1_of_1, Patch_227_1_of_1_eq, model_map_227_1_of_1, ratpts_227_1_of_1_eq, good_primes_for_single_model_227_1_of_1, infinite_point_primes_data_model_227_1_of_1, bad_point_primes_data_model_227_1_of_1 := find_and_test_model(7*(W + X - Y + 5*Z), 2*V + 6*W + 5*X - Y + 17*Z, 7*17*Z, 2*V + 6*W + 5*X - Y + 17*Z, X0_plus_227_galbraith, x, y, ratpts_227_galbraith);


load "qc_modular.m";*/
p := 23;
Patch_227_1_of_1_eq := 
y^3 + (-16*x^3 + 49*x^2 - 41*x + 5)*y^2 + (102*x^6 - 668*x^5 + 1787*x^4 -
    2536*x^3 + 2094*x^2 - 1027*x + 252)*y - 289*x^9 + 3111*x^8 - 15134*x^7 +
    44025*x^6 - 85025*x^5 + 113565*x^4 - 104896*x^3 + 64299*x^2 - 23544*x + 3888;
good_pts_227_1_of_1, bool_227_1_of_1, bad_pts_227_1_of_1, data_227_1_of_1, fake_rat_pts_227_1_of_1, bad_disks_227_1_of_1 := QCModAffine(Patch_227_1_of_1_eq, p : N:=20);
printf "good_pts_227_1_of_1: %o\n", good_pts_227_1_of_1;
printf "bool_227_1_of_1: %o\n", bool_227_1_of_1;
printf "bad_pts_227_1_of_1: %o\n", bad_pts_227_1_of_1;
printf "fake_rat_pts_227_1_of_1: %o\n", fake_rat_pts_227_1_of_1;
printf "bad_disks_227_1_of_1: %o\n\n", bad_disks_227_1_of_1;
/*
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_227_1_of_1: [
[ 3/2, -3/8 ],
[ 3, 42 ],
[ 1, 0 ],
[ 2, 3 ]
]
bool_227_1_of_1: true
bad_pts_227_1_of_1: []
fake_rat_pts_227_1_of_1: [* *]
bad_disks_227_1_of_1: []

*/

/*
X0_plus_263_1_of_1, Patch_263_1_of_1_eq, model_map_263_1_of_1, ratpts_263_1_of_1_eq, good_primes_for_single_model_263_1_of_1, infinite_point_primes_data_model_263_1_of_1, bad_point_primes_data_model_263_1_of_1 := find_and_test_model(V - 2*X + Z, 2*V - W - X, W + X - Z, 2*V - W - X, X0_plus_263, x, y, ratpts_263);



load "qc_modular.m";*/
p := 23;
Patch_263_1_of_1_eq := 
y^6 + (-2*x + 9/2)*y^5 + (7/2*x^2 - 12*x + 12)*y^4 + (-3*x^3 + 19/2*x^2 - 14*x +
    14)*y^3 + (7/2*x^4 - 31/2*x^3 + 41/2*x^2 - 19/2*x + 13/2)*y^2 + (-x^5 + 9/2*x^4 -
    21/2*x^3 + 25/2*x^2 - 9/2*x + 1)*y + x^6 - 9/2*x^5 + 7*x^4 - 9/2*x^3 + 2*x^2 -
    1/2*x;
good_pts_263_1_of_1, bool_263_1_of_1, bad_pts_263_1_of_1, data_263_1_of_1, fake_rat_pts_263_1_of_1, bad_disks_263_1_of_1 := QCModAffine(Patch_263_1_of_1_eq, p : N:=20);
printf "good_pts_263_1_of_1: %o\n", good_pts_263_1_of_1;
printf "bool_263_1_of_1: %o\n", bool_263_1_of_1;
printf "bad_pts_263_1_of_1: %o\n", bad_pts_263_1_of_1;
printf "fake_rat_pts_263_1_of_1: %o\n", fake_rat_pts_263_1_of_1;
printf "bad_disks_263_1_of_1: %o\n\n", bad_disks_263_1_of_1;
/*
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_263_1_of_1: [
[ -17/110, -13/110 ],
[ 0, 0 ],
[ 0, -1/2 ],
[ 0, -1 ],
[ 1/2, 0 ],
[ 3/2, -1 ]
]
bool_263_1_of_1: true
bad_pts_263_1_of_1: []
fake_rat_pts_263_1_of_1: [* *]
bad_disks_263_1_of_1: []

*/
