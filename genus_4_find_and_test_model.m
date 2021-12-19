load "model_equation_finder.m";
load "genus_4_canonical_models.m";

/* P4 variables are W,X,Y,Z (in that order)
 * selectively uncomment as needed
 */


n := 1; // taking n := 0 givs bigger coefficients, other values of n might work as well
// our trigonal model will be the result of taking the image under 
// x = Z / Y
// y = 42*Z / (W + X + 2*Y + n*Z)

X0_p_137, patch_137, map_137, ratpts_137, good_primes_137, infinite_points_137, bad_points_137
 := find_and_test_model(Z, Y, 42*Z, W + X + 2*Y + n*Z, X0_plus_137_galbraith, x, y, ratpts_137_galbraith);


// Output of the QCModAffine code for p := 5
// for our purposes to find all rational points on X_{0}^{+}(137).

load "qc_modular.m";
p := 5;
good_pts_137_1_of_1, bool_137_1_of_1, bad_pts_137_1_of_1, data_137_1_of_1, fake_rat_pts_137_1_of_1, bad_disks_137_1_of_1 := QCModAffine(patch_137, p : printlevel:=2);
printf "good_pts_137_1_of_1: %o\n", good_pts_137_1_of_1;
printf "bool_137_1_of_1: %o\n", bool_137_1_of_1;
printf "bad_pts_137_1_of_1: %o\n", bad_pts_137_1_of_1;
printf "fake_rat_pts_137_1_of_1: %o\n", fake_rat_pts_137_1_of_1;
printf "bad_disks_137_1_of_1: %o\n\n", bad_disks_137_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_137_1_of_1: [
[ -1/4, 1/8 ],
[ -1, 7 ],
[ -2/3, 4/27 ],
[ -1/2, 0 ],
[ -1/2, -1/4 ],
[ -1/2, -1/2 ],
[ 0, 2 ],
[ 0, 1 ],
[ 0, 0 ]
]
bool_137_1_of_1: true
bad_pts_137_1_of_1: []
fake_rat_pts_137_1_of_1: [* *]
bad_disks_137_1_of_1: []

*/

/*
X0_plus_173_1_of_1, Patch_173_1_of_1_eq, model_map_173_1_of_1, ratpts_173_1_of_1_eq, good_primes_for_single_model_173_1_of_1, infinite_point_primes_data_model_173_1_of_1, bad_point_primes_data_model_173_1_of_1 := find_and_test_model(W + X + 2*Y + 4*Z, W + Y, (3*7/2)*(X - 3*Y - 3*Z), 2*W + X + 9*Y + 15*Z, X0_plus_173_galbraith, x, y, ratpts_173_galbraith);
*/

/*
load "qc_modular.m";
p := 5;
Patch_173_1_of_1_eq := 
    y^4 + (33/2*x^4 - 17/4*x^3 + 2/3*x^2 - 165/4*x + 9/2)*y^3 + (231/2*x^8 - 439/8*x^7 +
    97/24*x^6 - 166127/288*x^5 + 9093/32*x^4 - 385/6*x^3 + 22523/36*x^2 - 44117/288*x +
    581/96)*y^2 + (3087/8*x^12 - 3087/16*x^11 - 4109/16*x^10 - 244931/96*x^9 +
    807937/384*x^8 - 2087801/6912*x^7 + 1685255/288*x^6 - 9214061/2304*x^5 +
    1561699/1728*x^4 - 27504295/6912*x^3 + 1462645/864*x^2 - 319291/2304*x +
    1225/384)*y + 9261/16*x^16 - 12789/32*x^15 - 33831/32*x^14 - 423275/128*x^13 +
    1602851/384*x^12 + 805901/512*x^11 + 156997955/13824*x^10 - 2336128517/165888*x^9 +
    302880799/110592*x^8 - 1445252699/82944*x^7 + 695411285/41472*x^6 -
    29956235/6912*x^5 + 1455024709/165888*x^4 - 168700357/27648*x^3 + 399595/512*x^2 -
    224665/6144*x + 2401/4096;
good_pts_173_1_of_1, bool_173_1_of_1, bad_pts_173_1_of_1, data_173_1_of_1, fake_rat_pts_173_1_of_1, bad_disks_173_1_of_1 := QCModAffine(Patch_173_1_of_1_eq, p : printlevel:=3);
printf "good_pts_173_1_of_1: %o\n", good_pts_173_1_of_1;
printf "bool_173_1_of_1: %o\n", bool_173_1_of_1;
printf "bad_pts_173_1_of_1: %o\n", bad_pts_173_1_of_1;
printf "fake_rat_pts_173_1_of_1: %o\n", fake_rat_pts_173_1_of_1;
printf "bad_disks_173_1_of_1: %o\n\n", bad_disks_173_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_173_1_of_1: [
[ 1, 13/2 ],
[ 1, 0 ],
[ 0, -7/8 ],
[ 0, -21/8 ],
[ 1/7, -704/343 ],
[ 2/3, 497/72 ]
]
bool_173_1_of_1: true
bad_pts_173_1_of_1: []
fake_rat_pts_173_1_of_1: [* *]
bad_disks_173_1_of_1: []

*/


/*
X0_plus_199_1_of_1, Patch_199_1_of_1_eq, model_map_199_1_of_1, ratpts_199_1_of_1_eq, good_primes_for_single_model_199_1_of_1, infinite_point_primes_data_model_199_1_of_1, bad_point_primes_data_model_199_1_of_1 := find_and_test_model(W + Y, W - X + 2*Y - 3*Z, 5*11*(Y - Z), X - 4*Y + 7*Z, X0_plus_199_galbraith, x, y, ratpts_199_galbraith);
*/

/*
load "qc_modular.m";
p := 7;
Patch_199_1_of_1_eq :=
y^4 + (188*x^4 - 531*x^3 + 730*x^2 - 585*x + 195)*y^3 + (6930*x^8 - 36467*x^7 +
    102801*x^6 - 201149*x^5 + 280108*x^4 - 273201*x^3 + 179227*x^2 - 70731*x +
    12474)*y^2 + (87725*x^12 - 678755*x^11 + 2812416*x^10 - 8175060*x^9 +
    17971418*x^8 - 30862293*x^7 + 42039858*x^6 - 45235151*x^5 + 37608990*x^4 -
    23234260*x^3 + 9991212*x^2 - 2649537*x + 323433)*y + 332750*x^16 -
    3454550*x^15 + 18928470*x^14 - 72428151*x^13 + 212751798*x^12 -
    501847128*x^11 + 975020037*x^10 - 1580007615*x^9 + 2144054322*x^8 -
    2430428083*x^7 + 2279776105*x^6 - 1737806496*x^5 + 1046077245*x^4 -
    476258451*x^3 + 153252000*x^2 - 30873150*x + 2910897;
good_pts_199_1_of_1, bool_199_1_of_1, bad_pts_199_1_of_1, data_199_1_of_1, fake_rat_pts_199_1_of_1, bad_disks_199_1_of_1 := QCModAffine(Patch_199_1_of_1_eq, p : printlevel:=3);
printf "good_pts_199_1_of_1: %o\n", good_pts_199_1_of_1;
printf "bool_199_1_of_1: %o\n", bool_199_1_of_1;
printf "bad_pts_199_1_of_1: %o\n", bad_pts_199_1_of_1;
printf "fake_rat_pts_199_1_of_1: %o\n", fake_rat_pts_199_1_of_1;
printf "bad_disks_199_1_of_1: %o\n\n", bad_disks_199_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_199_1_of_1: [
[ -1, -246 ],
[ 1/2, -87/16 ],
[ 1/2, -87/8 ],
[ 0, -33 ],
[ 0, -99 ],
[ 1, 0 ],
[ 1, -1 ],
[ 3/4, -15/128 ]
]
bool_199_1_of_1: true
bad_pts_199_1_of_1: []
fake_rat_pts_199_1_of_1: [* *]
bad_disks_199_1_of_1: []

*/

/*
X0_plus_251_1_of_1, Patch_251_1_of_1_eq, model_map_251_1_of_1, ratpts_251_1_of_1_eq, good_primes_for_single_model_251_1_of_1, infinite_point_primes_data_model_251_1_of_1, bad_point_primes_data_model_251_1_of_1 := find_and_test_model(W + 7*X + 5*Z, W - 6*X + 10*Y - 23*Z, (5/3)*Z, X - Z, X0_plus_251_galbraith, x, y, ratpts_251_galbraith);
*/

/*
load "qc_modular.m"; // Warning: manually change the PointSearch bound in qc_modular.m to at least 10000
p := 11;
Patch_251_1_of_1_eq := 
y^4 + (47/12*x^4 - 70/3*x^3 - 797/15*x^2 - 701/30*x - 11/2)*y^3 + (167/36*x^8 -
    8243/120*x^7 + 88039/1800*x^6 + 2948981/3600*x^5 + 2329433/1800*x^4 +
    3074269/3600*x^3 + 147621/400*x^2 + 6785/72*x + 107/9)*y^2 + (85/36*x^12 -
    1543/27*x^11 + 611267/2160*x^10 + 10865777/10800*x^9 - 26623409/13500*x^8 -
    62633249/6750*x^7 - 522218437/43200*x^6 - 370500179/43200*x^5 -
    519277507/108000*x^4 - 19034281/8640*x^3 - 385939/540*x^2 - 3685/27*x -
    320/27)*y + 25/54*x^16 - 2545/162*x^15 + 2717/18*x^14 - 6926/135*x^13 -
    68913259/21600*x^12 - 552889403/108000*x^11 + 4952390317/1080000*x^10 +
    45374981579/3240000*x^9 + 271327319/80000*x^8 - 4566839179/432000*x^7 -
    105387676043/12960000*x^6 + 1344611/19200*x^5 + 476339231/172800*x^4 +
    11053139/6912*x^3 + 388631/864*x^2 + 7225/108*x + 350/81;
good_pts_251_1_of_1, bool_251_1_of_1, bad_pts_251_1_of_1, data_251_1_of_1, fake_rat_pts_251_1_of_1, bad_disks_251_1_of_1 := QCModAffine(Patch_251_1_of_1_eq, p : printlevel:=3);
printf "good_pts_251_1_of_1: %o\n", good_pts_251_1_of_1;
printf "bool_251_1_of_1: %o\n", bool_251_1_of_1;
printf "bad_pts_251_1_of_1: %o\n", bad_pts_251_1_of_1;
printf "fake_rat_pts_251_1_of_1: %o\n", fake_rat_pts_251_1_of_1;
printf "bad_disks_251_1_of_1: %o\n\n", bad_disks_251_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_251_1_of_1: [
[ 5, 9365/12 ],
[ -1/4, 1003/5120 ],
[ 0, 5/6 ],
[ -7/6, 0 ],
[ -4, -4847/10 ],
[ 1, 0 ]
]
bool_251_1_of_1: true
bad_pts_251_1_of_1: []
fake_rat_pts_251_1_of_1: [* *]
bad_disks_251_1_of_1: []

*/

/*
X0_plus_311_1_of_1, Patch_311_1_of_1_eq, model_map_311_1_of_1, ratpts_311_1_of_1_eq, good_primes_for_single_model_311_1_of_1, infinite_point_primes_data_model_311_1_of_1, bad_point_primes_data_model_311_1_of_1 := find_and_test_model(W + X + Z, W + Y - Z, 3*7*Z, W + X + Y, X0_plus_311_galbraith, x, y, ratpts_311_galbraith);
*/

/*
load "qc_modular.m";
p := 5;
Patch_311_1_of_1_eq := 
y^4 + (48*x^4 - 140*x^3 + 74*x^2 + 23*x - 1)*y^3 + (903*x^8 - 5223*x^7 +
    10128*x^6 - 6291*x^5 - 1761*x^4 + 2334*x^3 + 90*x^2 - 360*x + 72)*y^2 +
    (7938*x^12 - 68229*x^11 + 224046*x^10 - 327762*x^9 + 135747*x^8 + 159903*x^7
    - 147780*x^6 - 33084*x^5 + 58905*x^4 + 1044*x^3 - 11340*x^2 + 144*x +
    1044)*y + 27783*x^16 - 318843*x^15 + 1488564*x^14 - 3490182*x^13 +
    3760074*x^12 + 40986*x^11 - 3866940*x^10 + 2042064*x^9 + 1845207*x^8 -
    1526499*x^7 - 580176*x^6 + 553014*x^5 + 152064*x^4 - 109188*x^3 - 30240*x^2
    + 9288*x + 3024;
good_pts_311_1_of_1, bool_311_1_of_1, bad_pts_311_1_of_1, data_311_1_of_1, fake_rat_pts_311_1_of_1, bad_disks_311_1_of_1 := QCModAffine(Patch_311_1_of_1_eq, p : printlevel:=3);
printf "good_pts_311_1_of_1: %o\n", good_pts_311_1_of_1;
printf "bool_311_1_of_1: %o\n", bool_311_1_of_1;
printf "bad_pts_311_1_of_1: %o\n", bad_pts_311_1_of_1;
printf "fake_rat_pts_311_1_of_1: %o\n", fake_rat_pts_311_1_of_1;
printf "bad_disks_311_1_of_1: %o\n\n", bad_disks_311_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_311_1_of_1: [
[ 2, 0 ],
[ 2, -9 ],
[ 7/3, 0 ],
[ 1, 0 ],
[ 12/7, 132/49 ]
]
bool_311_1_of_1: true
bad_pts_311_1_of_1: []
fake_rat_pts_311_1_of_1: [* *]
bad_disks_311_1_of_1: []

*/
