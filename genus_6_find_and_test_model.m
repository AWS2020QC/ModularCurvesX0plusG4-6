load "model_equation_finder.m";
load "genus_6_canonical_models.m";
load "get_more_gonal_maps.m";

/*
 * Gonal maps and minimal degree plane curve models defined over Q.
 *
 * 163 gonal map
 * 1. [X, Y + Z]
 * 2. [V, Y]
 * 163 plane curve model
 * 1. [U, V - Z, X]
 * 2. [U, V - Y, W + Z]
 *
 * 197 gonal map
 * 1. [X, Y]
 * 197 plane curve model
 * 1. [U + Z, V + W + Z, X - Y]
 *
 * 211 gonal map
 * 1. [X, V - W - 2*Y]
 * 2. [V + X - Z, W - X + Y]
 * 211 plane curve model
 * 1. [U + 2*Y + Z, V - W - 2*Y, X - Y - Z]
 * 2. [U - Z, V - Z, X]
 * 
 * 223 gonal map
 * 1. [V - Y, W - X - Z]
 * 223 plane curve model
 * 1. [U, V, W - X - Z]
 *
 * 269 gonal map
 * 1. [V + Z, W - Z]
 * 2. [Z, W - Y]
 * 3. [Y, Z]
 * 269 plane curve model
 * 1. [U + X + Y - 4*Z, V - X - 2*Y + 4*Z, W - X - 2*Y + 2*Z] 
 * 2. [V - Z, W - Y - Z, X - Z]
 * 3. [U - W + X, V - Z, Y - Z] 
 *
 * 271 gonal map: none
 * 271 plane curve model: none
 *
 * 359 gonal map
 * 1. [W + Y, V - Y + Z]
 * 359 plane curve model
 * 1. [U + V, W + Y, X]
 */

/*
 * selectively uncomment as needed
 */

/*
X0_plus_163_1_of_1, Patch_163_1_of_1_eq, model_map_163_1_of_1, ratpts_163_1_of_1_eq, good_primes_for_single_model_163_1_of_1, infinite_point_primes_data_model_163_1_of_1, bad_point_primes_data_model_163_1_of_1 := find_and_test_model(Y + Z, 2*X - Y - Z, V - Y, V + Y, X0_plus_163, x, y, ratpts_163 : printlevel := 3);

load "qc_modular.m";
p := 31;
Patch_163_1_of_1_eq := 
y^4 + (-1/40*x^4 - 2/5*x^3 - 11/20*x^2 - 1/5*x - 1/40)*y^3 + (-73/40*x^8 -
    9243/1600*x^7 - 13419/1600*x^6 - 11947/1600*x^5 - 7067/1600*x^4 -
    2873/1600*x^3 - 777/1600*x^2 - 129/1600*x - 9/1600)*y^2 + (33/40*x^12 +
    2643/800*x^11 + 414553/64000*x^10 + 10567/1280*x^9 + 484509/64000*x^8 +
    843/160*x^7 + 18241/6400*x^6 + 19383/16000*x^5 + 12861/32000*x^4 +
    409/4000*x^3 + 1213/64000*x^2 + 3/1280*x + 9/64000)*y - 3/8*x^16 -
    403/320*x^15 - 25921/12800*x^14 - 180533/102400*x^13 - 297327/640000*x^12 +
    1232487/1280000*x^11 + 212619/128000*x^10 + 3906469/2560000*x^9 +
    63421/64000*x^8 + 312841/640000*x^7 + 29941/160000*x^6 + 142093/2560000*x^5
    + 7981/640000*x^4 + 2583/1280000*x^3 + 27/128000*x^2 + 27/2560000*x;
good_pts_163_1_of_1, bool_163_1_of_1, bad_pts_163_1_of_1, data_163_1_of_1, fake_rat_pts_163_1_of_1, bad_disks_163_1_of_1 := QCModAffine(Patch_163_1_of_1_eq, p : printlevel:=3);
printf "good_pts_163_1_of_1: %o\n", good_pts_163_1_of_1;
printf "bool_163_1_of_1: %o\n", bool_163_1_of_1;
printf "bad_pts_163_1_of_1: %o\n", bad_pts_163_1_of_1;
printf "fake_rat_pts_163_1_of_1: %o\n", fake_rat_pts_163_1_of_1;
printf "bad_disks_163_1_of_1: %o\n\n", bad_disks_163_1_of_1;
+-------------------------------------------------------------------+
load "qc_modular.m";
p := 31;
Patch_163_1_of_1_eq := 
y^4 + (-1/40*x^4 - 2/5*x^3 - 11/20*x^2 - 1/5*x - 1/40)*y^3 + (-73/40*x^8 -
    9243/1600*x^7 - 13419/1600*x^6 - 11947/1600*x^5 - 7067/1600*x^4 -
    2873/1600*x^3 - 777/1600*x^2 - 129/1600*x - 9/1600)*y^2 + (33/40*x^12 +
    2643/800*x^11 + 414553/64000*x^10 + 10567/1280*x^9 + 484509/64000*x^8 +
    843/160*x^7 + 18241/6400*x^6 + 19383/16000*x^5 + 12861/32000*x^4 +
    409/4000*x^3 + 1213/64000*x^2 + 3/1280*x + 9/64000)*y - 3/8*x^16 -
    403/320*x^15 - 25921/12800*x^14 - 180533/102400*x^13 - 297327/640000*x^12 +
    1232487/1280000*x^11 + 212619/128000*x^10 + 3906469/2560000*x^9 +
    63421/64000*x^8 + 312841/640000*x^7 + 29941/160000*x^6 + 142093/2560000*x^5
    + 7981/640000*x^4 + 2583/1280000*x^3 + 27/128000*x^2 + 27/2560000*x;
good_pts_163_1_of_1, bool_163_1_of_1, bad_pts_163_1_of_1, data_163_1_of_1, fake_rat_pts_163_1_of_1, bad_disks_163_1_of_1 := QCModAffine(Patch_163_1_of_1_eq, p : printlevel:=3);
printf "good_pts_163_1_of_1: %o\n", good_pts_163_1_of_1;
printf "bool_163_1_of_1: %o\n", bool_163_1_of_1;
printf "bad_pts_163_1_of_1: %o\n", bad_pts_163_1_of_1;
printf "fake_rat_pts_163_1_of_1: %o\n", fake_rat_pts_163_1_of_1;
printf "bad_disks_163_1_of_1: %o\n\n", bad_disks_163_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_163_1_of_1: [
[ 1, 6/5 ],
[ 1, 0 ],
[ 0, 3/40 ],
[ 0, 1/40 ],
[ 0, 0 ],
[ 0, -3/40 ],
[ -1, 2/5 ],
[ -1, -2/5 ],
[ -1/5, 0 ],
[ -1/3, 0 ],
[ 7/6, 28981/3456 ]
]
bool_163_1_of_1: true
bad_pts_163_1_of_1: []
fake_rat_pts_163_1_of_1: [* *]
bad_disks_163_1_of_1: []


Total time: 2008.549 seconds, Total memory usage: 1737.50MB
 */
 
 /*
 X0_p_197_scheme_1_of_2, Q1, X0p197_model_map_1_of_2, ratpts197, X0p197_good_primes_1, inf_pts_1, bad_pts_1 := find_and_test_model(U+Z, 2*X-2*Y, V+W+Z, X-Y, X0_plus_197, x, y, ratpts_197);
 X0_p_197_scheme_2_of_2, Q2, X0p197_model_map_2_of_2, ratpts197, X0p197_good_primes_2, inf_pts_2, bad_pts_2 := find_and_test_model(U+Z, V + W + Z, X-Y,V+W+Z, X0_plus_197, x, y, ratpts_197);

load "qc_modular.m";
p := 23;
good_pts_197_1_of_2, bool_197_1_of_2, bad_pts_197_1_of_2, data_197_1_of_2, fake_rat_pts_197_1_of_2, bad_disks_197_1_of_2 := QCModAffine(Q1, p : printlevel:=3);
printf "good_pts_197_1_of_2: %o\n", good_pts_197_1_of_2;
printf "bool_197_1_of_2: %o\n", bool_197_1_of_2;
printf "bad_pts_197_1_of_2: %o\n", bad_pts_197_1_of_2;
printf "fake_rat_pts_197_1_of_2: %o\n", fake_rat_pts_197_1_of_2;
printf "bad_disks_197_1_of_2: %o\n\n", bad_disks_197_1_of_2;

good_pts_197_2_of_2, bool_197_2_of_2, bad_pts_197_2_of_2, data_197_2_of_2, fake_rat_pts_197_2_of_2, bad_disks_197_2_of_2 := QCModAffine(Q2, p : printlevel:=3);
printf "good_pts_197_2_of_2: %o\n", good_pts_197_2_of_2;
printf "bool_197_2_of_2: %o\n", bool_197_2_of_2;
printf "bad_pts_197_2_of_2: %o\n", bad_pts_197_2_of_2;
printf "fake_rat_pts_197_2_of_2: %o\n", fake_rat_pts_197_2_of_2;
printf "bad_disks_197_2_of_2: %o\n\n", bad_disks_197_2_of_2;

*/


/*
X0_plus_223_1_of_2, Patch_223_1_of_2_eq, model_map_223_1_of_2, ratpts_223_1_of_2_eq, good_primes_for_single_model_223_1_of_2, infinite_point_primes_data_model_223_1_of_2, bad_point_primes_data_model_223_1_of_2 := find_and_test_model(U - W, X + Y, 2*(V - Y), X + Y, X0_plus_223, x, y, ratpts_223 : printlevel := 3);

X0_plus_223_2_of_2, Patch_223_2_of_2_eq, model_map_223_2_of_2, ratpts_223_2_of_2_eq, good_primes_for_single_model_223_2_of_2, infinite_point_primes_data_model_223_2_of_2, bad_point_primes_data_model_223_2_of_2 := find_and_test_model(X + Y, U - W, 2*(V - Y), U - W, X0_plus_223, x, y, ratpts_223 : printlevel := 3);

load "qc_modular.m";
p := 19; 
Patch_223_1_of_2_eq := 
y^6 + (-x + 8)*y^5 + (8*x^2 + 30)*y^4 + (-4*x^3 + 88*x^2 + 40*x + 60)*y^3 +
    (8*x^4 - 32*x^3 + 264*x^2 + 128*x + 64)*y^2 + (-16*x^5 - 192*x^3 + 256*x^2 +
    176*x + 32)*y - 96*x^5 - 128*x^4 - 160*x^3 + 128*x;
Patch_223_2_of_2_eq :=
y^6 + (8*x - 1)*y^5 + (30*x^2 + 8)*y^4 + (60*x^3 + 40*x^2 + 88*x - 4)*y^3 +
    (64*x^4 + 128*x^3 + 264*x^2 - 32*x + 8)*y^2 + (32*x^5 + 176*x^4 + 256*x^3 -
    192*x^2 - 16)*y + 128*x^5 - 160*x^3 - 128*x^2 - 96*x;

good_pts_223_1_of_2, bool_223_1_of_2, bad_pts_223_1_of_2, data_223_1_of_2, fake_rat_pts_223_1_of_2, bad_disks_223_1_of_2 := QCModAffine(Patch_223_1_of_2_eq, p : printlevel:=3);
printf "good_pts_223_1_of_2: %o\n", good_pts_223_1_of_2;
printf "bool_223_1_of_2: %o\n", bool_223_1_of_2;
printf "bad_pts_223_1_of_2: %o\n", bad_pts_223_1_of_2;
printf "fake_rat_pts_223_1_of_2: %o\n", fake_rat_pts_223_1_of_2;
printf "bad_disks_223_1_of_2: %o\n\n", bad_disks_223_1_of_2;

good_pts_223_2_of_2, bool_223_2_of_2, bad_pts_223_2_of_2, data_223_2_of_2, fake_rat_pts_223_2_of_2, bad_disks_223_2_of_2 := QCModAffine(Patch_223_2_of_2_eq, p : printlevel:=3);
printf "good_pts_223_2_of_2: %o\n", good_pts_223_2_of_2;
printf "bool_223_2_of_2: %o\n", bool_223_2_of_2;
printf "bad_pts_223_2_of_2: %o\n", bad_pts_223_2_of_2;
printf "fake_rat_pts_223_2_of_2: %o\n", fake_rat_pts_223_2_of_2;
printf "bad_disks_223_2_of_2: %o\n\n", bad_disks_223_2_of_2;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_223_1_of_2: [
[ 0, 0 ],
[ 0, -2 ],
[ 1/2, -2 ],
[ -1, 0 ],
[ -1, -4 ],
[ 2/3, 0 ]
]
bool_223_1_of_2: true
bad_pts_223_1_of_2: []
fake_rat_pts_223_1_of_2: [* *]
bad_disks_223_1_of_2: []


Total time: 3858.010 seconds, Total memory usage: 2175.78MB
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_223_2_of_2: [
[ -1, 4 ],
[ -1, 0 ],
[ 3/2, 0 ],
[ 0, 0 ],
[ 2, -4 ]
]
bool_223_2_of_2: true
bad_pts_223_2_of_2: []
fake_rat_pts_223_2_of_2: [* *]
bad_disks_223_2_of_2: []


Total time: 2205.679 seconds, Total memory usage: 2143.75MB
*/


/*
load "qc_modular.m";
load "genus_6_canonical_models.m";
load "model_equation_finder.m";

X0_plus_269_1_of_1, Patch_269_1_of_1_eq, model_map_269_1_of_1, ratpts_269_1_of_1_eq, good_primes_for_single_model_269_1_of_1, infinite_point_primes_data_model_269_1_of_1, bad_point_primes_data_model_269_1_of_1 := find_and_test_model(2*Y -2*Z, 2*Y - Z, U, W- Y +2*Z, X0_plus_269, x, y, ratpts_269);
p:=29;

good_pts_269_1_of_1, bool_269_1_of_1, bad_pts_269_1_of_1, data_269_1_of_1, fake_rat_pts_269_1_of_1, bad_disks_269_1_of_1 := QCModAffine(Patch_269_1_of_1_eq, p : printlevel:=2);
printf "good_pts_269_1_of_1: %o\n", good_pts_269_1_of_1;
printf "bool_269_1_of_1: %o\n", bool_269_1_of_1;
printf "bad_pts_269_1_of_1: %o\n", bad_pts_269_1_of_1;
printf "fake_rat_pts_269_1_of_1: %o\n", fake_rat_pts_269_1_of_1;
printf "bad_disks_269_1_of_1: %o\n\n", bad_disks_269_1_of_1;
*/
/*
X0_plus_271_1_of_1, Patch_271_1_of_1_eq, model_map_271_1_of_1, ratpts_271_1_of_1_eq, good_primes_for_single_model_271_1_of_1, infinite_point_primes_data_model_271_1_of_1, bad_point_primes_data_model_271_1_of_1 := find_and_test_model(V + W, U - W + 2*X + Y, 80*(2*W - 2*X - Y + Z), U - W + 2*X + Y, X0_plus_271, x, y, ratpts_271 : printlevel := 3);

// Passing prec:=25 to the call to hodge_data in QCModAffine,
// prec:=20 was not enough.
load "qc_modular.m";
p := 13;
Patch_271_1_of_1_eq := 
y^5 + (-85*x^5 + 1127*x^4 - 1913*x^3 - 1182*x^2 - 163*x - 6)*y^4 + (10080*x^10 -
    140412*x^9 + 690654*x^8 - 1381518*x^7 + 483107*x^6 + 1493128*x^5 +
    629335*x^4 + 89925*x^3 + 672*x^2 - 820*x - 48)*y^3 + (51200*x^15 +
    5322240*x^14 - 62218208*x^13 + 256115488*x^12 - 497719884*x^11 +
    393876380*x^10 + 397941188*x^9 - 406993184*x^8 - 419868463*x^7 -
    107094101*x^6 + 10069328*x^5 + 10921601*x^4 + 2495279*x^3 + 280582*x^2 +
    16084*x + 376)*y^2 + (-77824000*x^20 + 842444800*x^19 - 1778347520*x^18 -
    10376567744*x^17 + 57185179088*x^16 - 105874923384*x^15 + 82483583448*x^14 +
    37040438788*x^13 - 128890371748*x^12 - 9768338538*x^11 + 103912556718*x^10 +
    65061175451*x^9 + 8351258644*x^8 - 6617560964*x^7 - 3771233527*x^6 -
    977933013*x^5 - 153099583*x^4 - 15265548*x^3 - 951824*x^2 - 33920*x - 528)*y
    + 1966080000*x^25 - 32800768000*x^24 + 195656908800*x^23 - 384548372480*x^22
    - 835604847872*x^21 + 5438540280384*x^20 - 10251070348080*x^19 +
    7222065612944*x^18 + 4271687117728*x^17 - 10456733694624*x^16 +
    2541993008328*x^15 + 7501013423496*x^14 - 4899818383384*x^13 -
    10352281169140*x^12 - 5911510862587*x^11 - 1394541574295*x^10 +
    91020107428*x^9 + 165638010722*x^8 + 56008509377*x^7 + 10935753977*x^6 +
    1399983446*x^5 + 119863400*x^4 + 6648960*x^3 + 216912*x^2 + 3168*x;
good_pts_271_1_of_1, bool_271_1_of_1, bad_pts_271_1_of_1, data_271_1_of_1, fake_rat_pts_271_1_of_1, bad_disks_271_1_of_1 := QCModAffine(Patch_271_1_of_1_eq, p : printlevel:=3);
printf "good_pts_271_1_of_1: %o\n", good_pts_271_1_of_1;
printf "bool_271_1_of_1: %o\n", bool_271_1_of_1;
printf "bad_pts_271_1_of_1: %o\n", bad_pts_271_1_of_1;
printf "fake_rat_pts_271_1_of_1: %o\n", fake_rat_pts_271_1_of_1;
printf "bad_disks_271_1_of_1: %o\n\n", bad_disks_271_1_of_1;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_271_1_of_1: [
[ -1, -414 ],
[ -3, 0 ],
[ 0, 6 ],
[ 0, 2 ],
[ 0, 0 ],
[ 1/3, 0 ]
]
bool_271_1_of_1: true
bad_pts_271_1_of_1: []
fake_rat_pts_271_1_of_1: [* *]
bad_disks_271_1_of_1: []


Total time: 123993.730 seconds, Total memory usage: 2656.25MB
*/


/*
X0_plus_359_1_of_2, Patch_359_1_of_2_eq, model_map_359_1_of_2, ratpts_359_1_of_2_eq, good_primes_for_single_model_359_1_of_2, infinite_point_primes_data_model_359_1_of_2, bad_point_primes_data_model_359_1_of_2 := find_and_test_model(W + Y, X, U + V, X, X0_plus_359, x, y, ratpts_359 : printlevel := 3);

X0_plus_359_2_of_2, Patch_359_2_of_2_eq, model_map_359_2_of_2, ratpts_359_2_of_2_eq, good_primes_for_single_model_359_2_of_2, infinite_point_primes_data_model_359_2_of_2, bad_point_primes_data_model_359_2_of_2 := find_and_test_model(X, U + V, W + Y, U + V, X0_plus_359, x, y, ratpts_359 : printlevel := 3);

load "qc_modular.m";
p := 7; 
Patch_359_1_of_2_eq := 
y^6 + (2*x + 1)*y^5 + (4*x^2 - 10*x - 6)*y^4 + (2*x^3 - 15*x^2 - 7*x + 9)*y^3 +
    (-4*x^4 - 10*x^3 - 9*x^2 + 12*x - 6)*y^2 + (-4*x^5 - 4*x^4 - 6*x^3 + 8*x^2 -
    3*x + 1)*y - x^6 - x^5 - 2*x^4 + x^3 - x^2;
Patch_359_2_of_2_eq :=
y^6 + (x + 4)*y^5 + (2*x^2 + 4*x + 4)*y^4 + (-x^3 + 6*x^2 + 10*x - 2)*y^3 + (x^4
    - 8*x^3 + 9*x^2 + 15*x - 4)*y^2 + (3*x^4 - 12*x^3 + 7*x^2 + 10*x - 2)*y -
    x^5 + 6*x^4 - 9*x^3 + 6*x^2 - x - 1;

good_pts_359_1_of_2, bool_359_1_of_2, bad_pts_359_1_of_2, data_359_1_of_2, fake_rat_pts_359_1_of_2, bad_disks_359_1_of_2 := QCModAffine(Patch_359_1_of_2_eq, p : printlevel:=3);
printf "good_pts_359_1_of_2: %o\n", good_pts_359_1_of_2;
printf "bool_359_1_of_2: %o\n", bool_359_1_of_2;
printf "bad_pts_359_1_of_2: %o\n", bad_pts_359_1_of_2;
printf "fake_rat_pts_359_1_of_2: %o\n", fake_rat_pts_359_1_of_2;
printf "bad_disks_359_1_of_2: %o\n\n", bad_disks_359_1_of_2;

good_pts_359_2_of_2, bool_359_2_of_2, bad_pts_359_2_of_2, data_359_2_of_2, fake_rat_pts_359_2_of_2, bad_disks_359_2_of_2 := QCModAffine(Patch_359_2_of_2_eq, p : printlevel:=3);
printf "good_pts_359_2_of_2: %o\n", good_pts_359_2_of_2;
printf "bool_359_2_of_2: %o\n", bool_359_2_of_2;
printf "bad_pts_359_2_of_2: %o\n", bad_pts_359_2_of_2;
printf "fake_rat_pts_359_2_of_2: %o\n", fake_rat_pts_359_2_of_2;
printf "bad_disks_359_2_of_2: %o\n\n", bad_disks_359_2_of_2;
+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_359_1_of_2: [
[ -2, 1 ],
[ 0, 1 ],
[ 0, 0 ],
[ -1, 1 ],
[ -1/3, 1/3 ]
]
bool_359_1_of_2: true
bad_pts_359_1_of_2: []
fake_rat_pts_359_1_of_2: [* *]
bad_disks_359_1_of_2: []

+-------------------------------------------------------------------+
Since their number matches the number of good affine rational points, this 
proves that this curve has only the known rational points outside the bad 
residue disks and disks at infinity
good_pts_359_2_of_2: [
[ 1, 0 ],
[ 1, -1 ],
[ 1, -2 ],
[ 0, 1 ],
[ 0, -1 ],
[ 3, -1 ]
]
bool_359_2_of_2: true
bad_pts_359_2_of_2: []
fake_rat_pts_359_2_of_2: [* *]
bad_disks_359_2_of_2: []

*/
