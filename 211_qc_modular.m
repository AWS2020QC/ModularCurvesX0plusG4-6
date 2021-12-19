load "qc_modular.m";

// patch 1

Q := y^6 + (9*x + 1)*y^5 + (32*x^2 + 12*x - 2)*y^4 + (57*x^3 + 49*x^2 - 6*x 
- 4)*y^3 + (53*x^4 + 89*x^3 + 5*x^2 - 19*x - 1)*y^2 + (24*x^5 + 72*x^4 + 30*x^3 
- 26*x^2 - 8*x + 2)*y + 4*x^6 + 20*x^5 + 24*x^4 - 8*x^3 - 14*x^2 + 3*x + 1;

p := 31;

good_pts_1, bool_1, bad_pts_1, data_1, fake_pts_1, bad_disks_1 := QCModAffine(Q, p : printlevel := 4);
// The curve has only the known rational points outside the bad residue disks and
// disks at infinity

good_pts_1;
/*[
    [ 1, -1 ],
    [ 1, -3 ],
    [ 0, -1 ],
    [ 1/2, -1 ],
    [ -1, 2 ],
    [ -1, 0 ],
    [ 2, -3 ]
]*/
bool_1;
//true
fake_pts_1;
// [* *]
bad_disks_1;
// []

// patch 2

Q := y^6 + (2 + 7*x)*y^5 + (-1 + 14*x + 18*x^2)*y^4 + (-4 - 9*x + 38*x^2 + 22*x^3)*y^3
+ (-2 - 22*x - 28*x^2 + 48*x^3 + 13*x^4)*y^2 + (1 - 6*x - 39*x^2 - 36*x^3 + 28*x^4 + 3*x^5)*y
+ 1 + 4*x - 2*x^2 - 22*x^3 - 17*x^4 + 6*x^5;

p := 31;

good_pts_2, bool_2, bad_pts_2, data_2, fake_pts_2, bad_disks_2 := QCModAffine(Q, p : printlevel := 4);
// The curve has only the known rational points outside the bad residue disks and
// disks at infinity

good_pts_2;
/* [
    [ 1, -2 ],
    [ 0, -1 ],
    [ -2, 3 ],
    [ -1, 2 ],
    [ -1, 1 ],
    [ -1/2, 0 ]
] */
bool_2;
//true
fake_pts_2;
// [* *]
bad_disks_2;
// []