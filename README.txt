This is Magma code to compute rational points on X_0^+(p) of genus 4-6 (where p is prime).

To understand this code, one needs to go through Section 3 of our paper
https://arxiv.org/abs/2105.04811

This code relies on Quadratic Chabauty package
https://github.com/steffenmueller/QCMod
developed by J. Balakrishnan, A. Best, F. Bianchi, N. Dogra, B. Lawrence, J. S. Müller, N. Triantafillou, J. Tuitman, and J. Vonk.

The intended way to use it is to download QCMod package and paste our files into the same main folder, so that calls such as
load "qc_modular.m"
can work.

Here is a short description of our files:
- genus_*_canonical_models.m contain canonical models for X0(N)+
- model_equation_finder has one main function, find_and_test_model, which, given tau_x and tau_y, computes the image of X0(N)+ in projective plane P^2, as well as bad and infinite F_p-points for primes p up to 50
- genus_*_find_and_test_model.m contains our maps and code to compute both the models and the rational points on them
- get_more_gonal_maps.m is an auxiliary file. Magma actually computes more than one gonal map, but only returns one. Since we needed more of these maps for genus 6 curves, this file was written for genus 6 case.
- preimage_checker.m checks that there are no additional rational points in the preimage of the image of the known rational points on our plane affine patch

- model_equation_finder_211.m is used only for the second patch of 211 (because it required an additional P^2-automorphism)
- model_equation_finder_new.m can be used when one wants to move the bad points to infinity
- the file 211_compute_model_1 computes the 1st model, 211_compute_model_2.m computes the 2nd model, 211_qc_modular.m performs QCModAffine on both models

- the file 359_example.m performs everything needed for level 359

Example usage:
load "qc_modular.m";
//this loads a number of other files from QCMod package

load "genus_4_canonical_models.m";
load "model_equation_finder.m";
X0p137_scheme, Q, X0p137_model_map, ratpts137, X0p137_good_primes, inf_pts, bad_pts := find_and_test_model(Z, Y, 42*Z, W + X + 2*Y + Z, X0_plus_137_galbraith, x, y, ratpts_137_galbraith);
X0p137_good_primes;
// [ 5 ]
// (you can print bad_pts to check that there are no bad points for p = 5, but then again you need to run QC to check if everything goes well)
// Important: there are no F_5-points on the canonical model of X_0^+(137) that map to infinite points on this model!
// Therefore, no rational points on X_0^+(137) can map to infinity on this model.


// For more information about quadratic Chabauty computations, you can raise the printlevel in the next command.
good_pts, bool, bad_pts, data, fake_pts, bad_disks := QCModAffine(Q, 5 : N:=20, printlevel:=1);
// The curve has only the known rational points outside the bad residue disks and
// disks at infinity

good_pts;
//[
//    [ -1/4, 1/8 ],
//    [ -1, 7 ],
//    [ -2/3, 4/27 ],
//    [ -1/2, 0 ],
//    [ -1/2, -1/4 ],
//    [ -1/2, -1/2 ],
//    [ 0, 2 ],
//    [ 0, 1 ],
//    [ 0, 0 ]
//]

bool;
// true
// Meaning that the computation proves that good_pts are all affine rational points in good residue disks.

bad_disks;
// []
// There are no bad disks.

// For more information about other outputs, one can check QCMod documentation (in qc_modular.m).

------------- PREIMAGE CHECK ---------------------

load "preimage_checker.m";
//load "genus_4_canonical_models.m"; if doing independently from above

additional_preimage(Z, Y, 42*Z, W + X + 2*Y + Z, X0_plus_137_galbraith, x, y, ratpts_137_galbraith);

//output is {@ @}
//there are no additional pre-images (additional rational points with the same image as the known rational points on our plane affine patch)






---------- COMMENT ON NEW QCMOD VERSION -----------
18-Jan-22: Recently, QCMod package was updated. One needs to pass additional precision argument for our examples to work, as it is done in this readme.file. Usually N:=20 suffices, but for some curves (levels 173, 197, 269 and 359), one can set precision to N:=30.
