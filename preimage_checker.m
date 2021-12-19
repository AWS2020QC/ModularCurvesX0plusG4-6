load "qc_modular.m";


//this function is used to check if there is any rational point lost in the steps in \varphi' -- when (x_num, x_denom), or (y_num,y_denom)=(0,0), or (x_denom,y_denom)=(0,0). I built this just to see what's going on, this step is also included in the additional_preimage function below, so we don't need to use this simple function.

function missing_in_projection(x_num, x_denom, y_num, y_denom, domain_curve)
map_to_P3 := map<domain_curve -> ProjectiveSpace(Rationals(),3) | [x_num*y_num,x_num*y_denom, x_denom*y_num, x_denom*y_denom]>;
missing_points:=RationalPoints([1,0,0,0]@@map_to_P3 , Rationals());
//when magma try to find the preimage of [1,0,0,0], it does not impose any condition on the first coordinate, so we are getting a possibly bigger set than the preimage of [1,0,0,0], but this is to our advantage -- the additional points are precisely the points lost in the map \tau_x\times \tau_y, so we covered all the missing points in the \varphi' step.
return missing_points;
end function;



//this function is what we need (except for 2nd model for 211), the input is the same as find_and_test_model, and the output is the set of additional preimages (which we expect to be empty):
function additional_preimage(x_num, x_denom, y_num, y_denom, domain_curve,x,y, known_rat_pts)

//Step 1: include all the points lost in the \varphi' map. More detailed explaination see the above function
map_to_P3 := map<domain_curve -> ProjectiveSpace(Rationals(),3) | [x_num*y_num,x_num*y_denom, x_denom*y_num, x_denom*y_denom]>;
missing_points:=RationalPoints([1,0,0,0]@@map_to_P3 , Rationals());

//Initial value for preimage_rat_pts set to be missing_points in the \varphi' step.
preimage_rat_pts:=missing_points;

//Step 2: Realize the composite map \varphi in codes -- modified from model_equation_finder.m
map_to_P1P1 := map<domain_curve -> ProjectiveSpace(Rationals(),2) | [x_num*y_denom, x_denom*y_num, x_denom*y_denom]>;//the map phi'
image_curve_non_monic_eq_xyz := DefiningEquation(Image(map_to_P1P1)); //the equation for C_N', before deviding through by leading coefficient
image_curve_non_monic_eq_xy := Evaluate(DefiningEquation(Image(map_to_P1P1)), [x, y, 1]);//dehomogenize the above
image_curve_non_monic_eq_xy := image_curve_non_monic_eq_xy / LeadingCoefficient(LeadingCoefficient(image_curve_non_monic_eq_xy));//make Q_0(x,1) monic
image_curve_non_monic_eq_xy_lc := LeadingCoefficient(image_curve_non_monic_eq_xy); //Q_0(x,1)
image_curve_monic_eq := Numerator(Evaluate(image_curve_non_monic_eq_xy, y / image_curve_non_monic_eq_xy_lc)); //Q(x,y)
d_1:=Degree(image_curve_non_monic_eq_xyz);
d_2:=Degree(image_curve_monic_eq);
d:=d_1-d_2; //total degree of Q_0
phi:= map<domain_curve -> ProjectiveSpace(Rationals(),2) | [x_num*x_denom^d*y_denom,Evaluate(image_curve_non_monic_eq_xy_lc,x_num/x_denom)*x_denom^d*x_denom*y_num, x_denom^(d+1)*y_denom]>;

//Step 3: Add preimage of the image of known_rat_pts under \varphi to the set preimage_rat_pts. I commented out some printing statements that I used to check my code.
    for ratpt in known_rat_pts do
        im:=phi(domain_curve!ratpt);
        //printf "The point %o is mapped to %o.\n", ratpt, im;
        //printf "All the points that map to %o is %o.\n", im, RationalPoints(im@@phi, Rationals());
        preimage_rat_pts:=preimage_rat_pts join RationalPoints(im@@phi, Rationals());
    end for;

//Step 4: From preimage_rat_pts, remove all the known_rat_pts.
additional_rat_pts:=preimage_rat_pts;
    for ratpt in known_rat_pts do
        additional_rat_pts:=additional_rat_pts diff {domain_curve!ratpt};
//print additional_rat_pts;
    end for;

return additional_rat_pts;
end function;




//This was an ad-hoc modification for the additional_preimage function for 2nd model of 211 where we post-composed an P^2 automorphism. Only need to change map_to_P1P1 and phi
//Currently it is not used, but we leave it for possible other applications.

function additional_preimage_211_2nd(x_num, x_denom, y_num, y_denom, domain_curve,x,y, known_rat_pts)

//Step 1: include all the points lost in the \varphi' map. Same as the general additional_preimage map above
map_to_P3 := map<domain_curve -> ProjectiveSpace(Rationals(),3) | [x_num*y_num,x_num*y_denom, x_denom*y_num, x_denom*y_denom]>;
missing_points:=RationalPoints([1,0,0,0]@@map_to_P3 , Rationals());

//Initial value for preimage_rat_pts set to be missing_points in the \varphi' step.
preimage_rat_pts:=missing_points;

//Step 2: Realize the composite map \varphi in codes -- modified from model_equation_finder.m
map_to_P1P1 := map<domain_curve -> ProjectiveSpace(Rationals(),2) | [-x_denom*y_denom,  x_denom*y_num + x_denom*y_denom, x_num*y_denom + 2*x_denom*y_denom]>;//the map phi' composed with our Aut(P^2) transformation -- the place we need to change
image_curve_non_monic_eq_xyz := DefiningEquation(Image(map_to_P1P1)); //the equation for C_N', before deviding through by leading coefficient
image_curve_non_monic_eq_xy := Evaluate(DefiningEquation(Image(map_to_P1P1)), [x, y, 1]);//dehomogenize the above
image_curve_non_monic_eq_xy := image_curve_non_monic_eq_xy / LeadingCoefficient(LeadingCoefficient(image_curve_non_monic_eq_xy));//make Q_0(x,1) monic
image_curve_non_monic_eq_xy_lc := LeadingCoefficient(image_curve_non_monic_eq_xy); //Q_0(x,1)
image_curve_monic_eq := Numerator(Evaluate(image_curve_non_monic_eq_xy, y / image_curve_non_monic_eq_xy_lc)); //Q(x,y)
d_1:=Degree(image_curve_non_monic_eq_xyz);
d_2:=Degree(image_curve_monic_eq);
d:=d_1-d_2; //total degree of Q_0
phi:= map<domain_curve -> ProjectiveSpace(Rationals(),2) | [(-x_denom*y_denom)*(x_num*y_denom + 2*x_denom*y_denom)^d,Evaluate(image_curve_non_monic_eq_xy_lc,(-x_denom*y_denom)/(x_num*y_denom + 2*x_denom*y_denom))*(x_num*y_denom + 2*x_denom*y_denom)^d*(x_denom*y_num + x_denom*y_denom), (x_num*y_denom + 2*x_denom*y_denom)^(d+1)]>; // another step we needed to change.

//Step 3: Add preimage of the image of known_rat_pts under \varphi to the set preimage_rat_pts. I commented out some printing statements that I used to check my code.
    for ratpt in known_rat_pts do
        im:=phi(domain_curve!ratpt);
        //printf "The point %o is mapped to %o.\n", ratpt, im;
        //printf "All the points that map to %o is %o.\n", im, RationalPoints(im@@phi, Rationals());
        preimage_rat_pts:=preimage_rat_pts join RationalPoints(im@@phi, Rationals());
    end for;

//Step 4: From preimage_rat_pts, remove all the known_rat_pts.
additional_rat_pts:=preimage_rat_pts;
    for ratpt in known_rat_pts do
        additional_rat_pts:=additional_rat_pts diff {domain_curve!ratpt};
//print additional_rat_pts;
    end for;

return additional_rat_pts;
end function;
