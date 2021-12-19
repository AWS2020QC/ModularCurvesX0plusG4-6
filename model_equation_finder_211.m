load "qc_modular.m";
P1<H,I> := ProjectiveSpace(Rationals(), 1);

function find_and_test_model_Q_only(x_num, x_denom, y_num, y_denom, domain_curve, x, y)
  map_to_P1P1 := map<domain_curve -> ProjectiveSpace(Rationals(),2) | [x_num*y_denom, x_denom*y_num, x_denom*y_denom]>;
  image_curve_non_monic_eq_xy := Evaluate(DefiningEquation(Image(map_to_P1P1)), [x, y, 1]);
  image_curve_non_monic_eq_xy := image_curve_non_monic_eq_xy / LeadingCoefficient(LeadingCoefficient(image_curve_non_monic_eq_xy));
  image_curve_non_monic_eq_xy_lc := LeadingCoefficient(image_curve_non_monic_eq_xy);
  image_curve_monic_eq := Numerator(Evaluate(image_curve_non_monic_eq_xy, y / image_curve_non_monic_eq_xy_lc));

  promising_primes := [];
  R := SquarefreePart(Discriminant(image_curve_monic_eq));
  if forall{f : f in Factorization(R) | (Degree(f[1]) gt 1) and (f[2] eq 1)} then
    print "this model doesn't have bad Q-points!";
    // there will be promising primes p by Chebotarev
    print "looking for promising primes p ...";
    promising_primes := [];
    for p in PrimesUpTo(100) do
      R_p := Factorization(ChangeRing(PolynomialRing(Rationals())!R, FiniteField(p)));
      if forall{f : f in Factorization(R) | (Degree(f[1]) gt 1) and (f[2] eq 1)} then
        printf "no bad F_p-points for p = %o\n", p;
        promising_primes cat:= [p];
      end if;
    end for;
  else
    print "this model has bad Q-points!";
  end if;

  return image_curve_monic_eq, image_curve_non_monic_eq_xy, promising_primes;
end function;

function TranslationToLineAtInfinity(projective_plane, points)
  if #points notin {2} then error("The number of points must be equal to 2."); end if;
  denom := 1;
  matrix_data := [Eltseq(points[2]), Eltseq(points[1])];
  for x in matrix_data[1] do
    denom := Lcm(Denominator(x), denom);
  end for;
  for x in matrix_data[2] do
    denom := Lcm(Denominator(x), denom);
  end for;
  matrix_data := [[Integers()!(denom * x) : x in elt] : elt in matrix_data];
  points_matrix := Matrix(Rationals(), #points,3, matrix_data);
  // M is a unimodular matrix such that M * P = (1 : 0 : 0) and M * Q = (0 : 1 : 0) lie on the line {z = 0} at infinity
  _, _, M := SmithForm(points_matrix);
  //print M;
  aut := Automorphism(projective_plane, M);
  //print [aut(p) : p in points];
  return aut;
end function;

        
        function point_reduction(pt, p, projective_plane)
          return projective_plane![GF(p)!pt[1], GF(p)!pt[2], GF(p)!pt[3]];
        end function;

                    function lift_point(pt, projective_plane)
              return projective_plane![Integers()!pt[1], Integers()!pt[2], Integers()!pt[3]];
            end function;

forward find_and_test_model_211_2;

function find_and_test_model_new_211(x_num, x_denom, y_num, y_denom, domain_curve, x, y, known_rat_pts : primes_to_test := [], max_prime := 50, printlevel := 1)
  if primes_to_test eq [] then
    primes_to_test := PrimesUpTo(max_prime);
  end if;
  projective_plane<A,B,C> := ProjectiveSpace(Rationals(), 2);
  map_x := Extend(map<domain_curve -> P1 | [x_num, x_denom]>);
  map_y := Extend(map<domain_curve -> P1 | [y_num, y_denom]>);
  if printlevel gt 1 then printf "The degree of the map to P^1_x is %o.\nThe degree of the map to P^1_y is %o.\n\n", Degree(map_x), Degree(map_y); end if;
  function_field_domain_curve := FunctionField(domain_curve);

  if printlevel gt 2 then
    printf "You probably want the images of rational points to stay away from infinity, i.e., avoid x(P) = [1 : 0], y(P) = [1 : 0].\n";
    for ratpt in known_rat_pts do
      printf "For P = %o, x(P) = %o and y(P) = %o.\n", ratpt, map_x(domain_curve ! ratpt), map_y(domain_curve! ratpt);
    end for;
  end if;

  // Get alternative equations because sometimes the regular equations won't be enough
  map_all_eqs := [];
  for def_x in AllDefiningPolynomials(map_x) do
    for def_y in AllDefiningPolynomials(map_y) do
      Append(~map_all_eqs, [ def_x[1], def_x[2], def_y[1], def_y[2] ]);
    end for;
  end for;

  // Get Magma to compute the model and force it to be monic
  map_components := [[def_xy[1] * def_xy[4], def_xy[3] * def_xy[2], def_xy[2] * def_xy[4]] : def_xy in map_all_eqs];

  map_to_P1P1 := map<domain_curve -> projective_plane | map_components>;
  map_to_P1P1 := Restriction(map_to_P1P1, Domain(map_to_P1P1), Image(map_to_P1P1));
  map_to_P1P1_degree := Degree(map_to_P1P1);

  if printlevel gt 1 then printf "Your map to P^2 has degree %o. You want this to be 1, or else it won't be an embedding.\n", map_to_P1P1_degree; end if;
  if map_to_P1P1_degree gt 1 then
    error("Your map to P^2 is not an embedding.");
  end if;

  image_curve_non_monic_eq := DefiningEquation(Codomain(map_to_P1P1));
  image_curve_non_monic_eq := image_curve_non_monic_eq / LeadingCoefficient(LeadingCoefficient(Evaluate(image_curve_non_monic_eq, [x, y, 1])));
  if printlevel gt 1 then printf "\nThe equation of the curve in P^1_x \\times P^1_y is \n%o.\n", Evaluate(image_curve_non_monic_eq, [x,y,1]); end if;
  image_curve_non_monic_leading_coefficient := LeadingCoefficient(Evaluate(image_curve_non_monic_eq, [x, y, 1]));

  image_curve_non_monic_degree := Degree(image_curve_non_monic_leading_coefficient);
  image_curve_non_monic_leading_term := Evaluate(image_curve_non_monic_leading_coefficient, A/C);
  image_curve := Scheme(projective_plane, [Numerator( Evaluate(image_curve_non_monic_eq, [A, B/image_curve_non_monic_leading_term, C]))]);
  map_leading_coefficients := Coefficients(image_curve_non_monic_leading_coefficient);

  // Write out equations explicitly
  if #map_leading_coefficients ge 2 then
    map_direct_components := [[
      def_xy[1] * def_xy[2]^(#map_leading_coefficients - 2) * def_xy[4],
      def_xy[3] * (&+[map_leading_coefficients[i] * (def_xy[1]^(i - 1)) * (def_xy[2]^(#map_leading_coefficients - i))  : i in [1 .. #map_leading_coefficients] ]),
      def_xy[2]^(#map_leading_coefficients - 1) * def_xy[4]
    ] : def_xy in map_all_eqs];
  else //#map_leading_coefficients eq 1 // already monic
    map_direct_components := [[
      def_xy[1] * def_xy[4],
      map_leading_coefficients[1] * def_xy[3] * def_xy[2],
      def_xy[2] * def_xy[4]
    ] : def_xy in map_all_eqs];
  end if;

  image_curve_eq := Evaluate(DefiningEquation(image_curve), [x, y, 1]);
  image_curve_eq := Parent(y) ! (image_curve_eq / LeadingCoefficient(image_curve_eq));
  if printlevel gt 1 then printf "\nThe equation of the affine model in A^2_{x,y} (set z = 1) is \n%o.\n", Evaluate(image_curve_non_monic_eq, [x,y,1]); end if;
  if not IsIrreducible(image_curve_eq) then
    error "This model for the curve is somehow not irreducible. (Changing the prime won't help.)";
  end if;

  model_map := map<domain_curve -> projective_plane | map_direct_components>;
  image_ratpts := [model_map(domain_curve ! ratpt) : ratpt in known_rat_pts];

  scheme_at_infinity_x := Scheme(domain_curve, [dp[2] : dp in AllDefiningPolynomials(map_x)]);
  scheme_at_infinity_y := Scheme(domain_curve, [dp[2] : dp in AllDefiningPolynomials(map_y)]);
  ///bound := 5000;
  for ratpt in known_rat_pts do
    if (model_map(domain_curve ! ratpt))[3] eq 0 and printlevel gt 1 then printf "The image of the Q-point %o lies at infinity for this model. (Changing the prime won't help.)\n", ratpt; end if;
  end for;

  defining_eq_bad_disks := SquarefreePart(Discriminant(image_curve_eq));
  roots_in_bad_disks := Roots(defining_eq_bad_disks);
  image_bad_Q_pts := {@ @};
  for i in [1..#known_rat_pts] do
    ratx := map_x(domain_curve ! known_rat_pts[i]);
    if ratx[2] eq 0 then
      continue;
    end if;
    if Evaluate(defining_eq_bad_disks, ratx[1] / ratx[2]) eq 0 then
      printf "The image of the Q-point %o is bad for this model. (Changing the prime won't help.)\n", known_rat_pts[i];
      print ratx;
      image_bad_Q_pt := model_map(domain_curve ! known_rat_pts[i]);
      Include(~image_bad_Q_pts, image_bad_Q_pt);
    end if;
  end for;
  printf "There are %o bad Q-points %o in the image for this model.\n", #image_bad_Q_pts, image_bad_Q_pts;


  if #image_bad_Q_pts eq 0 then
      print "There are no bad Q-point in the image. Try to move all pairs of points to infinity.";

      for p in primes_to_test do
        print "Trying p =", p;

        P2_Fp := Reduction(projective_plane, p);

        count := 0;
        //PQreds := [PQred : PQred in Subsets(Points(P2_Fp), 2) | Qred ne point_reduction(P, p, P2_Fp)];
        PQreds := Subsets(IndexedSetToSet(Points(P2_Fp)), 2);
        for PQred in PQreds do
          count +:= 1;
          PQred_Seq := SetToSequence(PQred);

          P := lift_point(PQred_Seq[1], projective_plane);
          Q := lift_point(PQred_Seq[2], projective_plane);

          printf "\n#%o/%o Choosing P = %o, Q = %o.\n", count, #PQreds, P, Q;
          aut := TranslationToLineAtInfinity(projective_plane, [P, Q]);
          printf "Apply the automorphism %o to the image curve.\n", aut;
          new_model_map := model_map * aut;
          //MinimizeReduce(Image(new_model_map));
          // repeat with new model
          //time image_curve2, image_curve_eq2, model_map2, image_ratpts2, good_primes_for_single_model2, \
          //infinite_point_primes_data2, bad_point_primes_data2, defining_eq_bad_disks2, image_curve_non_monic_leading_coefficient2, image_curve_non_monic_eq2 \
          //:= find_and_test_model_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
          return find_and_test_model_211_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
        end for;
      end for;
      print "\n\nDone.";
      return [];
  end if;
  

  if #image_bad_Q_pts gt 0 then
    if #image_bad_Q_pts eq 1 then // TODO: do this for all bad *F_p*-points
      print "There is exactly one bad Q-point in the image. Try to combine with all other (good) points.";
      P := image_bad_Q_pts[1];

     /* for good_point in image_ratpts do
        if good_point eq chosen_point then
          continue;
        end if;
        aut := TranslationToLineAtInfinity(projective_plane, [chosen_point, good_point]);
        printf "Apply the automorphism %o to the image curve.\n", aut;
        new_model_map := model_map * aut;
        // repeat with new model
          time image_curve2, image_curve_eq2, model_map2, image_ratpts2, good_primes_for_single_model2, \
        infinite_point_primes_data2, bad_point_primes_data2, defining_eq_bad_disks2, image_curve_non_monic_leading_coefficient2, image_curve_non_monic_eq2 \
        := find_and_test_model_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
      end for;

      print "\n\nDone.";
      return image_curve2, image_curve_eq2, model_map2, image_ratpts2, good_primes_for_single_model2, \
      infinite_point_primes_data2, bad_point_primes_data2, defining_eq_bad_disks2, image_curve_non_monic_leading_coefficient2, image_curve_non_monic_eq2;*/

      for p in primes_to_test do
        print "Trying p =", p;

        P2_Fp := Reduction(projective_plane, p);

        count := 0;
        Qreds := [Qred : Qred in Points(P2_Fp) | Qred ne point_reduction(P, p, P2_Fp)];
        //PQreds := Subsets(Points(P2_Fp), 2);
        for Qred in Qreds do
          count +:= 1;

          //P := lift_point(PQred[1], projective_plane);
          Q := lift_point(Qred, projective_plane);

          printf "\n#%o/%o Choosing P = %o, Q = %o.\n", count, #Qreds, P, Q;
          aut := TranslationToLineAtInfinity(projective_plane, [P, Q]);
          printf "Apply the automorphism %o to the image curve.\n", aut;
          new_model_map := model_map * aut;
          //MinimizeReduce(Image(new_model_map));
          // repeat with new model
          //time image_curve2, image_curve_eq2, model_map2, image_ratpts2, good_primes_for_single_model2, \
          //infinite_point_primes_data2, bad_point_primes_data2, defining_eq_bad_disks2, image_curve_non_monic_leading_coefficient2, image_curve_non_monic_eq2 \
          //:= find_and_test_model_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
          return find_and_test_model_211_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
        end for;
      end for;
      print "\n\nDone.";
      return [];
    else
      if #image_bad_Q_pts eq 2 then
        print "There are exactly 2 bad Q-points in the image. Try to move then to infinity.";
        chosen_points := [image_bad_Q_pts[i] : i in [1..Min(#image_bad_Q_pts, 2)]];
        aut := map < projective_plane -> projective_plane | [-C, B+C, A+2*C] >;
        printf "Apply the automorphism %o to the image curve.\n", aut;
        new_model_map := model_map * aut;
        // repeat with new model
        return find_and_test_model_211_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
      end if;
    end if;
  end if;

  bad_denominator_primes := LCM([Denominator(c) : c in Coefficients(defining_eq_bad_disks)]);

  infinite_point_data_format := recformat<p, infinite_points>;
  bad_point_data_format := recformat<p, bad_points>;

  if primes_to_test eq [] then
    primes_to_test := PrimesUpTo(max_prime);
  end if;
  bad_point_primes_data := [];
  infinite_point_primes_data := [];
  good_primes_for_single_model := [];

  for p in primes_to_test do
    Fp:=FiniteField(p);
    if printlevel gt 1 then print "\n"; end if;
    try
      if (bad_denominator_primes mod p eq 0) then
        error "bad prime because r(x) is not p-adically integral";
      end if;
      d:=Degree(image_curve_eq);

      A2:=AffineSpace(Fp,2);
      Fpxy:=CoordinateRing(A2);
      image_curve_modp:=Fpxy!0;
      C:=Coefficients(image_curve_eq);
      for i:=1 to #C do
        D:=Coefficients(C[i]);
        for j:=1 to #D do
            image_curve_modp:=image_curve_modp+(Fp!D[j])*Fpxy.1^(j-1)*Fpxy.2^(i-1);
        end for;
      end for;
      if not IsIrreducible(image_curve_modp) then
        error "bad prime because the mod p reduction of the model for the curve is somehow not irreducible";
      end if;

      g:=Genus(Curve(Scheme(A2,image_curve_modp)));
      r,Delta,s:=auxpolys(image_curve_eq);
      W0:=mat_W0(image_curve_eq);
      W0inv:=W0^(-1);
      Winf:=mat_Winf(image_curve_eq);
      Winfinv:=Winf^(-1);
      W:=Winf*W0^(-1);

      if (Fp!LeadingCoefficient(Delta) eq 0) or (Degree(r) lt 1) or (not smooth(r,p)) or (not (is_integral(W0,p) and is_integral(W0inv,p) and is_integral(Winf,p) and is_integral(Winfinv,p))) then
        if printlevel gt 1 then
	  print "At least one of the following is true, so we get a bad prime.";
	  printf "FpLeadingCoefficientDelta eq 0: %o\n", Fp!LeadingCoefficient(Delta) eq 0;
          Fpx:=PolynomialRing(Fp);
	  printf "Degree(r) < 1: %o\n", (Degree(r) lt 1);
	  printf "Degree(r) drops after reduction: %o\n", (not Degree(Fpx!r) eq Degree(r));
	  printf "r becomes inseparable after reduction: %o\n", (not IsSeparable(Fpx!r));
	  printf "Data not integral for W0, p: %o\n", (not is_integral(W0,p));
	  printf "Data not integral for W0inv, p: %o\n", (not is_integral(W0inv,p));
	  printf "Data not integral for Winf, p: %o\n", (not is_integral(Winf,p));
	  printf "Data not integral for Winfinv, p: %o\n", (not is_integral(Winfinv,p));
	end if;
        error "bad prime according to qc_modular";
      end if;
    catch e
      if printlevel gt 1 then printf "p = %o: %o.\n", p, e; end if;
      //p := NextPrime(p);
      continue;
    end try;
    domain_curve_mod_p := Reduction(domain_curve, p);
    Fp_points_at_infinity_x := [domain_curve_mod_p ! Eltseq(pt) : pt in Points(Reduction(scheme_at_infinity_x, p))];
    Fp_points_at_infinity_y := [domain_curve_mod_p ! Eltseq(pt) : pt in Points(Reduction(scheme_at_infinity_y, p))];
    Fp_points_at_infinity := Fp_points_at_infinity_x cat Fp_points_at_infinity_y;
    infinite_points_p_data := rec<infinite_point_data_format|>;
    infinite_points_p_data`p := p; infinite_points_p_data`infinite_points := Fp_points_at_infinity;
    if #Fp_points_at_infinity gt 0 then
      if printlevel gt 1 then printf "p = %o: there are %o points at infinity for the reduction modulo p.\n", p, #Fp_points_at_infinity; end if;
      Append(~infinite_point_primes_data, infinite_points_p_data);
    end if;
    Fp_bad_roots := Roots(ChangeRing(defining_eq_bad_disks, Fp));
    Fp_bad_points := [];
    for bad_x_mod_p in Fp_bad_roots do
      bad_points_upstairs := Points(Reduction(Scheme(domain_curve, [xm[1] - (Integers() ! bad_x_mod_p[1]) * xm[2] : xm in AllDefiningPolynomials(map_x)]), p));
      Fp_bad_points := Fp_bad_points cat [domain_curve_mod_p ! Eltseq(pt) : pt in bad_points_upstairs];
    end for;
    bad_points_p_data := rec<bad_point_data_format |>;
    bad_points_p_data`p := p; bad_points_p_data`bad_points := Fp_bad_points;
    if #Fp_bad_points gt 0 then
      if printlevel gt 1 then
        printf "p = %o: there are %o bad points for the reduction modulo p: %o.\n", p, #Fp_bad_points, Fp_bad_points;

        //lifts_of_bad_Fp_points := [pt : pt in image_ratpts | point_reduction(pt, p) in Fp_bad_points];
        //printf "reductions of %o\n.", lifts_of_bad_Fp_points;

       /* if #Fp_bad_points eq 1 then
          print "Moving this point via automorphisms to infinity.";
          P := Fp_bad_points[1];

          P2 := AmbientSpace(image_curve_mod_p);
          P2Q := AmbientSpace(image_curve);
          print P2;
          for Q in Points(P2) do
            if P eq Q then continue; end if;



            aut := TranslationToLineAtInfinity(P2Q, [lift_point(P), lift_point(Q)]);
            print "trying automorphism", aut;
          end for;
        end if;*/
      end if;
      Append(~bad_point_primes_data, bad_points_p_data);
    end if;
    if #Fp_points_at_infinity + #Fp_bad_points eq 0 then
      if printlevel gt 1 then printf "p = %o: no issues, possible candidate for a single model.\n", p; end if;
      Append(~good_primes_for_single_model, p);
    end if;
    //p := NextPrime(p);
  end for;

  return image_curve, image_curve_eq, model_map, image_ratpts, good_primes_for_single_model, infinite_point_primes_data, bad_point_primes_data, defining_eq_bad_disks, image_curve_non_monic_leading_coefficient, image_curve_non_monic_eq;
end function;

function find_and_test_model_211_2(model_map, known_rat_pts : primes_to_test := [], max_prime := 50, printlevel := 1)
  domain_curve := Domain(model_map);
  image_curve := Image(model_map);
  projective_plane<x,y,z> := AmbientSpace(image_curve);
  if printlevel gt 2 then
    printf "You probably want the images of rational points to stay away from infinity, i.e., avoid P = [x : y : 0].\n";
    for ratpt in known_rat_pts do
      printf "For P = %o, z(P) = %o.\n", ratpt, model_map(domain_curve ! ratpt)[3];
    end for;
  end if;

  //Qx<x> := PolynomialRing(Rationals(),1);
  //Qxy<y> := PolynomialRing(Qx, 1);
  //Qxyz<z> := PolynomialRing(Qxy, 1);
  //Qxyz<x,y,z> := PolynomialRing(Rationals(),3);
  image_curve_non_monic_eq := DefiningEquation(image_curve);
  image_curve_non_monic_leading_coefficient := LeadingCoefficient(Evaluate(image_curve_non_monic_eq, [x, y, 1]));
  if image_curve_non_monic_leading_coefficient ne 1 then
    image_curve_non_monic_eq := image_curve_non_monic_eq / LeadingCoefficient(image_curve_non_monic_leading_coefficient);
  end if;
  if printlevel gt 1 then printf "\nThe equation of the curve in P^2_{x,y,z} is \n%o.\n", Evaluate(image_curve_non_monic_eq, [x,y,1]); end if;
  image_curve_non_monic_leading_coefficient := LeadingCoefficient(Evaluate(image_curve_non_monic_eq, [x, y, 1]));

  QXY<X,Y> := PolynomialRing(Rationals(), 2);
  image_curve_eq := Evaluate(DefiningEquation(image_curve), [X, Y, 1]);
  leading_coefficient_in_Y := LeadingCoefficient(image_curve_eq, Y);
  if Degree(leading_coefficient_in_Y) gt 0 then
    print "This model is not monic in Y.";
    return [], [], [], [], [], [], [], [], [], [];
  end if;
  image_curve_eq := Parent(Y) ! (image_curve_eq / leading_coefficient_in_Y);
  if printlevel gt 1 then printf "\nThe equation of the affine model in A^2_{X,Y} (set z = 1) is \n%o.\n", image_curve_eq; end if;
  if not IsIrreducible(image_curve_eq) then
    error "This model for the curve is somehow not irreducible. (Changing the prime won't help.)";
  end if;

  image_ratpts := [model_map(domain_curve ! ratpt) : ratpt in known_rat_pts];

  scheme_at_infinity_z := Scheme(projective_plane, [z])@@model_map;

  for ratpt in known_rat_pts do
    if (model_map(domain_curve ! ratpt))[3] eq 0 and printlevel gt 1 then printf "The image of the Q-point %o lies at infinity for this model. (Changing the prime won't help.)\n", ratpt; end if;
  end for;

  QT<T> := PolynomialRing(Rationals());
  polynomial_in_T := hom<Parent(image_curve_eq) -> QT | T, 1>;
  defining_eq_bad_disks := polynomial_in_T(SquarefreePart(Discriminant(image_curve_eq, Y)));
  //print defining_eq_bad_disks;
  roots_in_bad_disks := Roots(defining_eq_bad_disks);
  number_of_bad_points_in_image := 0;
  image_bad_Q_pts := [];
  for i in [1..#known_rat_pts] do
    rat := model_map(domain_curve ! known_rat_pts[i]);
    if rat[2] eq 0 then // [x : 0 : 1]
      continue;
    end if;
    if Evaluate(defining_eq_bad_disks, rat[1] / rat[2]) eq 0 then
      printf "The image %o of the Q-point %o is bad for this model. (Changing the prime won't help.)\n", rat, known_rat_pts[i];
      number_of_bad_points_in_image +:= 1;
      Append(~image_bad_Q_pts, rat);
    end if;
  end for;
  if number_of_bad_points_in_image gt 0 then
    if number_of_bad_points_in_image eq 2 then
      printf "There are exactly 2 bad Q-points. Moving them to infinity.";
      chosen_points := [image_bad_Q_pts[i] : i in [1..Min(#image_bad_Q_pts, 2)]];
      aut := TranslationToLineAtInfinity(projective_plane, chosen_points);
      printf "Apply the automorphism %o to the image curve.\n", aut;
      new_model_map := model_map * aut;
      // repeat with new model
      return find_and_test_model_211_2(new_model_map, known_rat_pts : primes_to_test := primes_to_test, max_prime := max_prime, printlevel := printlevel);
    end if;
    printf "There are %o bad Q-points, this patch is not going to work.\n", number_of_bad_points_in_image;
    return [], [], [], [], [], [], [], [], [], [];
  else
    printf "There are no bad Q-points. Great!\n";
  end if;

  bad_denominator_primes := LCM([Denominator(c) : c in Coefficients(defining_eq_bad_disks)]);

  infinite_point_data_format := recformat<p, infinite_points>;
  bad_point_data_format := recformat<p, bad_points>;

  if primes_to_test eq [] then 
    primes_to_test := PrimesUpTo(max_prime);
  end if;
  bad_point_primes_data := [];
  infinite_point_primes_data := [];
  good_primes_for_single_model := [];

  for p in primes_to_test do
    Fp := FiniteField(p);
    if printlevel gt 1 then print "\n"; end if;
    try
      if (bad_denominator_primes mod p eq 0) then
        error "bad prime because r(x) is not p-adically integral";
      end if;
      d := Degree(image_curve_eq);

      A2 := AffineSpace(Fp,2);
      Fpxy := CoordinateRing(A2);
    /*  image_curve_modp := Fpxy!0;
      C := Coefficients(image_curve_eq);
      for i := 1 to #C do
        D := Coefficients(C[i]);
        for j := 1 to #D do
            image_curve_modp +:= (Fp!D[j]) * Fpxy.1^(j-1) * Fpxy.2^(i-1);
        end for;
      end for;*/
      image_curve_modp := Fpxy!image_curve_eq;

      if not IsIrreducible(image_curve_modp) then
        error "bad prime because the mod p reduction of the model for the curve is somehow not irreducible";
      end if;

      g := Genus(Curve(Scheme(A2, image_curve_modp)));
      /*r,Delta,s:=auxpolys(image_curve_eq);
      W0:=mat_W0(image_curve_eq);
      W0inv:=W0^(-1);
      Winf:=mat_Winf(image_curve_eq);
      Winfinv:=Winf^(-1);
      W:=Winf*W0^(-1);

      if (Fp!LeadingCoefficient(Delta) eq 0) or (Degree(r) lt 1) or (not smooth(r,p)) or (not (is_integral(W0,p) and is_integral(W0inv,p) and is_integral(Winf,p) and is_integral(Winfinv,p))) then
        if printlevel gt 1 then
	  print "At least one of the following is true, so we get a bad prime.";
	  printf "FpLeadingCoefficientDelta eq 0: %o\n", Fp!LeadingCoefficient(Delta) eq 0;
          Fpx:=PolynomialRing(Fp);
	  printf "Degree(r) < 1: %o\n", (Degree(r) lt 1);
	  printf "Degree(r) drops after reduction: %o\n", (not Degree(Fpx!r) eq Degree(r));
	  printf "r becomes inseparable after reduction: %o\n", (not IsSeparable(Fpx!r));
	  printf "Data not integral for W0, p: %o\n", (not is_integral(W0,p));
	  printf "Data not integral for W0inv, p: %o\n", (not is_integral(W0inv,p));
	  printf "Data not integral for Winf, p: %o\n", (not is_integral(Winf,p));
	  printf "Data not integral for Winfinv, p: %o\n", (not is_integral(Winfinv,p));
	end if;
        error "bad prime according to qc_modular";
      end if;*/
    catch e
      if printlevel gt 1 then printf "p = %o: %o.\n", p, e; end if;
      //p := NextPrime(p);
      continue;
    end try;
    domain_curve_mod_p := Reduction(domain_curve, p);
    // Fp-points of the curve reducing to a point on the line at infinity on the image
    Fp_points_at_infinity_z := [domain_curve_mod_p ! Eltseq(pt) : pt in Points(Reduction(scheme_at_infinity_z, p))];
    Fp_points_at_infinity := Fp_points_at_infinity_z;
    infinite_points_p_data := rec<infinite_point_data_format|>;
    infinite_points_p_data`p := p; infinite_points_p_data`infinite_points := Fp_points_at_infinity;
    if #Fp_points_at_infinity gt 0 then
      if printlevel gt 1 then printf "p = %o: there are %o points at infinity for the reduction modulo p.\n", p, #Fp_points_at_infinity; end if;
      Append(~infinite_point_primes_data, infinite_points_p_data);
    end if;
    if IsZero(ChangeRing(defining_eq_bad_disks, Fp)) then
      printf "p = %o: equation defining bad discs is zero. Skipping.\n", p;
      //p := NextPrime(p);
      continue;
    end if;
    Fp_bad_roots := Roots(ChangeRing(defining_eq_bad_disks, Fp));
    Fp_bad_points := [];
    for bad_x_mod_p in Fp_bad_roots do
      //pts_in_domain_with_bad_x := Scheme(domain_curve, [xm[1] - (Integers() ! bad_x_mod_p[1]) * xm[2] : xm in AllDefiningPolynomials(map_x)]);
      bad_x := bad_x_mod_p[1];
      //bad_y := ;
      pts_in_domain_with_bad_x := Scheme(image_curve, [x - (Integers()!bad_x)*z/*, y - (Integers()!bad_y)*z*/])@@model_map;
      bad_points_upstairs := Points(Reduction(pts_in_domain_with_bad_x, p));
      Fp_bad_points cat:= [domain_curve_mod_p ! Eltseq(pt) : pt in bad_points_upstairs];
    end for;
    bad_points_p_data := rec<bad_point_data_format |>;
    bad_points_p_data`p := p; bad_points_p_data`bad_points := Fp_bad_points;
    if #Fp_bad_points gt 0 then
      if printlevel gt 1 then printf "p = %o: there are %o bad points for the reduction modulo p.\n", p, #Fp_bad_points; end if;
      Append(~bad_point_primes_data, bad_points_p_data);
    end if;
    if /*#Fp_points_at_infinity +*/ #Fp_bad_points eq 0 then
      if printlevel gt 1 then printf "p = %o: no issues (no bad points for the reduction modulo p), possible candidate for a single model.\n", p; end if;
      Append(~good_primes_for_single_model, p);
    end if;
    //p := NextPrime(p);
  end for;

  return image_curve, image_curve_eq, model_map, image_ratpts, good_primes_for_single_model, infinite_point_primes_data, bad_point_primes_data, defining_eq_bad_disks, image_curve_non_monic_leading_coefficient, image_curve_non_monic_eq;
end function;
