////////////////////////////////////////////////////////////////
//      Gonal maps and plane models for genus 6 curves        //
//            mch - 12/11                                     //
//                                                            //
// Example usage:                                             //
//  load "genus_6_canonical_models.m";                        //
//  load "get_more_gonal_maps.m";                             //
//  maps := AllGenus6PlaneCurveModels(X0_plus_163);           //
//  maps;                                                     //
////////////////////////////////////////////////////////////////

function prime_ideal_to_point(P)
// P is a maximal ideal of polynomial ring R=k[x1,..xn]
// return a field L=k(alpha) which is the residue field R/P
// and the L-vector [x1 mod P,..,xn mod P]. Don't want to use
// FldAC as this gives the full splitting field

//NB: Assumed that Degree(P) le 5 and that if Degree(P)=4
// then R/P isn't a Galois extension of k [actually Deg(P)=4
// should never happen here anyway]

    R := Generic(P);
    k := BaseRing(R);
    n := Rank(R);

    GB := GroebnerBasis(P);
    assert #GB eq n;
    subs := [k!0 : i in [1..n]];
    j := n;
    while Degree(GB[j],j) eq 1 do
	cs := Coefficients(GB[j],j);
	assert cs[2] eq 1;
	subs[j] := -(k!cs[1]);
	j -:= 1;
	if j eq 0 then break; end if;
    end while;
    if j eq 0 then return subs; end if;
    pol := UnivariatePolynomial(GB[j]);
    L<a> := ext<k|pol>;
    ChangeUniverse(~subs,L);
    subs[j] := a;
    while j gt 1 do
	j -:= 1;
	f := GB[j];
	cs := Coefficients(f,j);
	assert #cs eq 2;
	assert cs[2] eq 1;
	subs[j] := -Evaluate(cs[1],subs);   
    end while;
    return subs;

end function;

function get_Y(X,res)
// return the uniquely determined surface Y that contains the 4-gonal canonical
//  genus 6 curve X

    R := CoordinateRing(Ambient(X));
    mp0 := Matrix(BoundaryMap(res,1));
    mp1 := Matrix(BoundaryMap(res,2));
    mp2 := Matrix(BoundaryMap(res,3));
    deg1 := [Max([LeadingTotalDegree(e) : e in r | e ne 0]) : 
			r in RowSequence(mp1)];
    deg2 := [Max([LeadingTotalDegree(e) : e in r | e ne 0]) : 
			r in RowSequence(Transpose(mp2))];
    assert Sort(deg1) eq [1: i in [1..5]] cat [2: i in [1..5]];
    assert &and[deg1[i]+deg2[i] eq 3 : i in [1..10]];

    perm := [i : i in [1..10] | deg1[i] eq 1] cat 
	[i : i in [1..10] |deg1[i] eq 2];
    iperm := [Index(perm,i) : i in [1..10]];

    mat1 := Matrix([mp1[p] : p in perm]);
    mat11 := RowSubmatrix(mat1,1,5);
    matf := Matrix([[MonomialCoefficient(mat11[k,j],R.i) : j in [1..Ncols(mat11)]]
           :i in [1..6], k in [1..Nrows(mat11)]]);
    b := Basis(Kernel(Transpose(matf)))[1];
    V := Kernel(Transpose(Matrix(b)));
    assert Dimension(V) eq 5;
    matV := EchelonForm(Matrix(Basis(V))); //EchelonForm usually unnecessary!
    B := Eltseq(ChangeRing(matV,R)*Matrix(mp0));
    return Scheme(Ambient(X),B),matV;

end function;

function get_mat_from_ap_data(dat,mat)
// return the 2x3 matrix giving the scroll coming from ap_data solution dat.
//  mat is the 5x5 matrix of linear forms giving R(-3)^5 -> R(-2)^5

    I := dat[2];
    c1,c2,d1,d2,d3 := Explode(dat[3]);
    k := BaseRing(I);
    if dat[1] eq 1 then
      P := Generic(I);
      rs := [k!NormalForm(P.i,I) : i in [7..9]];
      ss := [k!NormalForm(P.i,I) : i in [10..12]];
      L := k;
    else
      vec := prime_ideal_to_point(I);
      rs := vec[7..9];
      ss := vec[10..12];
      L := Universe(vec);
    end if;
    M := ZeroMatrix(L,2,5);
    M[1,c1] := 1; M[2,c2] := 1;
    inds := [i : i in [1..5] | i notin [c1,c2]];
    for i in [1..3] do
	M[1,inds[i]] := rs[i]; M[2,inds[i]] := ss[i];
    end for;

    R := BaseRing(mat);
    if not (L cmpeq k) then
	R := ChangeRing(R,L);
	mat := ChangeRing(mat,R);	
    end if;
    matsol := Submatrix(ChangeRing(M,R)*mat,[1,2],[d1,d2,d3]);
    return matsol,(L cmpeq k);

end function;

function g_1_4_ap_data(X,res,matY)
// X is a gonality 4 canonical genus 6 curve (with Y a Del Pezzo), res is the 1 6 10 6 1
//  minimal free resolution of R/IX. The rows of matY give the generators of the
//  5D subspace of R(-2)^6 that map onto the space of quadrics generating Y in res. 

    // Restrict the map R(-3)^5 -> R(-2)^6 to R(-3)^5 -> R(-2)^5
    R := CoordinateRing(Ambient(X));
    mp0 := Matrix(BoundaryMap(res,1));
    mp1 := Matrix(BoundaryMap(res,2)); // R(-2)^5+R(-3)^5 -> R(-2)^6
    mp2 := Matrix(BoundaryMap(res,3)); // R(-4)^6 -> R(-2)^5+R(-3)^5

    // permute the rows/cols at the rank 10 part of res to split into R(-3)^5+R(-4)^5
    //  mat1 will be the 6x5 matrix linear form matrix giving R(-3)^5 -> R(-2)^6
    deg1 := [Max([LeadingTotalDegree(e) : e in r | e ne 0]) : 
			r in RowSequence(mp1)];
    deg2 := [Max([LeadingTotalDegree(e) : e in r | e ne 0]) : 
			r in RowSequence(Transpose(mp2))];
    assert Sort(deg1) eq [1: i in [1..5]] cat [2: i in [1..5]];
    assert &and[deg1[i]+deg2[i] eq 3 : i in [1..10]];

    perm := [i : i in [1..10] | deg1[i] eq 1] cat 
	[i : i in [1..10] |deg1[i] eq 2];
    iperm := [Index(perm,i) : i in [1..10]];
    mat1 := Matrix([mp1[p] : p in perm]);
    mat1 := RowSubmatrix(mat1,1,5);

    // restrict to R(-3)^5 -> R(-2)^5
    inds := [Min(Support(matY[i])) : i in [1..5]];
    assert #Seqset(inds) eq 5;
    mat1r := Transpose(Matrix([r[i] : i in inds])) where r is RowSequence(Transpose(mat1));

    // Determine pairs (rank 2 direct summand M, rank 3 direct summand N) s.t. mat1r maps M -> N
    // M will be generated by vectors [1 0 r1 r2 r3], [0 1 s1 s2 s3] with some permutation of the 
    // columns.
    // N will be generated by vectors [1 0 0 u1 u2], [0 1 0 v1 v2], [0 0 1 w1 w2] with some
    // permtation of the columns

    P := PolynomialRing(BaseRing(X),12);
    AssignNames(~P,
	["u" cat IntegerToString(i) : i in [1..2]] cat
	["v" cat IntegerToString(i) : i in [1..2]] cat
	["w" cat IntegerToString(i) : i in [1..2]] cat
	["r" cat IntegerToString(i) : i in [1..3]] cat
	["s" cat IntegerToString(i) : i in [1..3]]);
    u1 := P.1; u2 := P.2;
    us := [P.i : i in [1,2]];
    v1 := P.3; v2 := P.4;
    vs := [P.i : i in [3,4]];
    w1 := P.5; w2 := P.6;
    ws := [P.i : i in [5,6]];
    r1 := P.7; r2 := P.8; r3 := P.9;
    rs := [P.i : i in [7..9]];
    s1 := P.10; s2 := P.11; s3 := P.12;
    ss := [P.i : i in [10..12]];
    PR := PolynomialRing(P,6);
    mat1R := ChangeRing(mat1r,PR);

    ap_data := [**];
    bdone := false;
    tot_deg := 0;

    for d1 in [1..3] do 
    for d2 in [d1+1..4] do 
    for d3 in [d2+1..5] do
    for c1 in [1..4] do 
    for c2 in [c1+1..5] do

      rels0 := [P|];
      // build the matrices for M and N with permutation of cols echelonising them and
      //  based on ci di, echelonise them and set various ri,si.. =0 in rels0
      matrs := ZeroMatrix(P,2,5);
      matuvw := ZeroMatrix(P,3,5);
      rels0 cat:= rs[1..c1-1]; matrs[1,c1] := 1; j := c1;
      for i in [k : k in [c1+1..5] | k ne c2] do matrs[1,i] := rs[j]; j +:= 1; end for;
      rels0 cat:= ss[1..c2-2]; matrs[2,c2] := 1; j := c2-1;
      for i in [c2+1..5] do matrs[2,i] := ss[j]; j +:= 1; end for;
      rels0 cat:= us[1..d1-1]; matuvw[1,d1] := 1; j := d1;
      for i in [k : k in [d1+1..5] | (k ne d2) and (k ne d3)] do matuvw[1,i] := us[j]; j +:= 1; end for;
      rels0 cat:= vs[1..d2-2]; matuvw[2,d2] := 1; j := d2-1;
      for i in [k : k in [d2+1..5] | k ne d3] do matuvw[2,i] := vs[j]; j +:= 1; end for;
      rels0 cat:= ws[1..d3-3]; matuvw[3,d3] := 1; j := d3-2;
      for i in [d3+1..5] do matuvw[3,i] := ws[j]; j +:= 1; end for;
//Append(~ap_data,<matrs,matuvw,rels0,[c1,c2,d1,d2,d3]>); continue;
      // matrix whose rows give the image of M in R(-2)^5
      Ls := ChangeRing(matrs,PR)*mat1R;

      // get equations for inclusion of M into N
      seq := [d1,d2,d3];
      rels := &cat[[r[i]-&+[matuvw[j,i]*r[seq[j]] : j in [1..3]]
	: i in [1..5] | i notin seq] : r in RowSequence(Ls)];
      rels1 := Setseq( Exclude( Seqset(&cat[[MonomialCoefficient(e,PR.i) : 
	i in [1..6]] : e in rels]), P!0) );

      I := ideal<P|rels0 cat rels1>;
      dim := Dimension(I);
      assert dim le 0;
      if dim eq 0 then
        deg := Degree(Homogenization(I));
	tot_deg +:= deg;
	if tot_deg ge 5 then bdone := true; end if;
	Append(~ap_data,<deg,I,[c1,c2,d1,d2,d3]>);
      end if;
      if bdone then break; end if;
    end for;
      if bdone then break; end if;
    end for;
      if bdone then break; end if;
    end for;
      if bdone then break; end if;
    end for;
      if bdone then break; end if;
    end for;

    apd := [dat : dat in ap_data | dat[1] eq 1];
    if true or (#apd eq 0) then
      apd := [];
      for dat in ap_data do
	if dat[1] eq 1 then
	    Append(~apd,dat);
	else
	    pcs := RadicalDecomposition(dat[2]);
	    for p in pcs do
		Append(~apd,<Degree(Homogenization(p)),p,dat[3]>);
	    end for;
	end if;
      end for;
    end if;
    ap_data := apd;
    return ap_data,mat1r;

end function;

function clean_up_equations(deqs)
  ambient_cr := Parent(deqs[1]);
  rank_ambient_cr := Rank(ambient_cr);
  deqs_mat := [[MonomialCoefficient(deq, ambient_cr.i) : i in [1..rank_ambient_cr]] : deq in deqs];
  deqs_clear_denom := LCM([Denominator(coef) : coef in Flat(deqs_mat)]);
  deqs_reduced := Rows(BasisMatrix(LLL(PureLattice(Lattice(ChangeRing(Matrix(deqs_mat) * deqs_clear_denom, Integers()))))));
  deqs_reduced_in_ring := [ &+[derow[i]*ambient_cr.i : i in [1..rank_ambient_cr]] : derow in deqs_reduced];
  deqs_reduced_sign := [ (LeadingCoefficient(deq) lt 0) select -deq else deq : deq in deqs_reduced_in_ring];
  return deqs_reduced_sign;
end function;

function AllGenus6GonalMaps(C)
  //gonality 4 (general) case
  g,bHyp,cmap := GenusAndCanonicalMap(C);
  Cc := Codomain(cmap);//canonical image
  res := MinimalFreeResolution(QuotientModule(Ideal(Cc)));
  Y,matY := get_Y(Cc, res);
  ap_data,mat := g_1_4_ap_data(Cc,res,matY);
  k := BaseRing(C);
 
  all_g6_gonal_maps := <>;
  for dat in ap_data do
    if dat[1] gt 1 then continue; end if;
    val := get_mat_from_ap_data(dat,mat);
    L := BaseRing(BaseRing(val));
    Append(~all_g6_gonal_maps, map<(k cmpeq L) select C else BaseChange(C, L)->Curve(ProjectiveSpace(L,1))|clean_up_equations([val[1,1],val[2,1]])>);
  end for;

  return all_g6_gonal_maps;
end function;

function AllGenus6PlaneCurveModels(C)
  //gonality 4 (general) case
  g,bHyp,cmap := GenusAndCanonicalMap(C);
  Cc := Codomain(cmap);//canonical image
  res := MinimalFreeResolution(QuotientModule(Ideal(Cc)));
  Y,matY := get_Y(Cc, res);
  ap_data,mat := g_1_4_ap_data(Cc,res,matY);
  k := BaseRing(C);
 
  all_g6_pcm := <>;
  for dat in ap_data do
    if dat[1] gt 1 then continue; end if;
    val := get_mat_from_ap_data(dat,mat);
    L := BaseRing(BaseRing(val));
   
    if not (k cmpeq L) then
         C1 := BaseChange(C,L);
    else
         C1 := C;
    end if;
    R1 := CoordinateRing(Ambient(C1));

    l1 := val[1,1]; l2 := val[2,1];
    I1 := Ideal(C1)+ideal<R1|l1>;
    ID := ColonIdeal(I1,l2);
    if Degree(ID) ne 4 then
       ID := ColonIdeal(I1,ideal<R1|l2>);
       ID := Saturation(ID);
    end if;
    assert Degree(ID) eq 4;
    lins := [l : l in Basis(ID) | TotalDegree(l) eq 1];
    assert #lins eq 3;

    Append(~all_g6_pcm, map<(k cmpeq L) select C else BaseChange(C, L)->(ProjectiveSpace(L,2))|lins>);
  end for;

  return all_g6_pcm;
end function;
