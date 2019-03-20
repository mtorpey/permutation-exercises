OrbitStabiliserMT := function(G, a)
  # Input: a group (with generators) and a point
  # Output: a record
  local gens, orb, schreier, reps, stab_gens, pos, pt, g, new_pt, newgen, stab;
  gens := GeneratorsOfGroup(G);
  orb := [a];       # points in the orbit
  schreier := [0];  # refs to generators in Schreier vector
  reps := [()];     # coset reps of stabiliser (can be calculated from schreier)
  stab_gens := [];  # stabiliser by Schreier generators
  pos := 1;
  while pos <= Length(orb) do
    pt := orb[pos];
    for g in [1 .. Length(gens)] do
      new_pt := pt ^ gens[g];
      if not new_pt in orb then
        # New point in the orbit
        Add(orb, new_pt);
        Add(schreier, g);
        Add(reps, reps[pos] * gens[g]);
      else
        # New Schreier generator
        newgen := reps[pos] * gens[g] * reps[Position(orb, new_pt)]^-1;
        Add(stab_gens, newgen);
      fi;
    od;
    pos := pos + 1;
  od;
  #stab := Group(SmallGeneratingSet(Group(stab_gens)));
  stab := Group(stab_gens);
  return rec(gens := gens, orb := orb, schreier := schreier, stab := stab);
end;



StabiliserChain := function(G)
  local non_trivial_elm, pt, r;
  # Input: a perm group (with generators)
  # Output: a record
  non_trivial_elm := First(GeneratorsOfGroup(G), g-> not IsEmpty(MovedPoints(g)));
  if non_trivial_elm = fail then
    return rec(gens := []);
  fi;
  pt := MovedPoints(non_trivial_elm)[1];
  r := OrbitStabiliserMT(G, pt);
  r.stab := StabiliserChain(r.stab);
  return r;
end;



SchreierVectorCosetRep := function(v, a)
  local pos, elm, g;
  # Input: Schreier vector and int
  # Output: a perm (one that maps orb[1] to a)
  pos := Position(v.orb, a);
  if pos = fail then
    return fail;
  fi;
  elm := ();
  while v.schreier[pos] <> 0 do
    g := v.gens[v.schreier[pos]];
    elm := g * elm;
    pos := Position(v.orb, v.orb[pos] ^ (g^-1));
  od;
  return elm;
end;



FactoriseIntoReps := function(sc, g)
  local L, rep;
  # Input: Schreier vector (i.e. stab chain) and perm
  # Output: a list of perms (their product equals g)
  L := [];
  while not IsEmpty(sc.gens) do
    rep := SchreierVectorCosetRep(sc, sc.orb[1] ^ g);
    if rep = fail then
      return fail;
    fi;
    Add(L, rep);
    g := g * rep^-1;
    sc := sc.stab;
  od;
  if g <> () then
    return fail;
  fi;
  return Reversed(L);
end;



a_to_b := function(v, a, b)
  return SchreierVectorCosetRep(v, a)^-1 * SchreierVectorCosetRep(v, b);
end;
