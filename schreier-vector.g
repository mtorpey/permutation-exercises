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



SchreierMultiplier := function(v, a)
  local gens, pos, elm, g;
  # Input: Schreier vector and int
  # Output: a perm (one that maps the first elm in orbit to a)
  gens := GeneratorsOfGroup(v.group);
  if not a in v.orb then
    return fail;
  fi;
  pos := Position(v.orb, a);
  elm := ();
  while v.schreier[pos] <> 0 do
    g := gens[v.schreier[pos]];
    elm := g * elm;
    pos := Position(v.orb, v.orb[pos] ^ (g^-1));
  od;
  return elm;
end;



a_to_b := function(v, a, b)
  return SchreierMultiplier(v, a)^-1 * SchreierMultiplier(v, b);
end;
