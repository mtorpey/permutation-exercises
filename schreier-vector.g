SchreierVector := function(G, a)
  # Input: a group (with generators) and a point
  # Output: a record
  local gens, orb, schreier, pos, pt, g, new_pt;
  gens := GeneratorsOfGroup(G);
  orb := [a];
  schreier := [0];  # the identity
  pos := 1;
  while pos <= Length(orb) do
    pt := orb[pos];
    for g in [1 .. Length(gens)] do
      new_pt := pt ^ gens[g];
      if not new_pt in orb then
        Add(orb, new_pt);
        Add(schreier, g);
      fi;
    od;
    pos := pos + 1;
  od;
  return rec(group := G, orb := orb, schreier := schreier);
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
