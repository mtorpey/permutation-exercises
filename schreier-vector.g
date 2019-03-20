StabiliserChain := function(G)
  local sc, gen;
  # Input: a perm group (with generators)
  # Output: a record
  sc := rec(gens := []);
  for gen in GeneratorsOfGroup(G) do
    ExtendSC(sc, gen);
  od;
  return sc;
end;



ExtendSC := function(sc, gen)
  local f, pos, nr_old_pts, pt, gens_to_apply, i, new_pt, newgen;
  # Input: a stabiliser chain and a perm
  # Output: a stabiliser chain (which now contains gen)
  f := FactoriseIntoReps(sc, gen);
  if f = fail then
    # gen is new

    if IsEmpty(sc.gens) then
      # Bottom of chain: create new layer
      sc.orb := [MovedPoints(gen)[1]];
      sc.schreier := [0];
      sc.stab := rec(gens := []);
      # TODO: special procedure for just one generator
    fi;
    
    Add(sc.gens, gen);
    
    # Apply the chain to the existing orbit
    pos := 1;
    nr_old_pts := Length(sc.orb);
    while pos <= Length(sc.orb) do
      pt := sc.orb[pos];
      if pos <= nr_old_pts then
        gens_to_apply := [Length(sc.gens)];  # old point: only the new gen
      else
        gens_to_apply := [1 .. Length(sc.gens)];  # new point: all the gens
      fi;
      for i in gens_to_apply do
        new_pt := pt ^ sc.gens[i];
        if not new_pt in sc.orb then
          # New point in the orbit
          Add(sc.orb, new_pt);
          Add(sc.schreier, i);
        else
          # New Schreier generator
          newgen := SchreierVectorCosetRep(sc, sc.orb[pos]) 
                    * sc.gens[i]
                    * SchreierVectorCosetRep(sc, new_pt) ^ -1;
          ExtendSC(sc.stab, newgen);
        fi;
      od;
      pos := pos + 1;
    od;
  fi;
  return sc;
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



PrintSC := function(sc)
  local i, size;
  i := 1;
  size := 1;
  while not IsEmpty(sc.gens) do
    Print("Level ", i, ": fix ", sc.orb[1], " (", Length(sc.orb), " cosets, ");
    Print(Length(sc.gens), " gens)\n");
    size := size * Length(sc.orb);
    sc := sc.stab;
    i := i + 1;
  od;
  Print("Group size: ", size, "\n");
end;
