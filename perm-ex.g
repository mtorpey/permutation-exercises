# A permutation acts on the points [1..n] and nothing else
# Permutations of different degrees are incompatible

perm := function(list)
  # Input: a list (images of [1..n] under p)
  # Output: a record containing the list as img_list
  return rec(degree := Length(list), img_list := list);
end;

apply_perm := function(i, p)
  # Input: pos int and perm
  # Output: pos int (i^p) if i is not too high
  if IsPosInt(i) and i <= p.degree then
    return p.img_list[i];
  fi;
  return fail;
end;

inv_perm := function(p)
  # Input: a perm
  # Output: a perm (the inverse)
  local img_list, i;
  img_list := EmptyPlist(p.degree);
  for i in [1 .. p.degree] do
    img_list[apply_perm(i, p)] := i;
  od;
  return rec(degree := p.degree, img_list := img_list);
end;

mult_perms := function(p, q)
  # Input: two perms
  # Output: a perm (the product)
  # We compose maps from left to right
  local img_list;
  if p.degree <> q.degree then
    return fail;
  fi;
  img_list := List([1 .. p.degree], i -> apply_perm(apply_perm(i, p), q));
  return rec(degree := p.degree, img_list := img_list);
end;

random_perm := function(arg...)
  # Input: an int (the degree), or nothing
  # Output: a perm
  local deg, pts, img_list;
  if Length(arg) = 0 then
    deg := Random([1 .. 10]);
  elif Length(arg) = 1 then
    deg := arg[1];
  else
    return fail;
  fi;
  pts := [1 .. deg];
  img_list := List([1 .. deg], i -> Remove(pts, Random([1 .. deg - i + 1])));
  return rec(degree := deg, img_list := img_list);
end;

cycles_perm := function(p)
  # Input: a perm
  # Output: a list of lists (disjoint cycles of p)
  local cycles, pts;
  cycles := [];
  pts := BlistList([1 .. deg]);
  while not IsEmpty(pts) do
    i := Remove(pts, 1);
    # TODO
  od;
  return cycles;
end;

#display_perm := function(p)
  # Input: a perm
  # Output: a string (disjoint cycle notation)
