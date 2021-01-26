LiftedColouredSurface:=function(f)
    local 
        sub, clr, i, closure, j, k;

    sub:=RegularCWMapToCWSubcomplex(f)[2];
    clr:=List([1..3],x->List([1..Length(Target(f)!.boundaries[x])],y->[]));
    # list of colours of each cell of the tubular surface
    for i in sub[3] do # given in the indexing of the supercomplex
        Append(clr[3][i],f!.colour(2,i)*1);
        closure:=ClosureCWCell(Target(f),2,i)[2];
        for j in [1,2] do
            for k in [1..Length(closure[j])] do
                Append(clr[j][closure[j][k]],f!.colour(2,i)*1);
            od;
        od;
    od;
    for i in [1,2] do # clean-up the colours
        for j in [1..Length(clr[i])] do
            if IsBound(clr[i][j]) then
                clr[i][j]:=Set(clr[i][j]);
            fi;
        od;
    od;
    return clr;
end;
Read("~/a/2knots/arc-to-tubular-surface.g");
iota:=ArcDiagramToTubularSurface([[[2,4],[1,3],[2,4],[1,3]],[0,-1],[2]]);