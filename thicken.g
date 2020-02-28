Thicken:=function(iota)
    local
        M, Bn, dimM, dimBn, IsInterior;

    if not IsHapRegularCWMap(iota) then
        Error("input must be an inclusion of regular CW-complexes.");
    fi;

    M:=iota!.source*1;
    Bn:=iota!.target*1;

    dimM:=EvaluateProperty(M,"dimension");
    dimBn:=EvaluateProperty(Bn,"dimension");

    if not dimM + 2 = dimBn then
        Error("source complex must be a",dimBn-2"-complex.");
    fi;

    IsInterior:=function(map) # to ensure M < Int(Bn)
        local Sbnd, Tbnd, i, j, cell;

        Sbnd:=map!.source!.boundaries;
        Tbnd:=map!.target!.boundaries;

        for i in [1..Length(Sbnd)] do
            for j in [1..Length(Sbnd[i])] do
                cell:=[i,map!.mapping(i-1,j)];
                # check dem cobnds
            od;
        od;

        return true;
    end;

    if not IsInterior(iota) then
        Error(
   "the source complex is not a subspace of the interior of the target complex."
        );
    fi;
end;
