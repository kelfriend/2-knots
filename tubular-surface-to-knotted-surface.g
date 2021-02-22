LiftColouredSurface:=function(f)
    local
        src, trg, lens, colours,
        cbnd4_1cells, cbnd4_1cells_bnd,
        i, j, clr, closure, k,
        min, max, l_colours,
        cell, col1, col2, int,
        bnd, bbnd, l, cobnd4,
        l_src, prods;

    trg:=RegularCWMapToCWSubcomplex(f);
    src:=trg[2]*1;
    trg:=trg[1];
    lens:=List([1..Length(trg!.boundaries)-1],x->Length(trg!.boundaries[x]));
    colours:=List([1..Length(src)],x->List([1..trg!.nrCells(x-1)],y->[]));
    cbnd4_1cells:=[]; # list of all 1-cells of src whose coboundary consists of 4 2-cells
    cbnd4_1cells_bnd:=[]; # the boundary of the above 1-cells
    for i in [1..Source(f)!.nrCells(1)] do
        if Source(f)!.coboundaries[2][i][1]=4 then
            Add(cbnd4_1cells,src[2][i]);
            for j in [2,3] do
                Add(
                    cbnd4_1cells_bnd,
                    src[1][Source(f)!.boundaries[2][i][j]]
                );
            od;
        fi;
    od;
    # we want to compute the inclusion src* -> trg x I where
    # I is the interval [min,max] and where min and max are the
    # smallest/largest integers that f!.colour assigns to 2-cells
    #       we begin by associating to each cell e^n of src a list
    # of integers e.g. [-1,2] indicating that e^n x {-1} and
    # e^n x {2} will be in src x I
    for i in [1..Length(src[Length(src)])] do
        clr:=f!.colour(Length(src)-1,src[Length(src)][i]);
        closure:=ClosureCWCell(trg,Length(src)-1,src[Length(src)][i])[2];
        for j in [1..Length(closure)] do
            for k in [1..Length(closure[j])] do
                colours[j][closure[j][k]]:=Set(
                    Concatenation(
                        colours[j][closure[j][k]],
                        clr
                    )
                );
            od;
        od;
    od;
    min:=Minimum(List(Filtered(colours[Length(colours)],x->x<>[]),Minimum));
    # ^smallest & \/largest 'colours'
    max:=Maximum(List(Filtered(colours[Length(colours)],x->x<>[]),Maximum));
    # make a copy of trg for each colour in [min-1,min+1]
    trg:=trg!.boundaries*1; Add(trg,[]);
    l_colours:=List(
        [1..Length(trg)],
        x->List([1..Length(trg[x])],
            y->[[x-1,y],min-1]
        )
    );
    for i in [min..max+1] do
        for j in [1..lens[1]] do
            Add(trg[1],[1,0]);
            Add(l_colours[1],[[0,j],i]);
        od;
    od;
    for i in [1..Length([min..max+1])] do
        for j in [1..Length(lens)-1] do
            for k in [1..lens[j+1]] do
                cell:=trg[j+1][k]*1;
                cell:=Concatenation(
                    [cell[1]],
                    cell{[2..Length(cell)]}+lens[j]*i
                );
                Add(trg[j+1],cell);
                Add(l_colours[j+1],[[j,k],[min..max+1][i]]);
            od;
        od;
    od;
    # `join' each Target(f) x {i-1} to Target(f) x {i}
    for i in [1..Length(lens)] do
        for j in [1..lens[i]] do
            cell:=[i-1,j];
            for k in [1..Length([min-1..max])] do
                col1:=[min-1..max][k];
                col2:=[min-1..max+1][k+1];
                int:=[col1,col2];
                bnd:=[]; # boundary of cell x int
                Add(
                    bnd,
                    Position(l_colours[i],[cell,col1])
                );
                Add(
                    bnd,
                    Position(l_colours[i],[cell,col2])
                );
                if i>1 then
                    bbnd:=trg[i][j]*1;
                    bbnd:=bbnd{[2..bbnd[1]+1]};
                    bbnd:=List(
                        bbnd,
                        x->[i-2,x]
                    );
                    for l in [1..Length(bbnd)] do
                        Add(
                            bnd,
                            Position(l_colours[i],[bbnd[l],int])
                        );
                    od;
                fi;
                bnd:=Set(bnd); Add(bnd,Length(bnd),1);
                Add(trg[i+1],bnd);
                Add(l_colours[i+1],[cell,int]);
            od;
        od;
    od;
    # identify the correct subcomplex of X x [min-1,max+1]
    cobnd4:=[]; # there are some 1-cells which shouldn't be lifted to 
    l_src:=List(src,x->[]);
    for i in [1..Length(colours)] do
        for j in [1..Length(colours[i])] do
            if colours[i][j]<>[] then
                prods:=[];
                if Length(colours[i][j])=1 then
                    Add(prods,colours[i][j][1]);
                else
                    int:=[Minimum(colours[i][j]),Maximum(colours[i][j])];
                    for k in [2..Length(int)] do
                        Add(prods,int[k-1]);
                        if not (i=2 and j in cbnd4_1cells) and
                            not (i=1 and j in cbnd4_1cells_bnd) then
                        Add(prods,[int[k-1],int[k]]);
                        fi;
                        Add(prods,int[k]);
                    od;
                    prods:=Set(prods);
                fi;
                for k in [1..Length(prods)] do
                    if IsInt(prods[k]) then
                        Add(
                            l_src[i],
                            Position(
                                l_colours[i],
                                [[i-1,j],prods[k]]
                            )
                        );
                    else
                        Add(
                            l_src[i+1],
                            Position(
                                l_colours[i+1],
                                [[i-1,j],prods[k]]
                            )
                        );
                    fi;
                od;
            fi;
        od;
    od;
    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(trg),
            l_src
        ]
    );
end;
