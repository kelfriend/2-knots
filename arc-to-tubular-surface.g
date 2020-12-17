ArcDiagramToTubularSurface:=function(arc)
    local
        prs, crs, clr, grd, i, IsIntersection,
        CornerConfiguration, GRD, crossings, j,
        k, nr0cells, bnd, sub, hbars, hslice, cell,
        int, max, vbars, vslice, loops1, loops2,
        unchosen, neighbours, Clockwise,
        path, unselectedEdges, sourceORtarget,
        faceloops, x, ClockwiseTurn, 2cell, sORt,
        ori, e1, e0, loops, present_loops, vertices,
        check, l0, l1, l2, l0_, l1_, l2_, IsEdgeInDuplicate,
        copy1, hbars2, vbars2, copy2, 3cell, colour, lcap,
        ucap, floor, ceiling, cap, cap_, loop, colour_,
        leftovers, pos, HorizontalOrVertical, l, y,
        SubcapTo3cell, IntersectingCylinders, l0__, l1__,
        l2__, 0_ceiling, crs_int, h_0_ceiling, h_ceiling,
        h_cap, h_2_bnd, h_3_bnd, nr_crs, v_crossings,
        h_crossings, v_0_ceiling, v_ceiling, v_cap,
        minmax, v_2_bnd, v_3_bnd, perm, crs_2_cells,
        v_cap_sep, v_ceiling_sep, 3cells, int1, int2, m;

    if IsList(arc[1][1]) then
        prs:=arc[1]*1;
        crs:=arc[2]*1;
#                -1 |     +1 |      0 |
# crossing types: --|-- or ----- or --+--
#                   |        |        |
        if Length(arc)=3 then
            clr:=arc[3]*1;
# colours: 1; bgb, 2; bgr, 3; rgb, 4; rgr
        fi;
    else
        prs:=arc*1;
    fi;

# the 0-skeleton of the disk
####################################################################################
    grd:=List([1..Length(prs)],x->[1..Length(prs)]*0);
    for i in [0..Length(prs)-1] do
        grd[Length(prs)-i][prs[i+1][1]]:=1;
        grd[Length(prs)-i][prs[i+1][2]]:=1;
    od;

    IsIntersection:=function(i,j)
        if grd[i][j]=0 and
            1 in grd[i]{[1..j]} and
                1 in grd[i]{[j..Length(prs)]} and
                    1 in List([1..i],x->grd[x][j]) and
                        1 in List([i..Length(prs)],x->grd[x][j]) then
            return true;
        fi;
        return false;
    end;
    CornerConfiguration:=function(i,j);
        if grd[i][j]=1 then
            if Size(Positions(grd[i]{[j..Length(prs)]},1))=2 then
                if Size(Positions(List([i..Length(prs)],x->grd[x][j]),1))=2 then
                        # Corner type 1, i.e :  __
                    return 1; #                |
                elif Size(Positions(List([1..i],x->grd[x][j]),1))=2 then
                        # Corner type 3, i.e :
                    return 3; #                |__
                fi;
            elif Size(Positions(grd[i]{[1..j]},1))=2 then
                if Size(Positions(List([i..Length(prs)],x->grd[x][j]),1))=2 then
                        # Corner type 2, i.e :  __
                    return 2; #                   |
                elif Size(Positions(List([1..i],x->grd[x][j]),1))=2 then
                        # Corner type 4, i.e :
                    return 4; #                 __|
                fi;
            fi;
        fi;
        return 0;
    end;
    
    GRD:=List([1..4*Length(prs)],x->[1..4*Length(prs)]*0);
    crossings:=[];
# quadruple the size of grd to allow for the 0-skeleton
# to be displayed nicely without overlap in an array
    for i in [1..Length(prs)] do
        for j in [1..Length(prs)] do
            if CornerConfiguration(i,j) in [1,4] then
                GRD[4*i-3][4*j-3]:=1;
                GRD[4*i][4*j]:=1;
            elif CornerConfiguration(i,j) in [2,3] then
                GRD[4*i-3][4*j]:=1;
                GRD[4*i][4*j-3]:=1;
            elif IsIntersection(i,j) then
                for k in [0,3] do
                    GRD[4*i-3][4*j-3+k]:=1;
                    GRD[4*i][4*j-3+k]:=1;
                    Add(crossings,[4*i-3,4*j-3+k]);
                    Add(crossings,[4*i,4*j-3+k]);
                od;
            fi;
        od;
    od;
# label the 0-cells row by row
    nr0cells:=2;
    for i in [1..4*Length(prs)] do
        for j in [1..4*Length(prs)] do
            if GRD[i][j]=1 then
                GRD[i][j]:=nr0cells;
                nr0cells:=nr0cells+1;
            fi;
        od;
    od;
    crossings:=List(crossings,x->GRD[x[1]][x[2]]);
    crossings:=List([1..Length(crossings)/4],x->List([1..4]+4*x-4,y->crossings[y]));
    GRD:=FrameArray(GRD);
    GRD[1][1]:=1;
    GRD[4*Length(prs)+2][4*Length(prs)+2]:=nr0cells;

    bnd:=List([1..5],x->[]); # eventual boundary list of the 3-ball containing
    sub:=List([[],[],[]]); # the boundary of a knotted surface as a subcomplex
    bnd[1]:=List([1..nr0cells],x->[1,0]);
    sub[1]:=[2..nr0cells-1];
    if IsBound(crs) then
        for i in [1..Length(crs)] do
            if crs[i]=0 then
                sub[1]:=Difference(sub[1],crossings[i]);
            fi;
        od;
    fi;
####################################################################################

# the 1-skeleton of the disk
####################################################################################
    # add the horizontal 1-cells
    hbars:=[];
    for i in [2..Length(GRD)-1] do
        hslice:=Filtered(GRD[i],x->x<>0);
        if hslice<>hslice*0 then
            Add(hbars,hslice);
        fi;
        for j in [1..Length(hslice)-1] do
            cell:=[2,hslice[j],hslice[j+1]];
            Add(bnd[2],cell);
            int:=List(crossings,x->Length(Intersection(cell{[2,3]},x)));
            max:=PositionMaximum(int);
            if IsBound(crs) then
                if int[max]=1 and
                    crs[max]=-1 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2]));
                elif int[max]=2 and
                    crs[max]<>0 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2])); 
                fi;
            else
                if max<>fail then
                    if int[max]=2 then #  horizontal 1-cells default to the top
                        Add(sub[2],Length(bnd[2])); # so they're not all included here
                    fi;
                fi;
            fi;
        od;
    od;
    hbars:=List([1..Length(hbars)/2],x->Concatenation(hbars[2*x-1],hbars[2*x]));

    # add the vertical 1-cells
    vbars:=[];
    for i in TransposedMat(GRD){[2..Length(GRD)-1]} do
        vslice:=Filtered(i,x->x<>0);
        if vslice<>vslice*0 then
            Add(vbars,vslice);
        fi;
        for j in [1..Length(vslice)-1] do
            cell:=[2,vslice[j],vslice[j+1]];
            Add(bnd[2],cell);
            int:=List(crossings,x->Length(Intersection(cell{[2,3]},x)));
            max:=PositionMaximum(int);
            if IsBound(crs) then
                if int[max]=1 then
                    if crs[max]=1 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2]));
                    fi;
                elif int[max]=2 then
                    if crs[max]<>0 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2]));
                    fi;
                else
                    Add(sub[2],Length(bnd[2]));
                fi;
            else
                Add(sub[2],Length(bnd[2])); # vertical 1-cells default to the bottom
            fi;
        od;
    od;
    vbars:=List([1..Length(vbars)/2],x->Concatenation(vbars[2*x-1],vbars[2*x]));

    # add the loops
    for i in [2,6..Length(GRD)-4] do
        loops1:=Filtered(GRD[i],x->x<>0);
        loops2:=Filtered(GRD[i+3],x->x<>0);
        for j in [1,2] do
            Add(bnd[2],[2,loops1[1],loops2[1]]);
            Add(sub[2],Length(bnd[2])); # loops always in subcomplex
            Add(bnd[2],[2,loops1[Length(loops1)],loops2[Length(loops2)]]);
            Add(sub[2],Length(bnd[2]));
        od;
    od;

    # add the remaining four 1-cells to keep things regular
    Add(bnd[2],[2,1,2]); Add(bnd[2],[2,nr0cells-1,nr0cells]);
    Add(bnd[2],[2,1,nr0cells]); Add(bnd[2],[2,1,nr0cells]);
####################################################################################

# the 2-skeleton of the disk
####################################################################################
    unchosen:=List(bnd[2],x->x{[2,3]});
    neighbours:=List(bnd[1],x->[]);

    for i in [1..Length(bnd[1])] do
        for j in [1..Length(unchosen)] do
            if i in unchosen[j] then
                Add(neighbours[i],j);
            fi;
        od;
    od;

    Clockwise:=function(neighbours)
        local
            oriented, first0, last0,
            i, j, x, k, l, posi, posx;

        oriented:=List(neighbours,x->List([1..12],y->"pass"));
        first0:=SortedList(neighbours[1]);
        last0:=SortedList(neighbours[Length(neighbours)]);

        oriented[1][7]:=first0[1];
        oriented[1][6]:=first0[3];
        oriented[1][8]:=first0[2];
        oriented[Length(oriented)][1]:=last0[1]; 
        oriented[Length(oriented)][2]:=last0[3];
        oriented[Length(oriented)][12]:=last0[2];

        for i in [2..Length(neighbours)-1] do
            for j in [1..Length(neighbours[i])] do
                x:=bnd[2][neighbours[i][j]];
                x:=Filtered(x{[2,3]},y->y<>i)[1];
                for k in [1..Length(GRD)] do
                    for l in [1..Length(GRD[1])] do
                        if i=GRD[k][l] then
                            posi:=[k,l];
                        fi;
                        if x=GRD[k][l]then
                            posx:=[k,l];
                        fi;
                    od;
                od;
                # _\\|//_
                #  //|\\
                if posi[1]>posx[1] then
                    if posi[2]=posx[2] then
                        oriented[i][1]:=neighbours[i][j];
                    elif posi[2]<posx[2] then
                        if oriented[i][2]="pass" then
                            oriented[i][2]:=neighbours[i][j];
                        else
                            oriented[i][3]:=neighbours[i][j];
                        fi;
                    elif posi[2]>posx[2] then
                        if oriented[i][12]="pass" then
                            oriented[i][12]:=neighbours[i][j];
                        else
                            oriented[i][11]:=neighbours[i][j];
                        fi;
                    fi;
                elif posi[1]=posx[1] then
                    if posi[2]<posx[2] then
                        oriented[i][4]:=neighbours[i][j];
                    elif posi[2]>posx[2] then
                        oriented[i][10]:=neighbours[i][j];
                    fi;
                elif posi[1]<posx[1] then
                    if posi[2]=posx[2] then
                        oriented[i][7]:=neighbours[i][j];
                    elif posi[2]<posx[2] then
                        if oriented[i][5]="pass" then
                            oriented[i][5]:=neighbours[i][j];
                        else
                            oriented[i][6]:=neighbours[i][j];
                        fi;
                    elif posi[2]>posx[2] then
                        if oriented[i][9]="pass" then
                            oriented[i][9]:=neighbours[i][j];
                        else
                            oriented[i][8]:=neighbours[i][j];
                        fi;
                    fi;
                fi;
            od;
        od;
        
        return oriented;
    end;

    path:=Clockwise(neighbours);

    unselectedEdges:=List([1..Length(bnd[2])-2]);
    unselectedEdges:=Concatenation(unselectedEdges,unselectedEdges); # use Append() instead
    Add(unselectedEdges,Length(bnd[2])-1);
    Add(unselectedEdges,Length(bnd[2]));

    ClockwiseTurn:=function(p,e)
        local f;
        f:=(Position(p,e) mod 12)+1;
        while p[f]="pass" do
            f:=(f mod 12)+1;
        od;
        return p[f];
    end;

    sourceORtarget:=List([1..Length(bnd[2])],y->[3,2]);
    x:=1;
    while unselectedEdges<>[] do
        while (not x in unselectedEdges) and (not e1 in unselectedEdges) do
            x:=x+1;
        od;
        2cell:=[x];
        sORt:=sourceORtarget[x][Length(sourceORtarget[x])];
        Unbind(sourceORtarget[x][Length(sourceORtarget[x])]);

        ori:=path[bnd[2][x][sORt]];
        e0:=bnd[2][x][sORt];
        e1:=ClockwiseTurn(ori,x);
        while e1<>x do
            Add(2cell,e1);
            e0:=Filtered(bnd[2][e1]{[2,3]},y->y<>e0)[1];
            ori:=path[e0];
            e1:=ClockwiseTurn(ori,e1);
        od;
        Add(2cell,Length(2cell),1);
        if not Set(2cell) in List(bnd[3],Set) then
            for i in Filtered(2cell{[2..Length(2cell)]},y->y in unselectedEdges) do
                Unbind(unselectedEdges[Position(unselectedEdges,i)]);
            od;
            Add(bnd[3],2cell);
            if Length(Intersection(2cell{[2..2cell[1]+1]},sub[2]))=2cell[1] and
                2cell[1]<>2 then
                    Add(sub[3],Length(bnd[3]));
            fi;
        fi;
    od;
    bnd[3]:=List(bnd[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));
####################################################################################

# the direct product of the disk with [0,1]
####################################################################################
# make a duplicate of the disk
    l0:=Length(bnd[1]); l0_:=Length(sub[1]);
    l1:=Length(bnd[2]); l1_:=Length(sub[2]);
    l2:=Length(bnd[3]); l2_:=Length(sub[3]);

    IsEdgeInDuplicate:=function(k)
        local x;

        if Length(Positions(bnd[2],bnd[2][k]))=2 and
            bnd[2][k]<>[2,l0+1,2*l0] then
                return true;
        elif not Position(bnd[2],bnd[2][k]-[0,l0,l0]) in sub[2] and
            not bnd[2][k]-[0,l0,l0] in [[2,1,2],[2,1,l0],[2,l0-1,l0]] then
                return true;
        else
            x:=List(
                List(crossings,y->y+l0),
                z->Length(Intersection(z,bnd[2][k]{[2,3]}))=2
            );
            if true in x then
                if IsBound(crs) then
                    if crs[Position(x,true)] in [1,-1] then
                        return true;
                    fi;   
                else
                    return true;
                fi;
            fi;                 
        fi;
        return false;
    end;

    loops:=Filtered(bnd[2]*1,x->Length(Positions(bnd[2],x))=2);
    loops:=Set(Concatenation(List(loops,x->x{[2,3]})));
    bnd[1]:=Concatenation(bnd[1],bnd[1]);
    sub[1]:=Concatenation(sub[1],[l0+2..2*l0-1]); 
# sub contains all duplicate 0-cells except for those in the frame

    copy1:=List(bnd[2],x->x+[0,l0,l0]);
    bnd[2]:=Concatenation(bnd[2],copy1);
    for i in [l1+1..2*l1] do
        if IsEdgeInDuplicate(i) then
            Add(sub[2],i);
        fi;
    od;

    hbars2:=List(hbars,x->x+l0);
    vbars2:=List(vbars,x->x+l0);
    copy2:=List(bnd[3],x->Concatenation([x[1]],List(x{[2..Length(x)]},y->y+l1)));
    bnd[3]:=Concatenation(bnd[3],copy2);
    for i in [l2..2*l2] do
        cell:=bnd[3][i]*1;
        x:=[]; # all 0-cells in a given 2-cell
        for j in [2..Length(cell)] do
            for k in [2,3] do
                Add(x,bnd[2][cell[j]][k]);
            od;
        od;
        x:=Set(x);
        if Length(Intersection(cell{[2..cell[1]+1]},sub[2]))=cell[1] and
            (
                true in List(hbars2,y->IsSubset(y,x)) or
                true in List(vbars2,y->IsSubset(y,x))
            ) and
                cell[1]<>2 then
                    Add(sub[3],i);
        fi;
    od;
    l0__:=Length(bnd[1]); l1__:=Length(sub[2]); l2__:=Length(sub[3]);

# join the original disk to the copy
# each n-cell of the disk yields an (n+1)-cell which connects it
# to its duplicate
    for i in [1..l0] do
        Add(bnd[2],[2,i,i+l0]);
        if (i in loops) and (not i in [1,l0]) then
            Add(sub[2],Length(bnd[2]));
        fi;
    od;

    for i in [1..l1] do
        Add(
            bnd[3],
            [
                4,
                i,
                Position(bnd[2],[2,bnd[2][i][2],bnd[2][i+l1][2]]),
                Position(bnd[2],[2,bnd[2][i][3],bnd[2][i+l1][3]]),
                i+l1
            ]
        );
        if i in sub[2] and
            Length(Positions(bnd[2],bnd[2][i]))=2 then
                Add(sub[3],Length(bnd[3]));
        fi;       
    od;

    for i in [1..l2] do
        3cell:=[];
        for j in bnd[3][i]{[2..Length(bnd[3][i])]} do
            Add(
                3cell,
                Position(
                    bnd[3],
                    [
                        4,
                        j,
                        Position(bnd[2],[2,bnd[2][j][2],bnd[2][j+l1][2]]),
                        Position(bnd[2],[2,bnd[2][j][3],bnd[2][j+l1][3]]),
                        j+l1
                    ]
                )
            );
        od;
        Add(3cell,i);
        Add(3cell,i+l2);
        Add(bnd[4],Concatenation([bnd[3][i][1]+2],3cell));
    od;

    for i in [3,4] do
        bnd[i]:=List(bnd[i],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));
    od;
####################################################################################

# join the loops according to crs
####################################################################################
    colour:=List([1..4],x->[]);

    if not IsBound(crs) then
        crs:=List([1..Length(crossings)],x->1);
    fi;
    
# start with the lower caps, they're more straight forward
    lcap:=Filtered([1..l2],y->not y in sub[3] and bnd[3][y][1]<>2); # lower 'dome'
    cap:=[]; # 3-cells inside caps for each horizontal/vertical tube in lower dome
    floor:=[];
    for i in [1..l2_] do
        Add(
            floor,
            bnd[3][sub[3][i]]{[2..bnd[3][sub[3][i]][1]+1]}*1
        );
        Add(cap,[sub[3][i]]);
    od;
    leftovers:=Difference(sub[2]{[1..l1_]},Concatenation(floor));
    leftovers:=Filtered( # are there any stray loops
        leftovers, # that need to be capped off?
        x->Intersection(
            Positions(bnd[2],[2,bnd[2][x][2],bnd[2][x][3]]),
            Concatenation(floor)
        )=[]
    );
    for i in [1..Length(leftovers)-1] do
        for j in [i+1..Length(leftovers)] do
            pos:=Position( # if leftovers is non-empty
                bnd[3], # add those already existing
                [ # 2-cells (with two 1-cells in their
                    2, # boundaries) of bnd to sub
                    leftovers[i],
                    leftovers[j]
                ]
            );
            if pos<>fail then
                Add(sub[3],pos);
                Add(lcap,pos);
                break;
            fi;
        od;
    od;
    for i in [1..Length(floor)] do # floor is now a list whose length is how many
        for j in [i+1..Length(floor)] do # 2-cells we have to add to both sub & bnd
            if Intersection(floor[i],floor[j])<>[] then
                floor[i]:=Concatenation(floor[i],floor[j]);
                floor[j]:=[];
                cap[i]:=Concatenation(cap[i],cap[j]);
                cap[j]:=[];
            fi;
        od;
    od;
    floor:=Set(Filtered(List(floor,Set),x->x<>[]));
    cap:=Set(Filtered(List(cap,Set),x->x<>[]));
    for i in [1..Length(floor)] do # we need to add some more 1-cells if two
        for j in [1..Length(floor[i])] do # given bars intersect at a loop
            if floor[i][j]<>0 then
                pos:=Positions(
                    bnd[2],
                    [
                        2,
                        bnd[2][floor[i][j]][2],
                        bnd[2][floor[i][j]][3]
                    ]
                );
                if Length(pos)=2 then
                    for k in [i+1..Length(floor)] do
                        cell:=Filtered(pos,x->x<>floor[i][j])[1];
                        if cell in floor[k] then
                            Add(
                                bnd[2],
                                [
                                    2,
                                    bnd[2][floor[i][j]][2],
                                    bnd[2][floor[i][j]][3]
                                ]
                            );
                            Add(sub[2],Length(bnd[2]));
                            Add(floor[i],Length(bnd[2]));
                            Add(floor[k],Length(bnd[2]));
                            floor[i][j]:=0;
                            floor[k][Position(floor[k],cell)]:=0;
                        fi;
                    od;
                fi;
            fi;
        od;
    od;
    floor:=List(floor,x->Filtered(x,y->y<>0));
    for i in [1..Length(floor)] do # swap the loops should they occur
        for j in [1..Length(floor[i])] do # twice in bnd[2]
            pos:=Positions(bnd[2],bnd[2][floor[i][j]]);
            if Length(pos)=2 then
                floor[i][j]:=Filtered(pos,x->not x in floor[i])[1];
            fi;
        od;
    od;
    # filter out either the vertical/horizontal bars at each crossing
    # depending on if the given bar is horizontal/vertical respectively
    HorizontalOrVertical:=function(l)
        local 0_cells_l;
        0_cells_l:=Set(
            Concatenation(
                List(
                    l,
                    x->bnd[2][x]{[2,3]}
                )
            )
        );
        if true in List(hbars,x->IsSubset(x,0_cells_l)) then
            return "horizontal";
        elif true in List(vbars,x->IsSubset(x,0_cells_l)) then
            return "vertical";
        else
            return "neither";
        fi;
    end;
    # cap needs some more 2-cells, also to be sorted into fewer
    # sublists if any sub-element of cap intersects with another
    cap_:=List(cap,x->List(x,y->bnd[3][y]{[2..bnd[3][y][1]+1]}));
    cap_:=List(cap_,x->List(x,y->List(y,z->bnd[2][z])));
    cap_:=List(cap_,Concatenation);
    for i in [1..Length(cap_)] do
        for j in [1..Length(cap_[i])] do
            if Length(Positions(bnd[2],cap_[i][j]))>=2 then
                Add(
                    cap[i],
                    Filtered(
                        [1..Length(bnd[3])],
                        x->bnd[3][x][1]=2 and
                        cap_[i][j] in List(bnd[3][x]{[2,3]},y->bnd[2][y])
                    )[1]
                );
            fi;
        od;
    od;
    for i in [1..Length(cap)] do
        for j in [i+1..Length(cap)] do
            if Intersection(cap[i],cap[j])<>[] then
                cap[i]:=Concatenation(cap[i],cap[j]);
                cap[j]:=[];
            fi;
        od;
    od;
    cap:=Filtered(List(cap,Set),x->x<>[]);
    for i in [1..Length(floor)] do
        if HorizontalOrVertical(floor[i])="horizontal" then
            # remove vertical cells at crossing
            for j in [1..Length(floor[i])] do
                cell:=List(floor[i],x->bnd[2][x]{[2,3]})[j];
                for k in [1..Length(crossings)] do
                    if Length(Intersection(cell,crossings[k]))=2 and
                        cell[1]<>cell[2]-1 then
                            Unbind(floor[i][j]);
                    fi;
                od;
            od;
        elif HorizontalOrVertical(floor[i])="vertical" then
            # remove horizontal cells at crossing
            for j in [1..Length(floor[i])] do
                cell:=List(floor[i],x->bnd[2][x]{[2,3]})[j];
                for k in [1..Length(crossings)] do
                    if Length(Intersection(cell,crossings[k]))=2 and
                        cell[1]=cell[2]-1 then
                            Unbind(floor[i][j]);
                    fi;
                od;
            od;
        fi;
    od;
    floor:=List(floor,Set);
    # finally, we can create the caps as well as add some 3-cells to bnd
    SubcapTo3cell:=function(i)
        for j in [1..Length(cap)] do
            if Intersection(
                Concatenation(
                    List(
                        cap[j],
                        x->Concatenation(
                            List(
                                bnd[3][x]{[2..bnd[3][x][1]+1]},
                                y->bnd[2][y]{[2,3]}
                            )
                        )
                    )
                ),
                Concatenation(List(floor[i],x->bnd[2][x]{[2,3]}))
            )<>[] then
                return j;
            fi;
        od;
    end;
    for i in [1..Length(floor)] do
        cell:=floor[i]*1;
        Add(cell,Length(cell),1);
        Add(bnd[3],cell);
        Add(sub[3],Length(bnd[3]));
        Add(lcap,Length(bnd[3]));
        Add(cap[SubcapTo3cell(i)],Length(bnd[3]));
    od;
    for i in [1..Length(cap)] do
        cell:=cap[i]*1;
        Add(cell,Length(cell),1);
        Add(bnd[4],cell);
    od;

# now for the upper caps, 0 in crs leads to a very elaborate CW-structure
    ucap:=Filtered([l2+1..2*l2],y->not y in sub[3] and bnd[3][y][1]<>2);
    cap:=[];
    ceiling:=[];
    crs_int:=[];
    h_2_bnd:=[];
    h_3_bnd:=[];
    v_2_bnd:=[];
    v_3_bnd:=[];
    IntersectingCylinders:=function(a,b,c,d)
        local n, i, m, j, l;
# attaches to a 0 crossing some additional regular CW-structure
# to allow for a self-intersection to occur
        n:=1*Length(bnd[1])+1;
        for i in [1..8] do # 0-skeleton of intersection
            Add(bnd[1],[1,0]);
            Add(sub[1],Length(bnd[1]));
        od;
        # 1-skeleton of intersection
        m:=1*Length(bnd[2])+1;
        Add(bnd[2],[2,a,n]); Add(sub[2],Length(bnd[2])); # m
        Add(bnd[2],[2,b,n+1]); Add(sub[2],Length(bnd[2])); # m+1
        Add(bnd[2],[2,c,n+6]); Add(sub[2],Length(bnd[2])); # m+2
        Add(bnd[2],[2,d,n+7]); Add(sub[2],Length(bnd[2])); # m+3
        for i in [0..3] do
            for j in [1,2] do
                Add(bnd[2],[2,n+2*i,n+1+2*i]); # m+4, m+5, m+6, m+7, m+10, m+11, m+14, m+15
                # top first, then bottom (refer to drawing)
                Add(sub[2],Length(bnd[2]));
            od;
            if i>0 then
                Add(bnd[2],[2,n+2*i-2,n+2*i]); Add(sub[2],Length(bnd[2])); # m+8, m+12, m+16
                Add(bnd[2],[2,n+2*i-1,n+2*i+1]); Add(sub[2],Length(bnd[2])); # m+9, m+13, m+17
            fi;
        od;
        Add(h_2_bnd,[[a,b,c,d],[m,m+1,m+2,m+3,m+4,m+14]]);
        Add(v_2_bnd,[[a,b,c,d],[[m,m+1,m+4],[m+2,m+3,m+14]]]);
        # 2-skeleton of intersection
        l:=1*Length(bnd[3])+1;
        Add( # l
            bnd[3],
            [
                4,
                Position(bnd[2],[2,a,b]),
                m,
                m+1,
                m+5
            ]
        );
        Add(sub[3],Length(bnd[3]));
        Add( # l+1
            bnd[3],
            [
                4,
                Position(bnd[2],[2,c,d]),
                m+2,
                m+3,
                m+15
            ]
        );
        # these 2-cells are those which should be coloured #########################
        if IsBound(clr) then                                                      ##
            pos:=Position(List(crossings,x->Set(x)+l0),Set([a,b,c,d]));           ##
        fi;                                                                       ##
        Add(bnd[3],[4,m+4,m+6,m+8,m+9]); # l+2                                    ##
        Add(sub[3],Length(bnd[3]));                                               ##              
        Add(bnd[3],[4,m+5,m+7,m+8,m+9]); # l+3                                    ##
        Add(sub[3],Length(bnd[3]));                                               ## 
        for i in [0,1] do                                                         ##
            Add(bnd[3],[4,m+6+4*i,m+10+4*i,m+12+4*i,m+13+4*i]); # l+4, l+6        ##  
            Add(sub[3],Length(bnd[3]));                                           ##
            Add(bnd[3],[4,m+7+4*i,m+11+4*i,m+12+4*i,m+13+4*i]); # l+5, l+7        ##
            Add(sub[3],Length(bnd[3]));                                           ##
        od;                                                                       ##
        if IsBound(clr) then                                                      ##
            if clr[pos]=1 then                                                    ##
                colour[3][Length(bnd[3])-5]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-4]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            elif clr[pos]=2 then                                                  ##
                colour[3][Length(bnd[3])-5]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-4]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            elif clr[pos]=3 then                                                  ##
                colour[3][Length(bnd[3])-5]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-4]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            else                                                                  ##
                colour[3][Length(bnd[3])-5]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-4]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            fi;                                                                   ##
        fi;                                                                       ##
        ############################################################################
        # from this point onwards, cells added are only present in bnd, not sub
        Add(bnd[3],[2,m+4,m+5]); # l+8
        Add(bnd[3],[2,m+14,m+15]); # l+9
        Add(v_3_bnd,[[a,b,c,d],[[l,l+8],[l+1,l+9]]]);
        Add( # l+10
            bnd[3],
            [
                6,
                Position(bnd[2],[2,a,c]),
                m,
                m+2,
                m+8,
                m+12,
                m+16
            ]
        );
        Add( # l+11
            bnd[3],
            [
                6,
                Position(bnd[2],[2,b,d]),
                m+1,
                m+3,
                m+9,
                m+13,
                m+17
            ]
        );
        Add(h_3_bnd,[[a,b,c,d],[l+2,l+4,l+6,l+10,l+11]]);
        # 3-skeleton of intersection
        Add(
            bnd[4],
            [
                8,
                l+2,
                l+3,
                l+4,
                l+5,
                l+6,
                l+7,
                l+8,
                l+9
            ]
        );
        Add(
            bnd[4],
            [
                8,
                Position(
                    List(bnd[3],Set),
                    Set(
                        [
                            4,
                            Position(bnd[2],[2,a,b]),
                            Position(bnd[2],[2,c,d]),
                            Position(bnd[2],[2,a,c]),
                            Position(bnd[2],[2,b,d])
                        ]
                    )
                ),
                l,
                l+1,
                l+3,
                l+5,
                l+7,
                l+10,
                l+11
            ]
        );
    end;

    for i in [l2_+1..l2__] do
        Add(
            ceiling,
            bnd[3][sub[3][i]]{[2..bnd[3][sub[3][i]][1]+1]}*1
        );
        Add(cap,[sub[3][i]]);
    od;
    leftovers:=Difference(sub[2]{[l1_+1..l1__]},Concatenation(ceiling));
    leftovers:=Filtered( # are there any stray loops
        leftovers, # that need to be capped off?
        x->Intersection(
            Positions(bnd[2],[2,bnd[2][x][2],bnd[2][x][3]]),
            Concatenation(ceiling)
        )=[]
    );
    l:=[];
    for i in [1..Length(leftovers)] do
        for j in [i+1..Length(leftovers)] do
            if bnd[2][leftovers[i]]=bnd[2][leftovers[j]] then
                Add(l,[leftovers[i],leftovers[j]]);      
            fi;
        od;
    od;
    for i in [1..Length(l)] do
        pos:=Position( # if leftovers is non-empty
            bnd[3], # add those already existing
            [ # 2-cells (with two 1-cells in their
                2, # boundaries) of bnd to sub
                l[i][1],
                l[i][2]
            ]
        );
        Add(sub[3],pos);
        Add(ucap,pos);
    od;
    0_ceiling:=List(ceiling,x->List(x,y->bnd[2][y]{[2,3]}));
    0_ceiling:=List(0_ceiling,x->Set(Concatenation(x)));
    for i in [1..Length(0_ceiling)] do # for each self-intersection, add
        for j in [1..Length(crossings)] do # the additional structure
            if Length(Intersection(0_ceiling[i],crossings[j]+l0))=4 and
                crs[j]=0 then
                    IntersectingCylinders(
                        0_ceiling[i][1],
                        0_ceiling[i][2],
                        0_ceiling[i][3],
                        0_ceiling[i][4]
                    );
            fi;
        od;
    od;
    for i in [1..Length(0_ceiling)] do
        for j in [i+1..Length(0_ceiling)] do
            int:=Intersection(0_ceiling[i],0_ceiling[j]);
            if int<>[] then # if a horizontal and vertical bar intersect at a
            # loop in the corner, add a new loop
                pos:=Positions(bnd[2],[2,Minimum(int),Maximum(int)]);
                if Length(pos)=2 then
                    Add(bnd[2],[2,Minimum(int),Maximum(int)]);
                    Add(sub[2],Length(bnd[2]));
                    ceiling[i]
                    [
                        Position(
                            ceiling[i],
                            Intersection(ceiling[i],pos)[1]
                        )
                    ]:=Length(bnd[2]);
                    ceiling[j]
                    [
                        Position(
                            ceiling[j],
                            Intersection(ceiling[j],pos)[1]
                        )
                    ]:=Length(bnd[2]);
                fi;
            fi;
        od;
    od;
    for i in [1..Length(ceiling)] do
        for j in [1..Length(ceiling[i])] do
        # swap the loops at the end of a non-intersecting cap 
        # to make things join up properly
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][ceiling[i][j]][2],
                    bnd[2][ceiling[i][j]][3]
                ]
            );
            if Length(pos)=2 then
                ceiling[i][j]:=Filtered(pos,x->x<>ceiling[i][j])[1];
            fi;
        od;
    od;

    # now to finally start creating caps--we'll start with horizontal caps as they
    # do not need the extra 3-cells that the vertical caps will need at
    # self-intersection points
    hbars:=List(hbars,x->x+l0);
    vbars:=List(vbars,x->x+l0);
    x:=List(ceiling,x->HorizontalOrVertical(x)="horizontal");
    h_0_ceiling:=0_ceiling*1;
    h_ceiling:=ceiling*1;
    h_cap:=cap*1;
    for i in [1..Length(0_ceiling)] do
        if x[i]=false then
            h_0_ceiling[i]:=[];
            h_ceiling[i]:=[];
            h_cap[i]:=[];
        fi;
    od;
    for i in [1..Length(0_ceiling)] do
        for j in [i+1..Length(0_ceiling)] do
            if Intersection(h_0_ceiling[i],h_0_ceiling[j])<>[] and
                x[i]=true and
                    x[j]=true then
                        h_0_ceiling[i]:=Concatenation(
                            h_0_ceiling[i],
                            h_0_ceiling[j]
                        );
                        h_0_ceiling[j]:=[];
                        h_ceiling[i]:=Concatenation(
                            h_ceiling[i],
                            h_ceiling[j]
                        );
                        h_ceiling[j]:=[];
                        h_cap[i]:=Concatenation(
                            h_cap[i],
                            h_cap[j]
                        );
                        h_cap[j]:=[];
            fi;
        od;
    od;
    Apply(h_0_ceiling,Set);
    Apply(h_ceiling,Set);
    Apply(h_cap,Set);
    for i in [1..Length(h_cap)] do # delete any stray 2-cells that are just crossing pts
        if Length(h_cap[i])=1 then
            cell:=bnd[3][h_cap[i][1]]*1;
            cell:=cell{[2..cell[1]+1]};
            cell:=Set(Concatenation(List(cell,x->bnd[2][x]{[2,3]})));
            if true in List(crossings,x->Intersection(cell,x+l0)=cell) then
                h_0_ceiling[i]:=[];
                h_ceiling[i]:=[];
                h_cap[i]:=[];
            fi;
        fi;
    od;
    # add some degree 2 2-cells to the boundary of our new 3-cells
    for i in [1..Length(h_ceiling)] do
        for j in [1..Length(h_ceiling[i])] do
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][h_ceiling[i][j]][2],
                    bnd[2][h_ceiling[i][j]][3]
                ]
            );
            if Length(pos)>=2 then
                pos:=Set(Filtered(pos,x->x<=2*l1));
                Add(h_cap[i],Position(bnd[3],[2,pos[1],pos[2]]));
            fi;
        od;
    od;

    v_crossings:=[]; # at each crossing, there are 2 horizontal and vertical cells
    for i in [0..Length(Concatenation(crossings))/2-1] do
        Add(
            v_crossings,
            Position(
                bnd[2],
                [
                    2,
                    Concatenation(crossings)[2*i+1]+l0,
                    Concatenation(crossings)[2*(i+1)]+l0
                ]
            )
        );
    od;
    h_crossings:=[];
    for i in [0..Length(Concatenation(crossings))/2-1] do
        Add(
            h_crossings,
            Position(
                bnd[2],
                [
                    2,
                    Set(Concatenation(crossings))[2*i+1]+l0,
                    Set(Concatenation(crossings))[2*(i+1)]+l0
                ]
            )
        );
    od;

    for i in [1..Length(h_ceiling)] do # unbind the vertical 1-cells from the boundary
        if Intersection(h_ceiling[i],v_crossings)<>[] then # of the horizontal 2-cells
            for j in [1..Length(h_ceiling[i])] do
                if h_ceiling[i][j] in v_crossings then
                    Unbind(h_ceiling[i][j]);
                fi;
            od;
            h_ceiling[i]:=Set(h_ceiling[i]);
        fi;
    od;
    for i in [1..Length(h_0_ceiling)] do # adjust the boundaries of the 2-cell and
    # 3-cell which are to be added at points of self-intersection
        for j in [1..Length(h_2_bnd)] do
            if Length(Intersection(h_0_ceiling[i],h_2_bnd[j][1]))=4 then
                for k in [1..Length(h_ceiling[i])] do
                    if h_ceiling[i][k] in h_crossings then
                        Unbind(h_ceiling[i][k]);
                    fi;
                od;
                Append(h_ceiling[i],h_2_bnd[j][2]);
                h_ceiling[i]:=Set(h_ceiling[i]);
                Append(h_cap[i],h_3_bnd[j][2]);
                h_cap[i]:=Set(h_cap[i]);
            fi;
        od;
    od;
    # everything now matches up, h_ceiling's non empty entries correspond
    # to the boundaries of 2-cells which are to be added
    # note: h_cap is not yet complete, we need to do the vertical bars before that's
    # the case
    for i in [1..Length(h_ceiling)] do
        if h_ceiling[i]<>[] then
            cell:=h_ceiling[i]*1;
            cell:=Set(cell);
            Add(cell,Length(cell),1);
            Add(bnd[3],cell);
            Add(sub[3],Length(bnd[3]));
            Add(ucap,Length(bnd[3]));
            Add(h_cap[i],Length(bnd[3]));
        fi;
    od;
    # time to repeat a similar process but for the vertical bars
    HorizontalOrVertical:=function(l)
        local 0_cells_l;
        0_cells_l:=Set(
            Concatenation(
                List(
                    l,
                    x->bnd[2][x]{[2,3]}
                )
            )
        );
        if true in List(vbars,x->IsSubset(x,0_cells_l)) then
            return "vertical";
        elif true in List(hbars,x->IsSubset(x,0_cells_l)) then
            return "horizontal";
        else
            return "neither";
        fi;
    end;
    x:=List(ceiling,x->HorizontalOrVertical(x)="vertical");
    v_0_ceiling:=0_ceiling*1;
    v_ceiling:=ceiling*1;
    v_cap:=cap*1;
    #return List(v_cap,x->List(bnd[3][x[1]]{[2..bnd[3][x[1]][1]+1]},y->bnd[2][y]{[2,3]}));
    #something's messing up here
    for i in [1..Length(0_ceiling)] do
        if x[i]=false then
            v_0_ceiling[i]:=[];
            v_ceiling[i]:=[];
            v_cap[i]:=[];
        fi;
    od;
    for i in [1..Length(0_ceiling)] do
        for j in [i+1..Length(0_ceiling)] do
            if Intersection(v_0_ceiling[i],v_0_ceiling[j])<>[] and
                x[i]=true and
                    x[j]=true then
                        Append(
                            v_0_ceiling[i],
                            v_0_ceiling[j]
                        );
                        v_0_ceiling[j]:=[];
                        Append(
                            v_ceiling[i],
                            v_ceiling[j]
                        );
                        v_ceiling[j]:=[];
                        Append(
                            v_cap[i],
                            v_cap[j]
                        );
                        v_cap[j]:=[];
            fi;
        od;
    od;
    Apply(v_0_ceiling,Set);
    Apply(v_ceiling,Set);
    Apply(v_cap,Set);
    for i in [1..Length(v_cap)] do # delete any stray 2-cells that are just crossing pts
        if Length(v_cap[i])=1 then
            cell:=bnd[3][v_cap[i][1]]*1;
            cell:=cell{[2..cell[1]+1]};
            cell:=Set(Concatenation(List(cell,x->bnd[2][x]{[2,3]})));
            if true in List(crossings,x->Intersection(cell,x+l0)=cell) then
                v_0_ceiling[i]:=[];
                v_ceiling[i]:=[];
                v_cap[i]:=[];
            fi;
        fi;
    od;
    # add some degree 2 2-cells to the boundary of our new 3-cells
    for i in [1..Length(v_ceiling)] do
        for j in [1..Length(v_ceiling[i])] do
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][v_ceiling[i][j]][2],
                    bnd[2][v_ceiling[i][j]][3]
                ]
            );
            if Length(pos)>=2 then
                pos:=Set(Filtered(pos,x->x<=2*l1));
                Add(v_cap[i],Position(bnd[3],[2,pos[1],pos[2]]));
            fi;
        od;
    od;
    for i in [1..Length(v_ceiling)] do # unbind the horizontal 1-cells from the boundary
        if Intersection(v_ceiling[i],h_crossings)<>[] then # of the vertical 2-cells
            for j in [1..Length(v_ceiling[i])] do
                if v_ceiling[i][j] in h_crossings then
                    Unbind(v_ceiling[i][j]);
                fi;
            od;
            v_ceiling[i]:=Set(v_ceiling[i]);
        fi;
    od;
    # we need to order the 2-cells in each vertical bar so that they can be
    # separated into entry/exit points for each 0 crossing
    # this is done by finding the sum of the min and max labeled 0-cells in each 2-cell
    minmax:=[];
    for i in [1..Length(v_cap)] do
        if v_cap[i]<>[] then
            minmax[i]:=List(
                v_cap[i],
                x->Set(
                    Concatenation(
                        List(
                            bnd[3][x]{[2..bnd[3][x][1]+1]},
                            y->bnd[2][y]{[2,3]}
                        )
                    )
                )
            );
            minmax[i]:=List(minmax[i],x->Minimum(x)+Maximum(x));
            minmax[i]:=List(minmax[i],x->Position(Set(minmax[i]),x));
            perm:=[];
            for j in [1..Length(v_cap[i])] do
                perm[minmax[i][j]]:=v_cap[i][j]*1;
            od;
            # minmax[i] is now a permutation which orders the 2-cells
            # of v_cap[i] vertically from top to bottom
            v_cap[i]:=perm*1;
        fi;
    od;
    # now to take each vbar in v_cap and separate it into multiple components
    # if there are any self-intersections involved:
    # n self intersections => n+1 components
    crs_2_cells:=Filtered(
        [l2+1..2*l2],
        x->Set(
            Concatenation(
                List(
                    bnd[3][x]{[2..bnd[3][x][1]+1]},
                    y->bnd[2][y]{[2,3]}
                )
            )
        ) in List(crossings,x->Set(x)+l0)
    );
    for i in [1..Length(crs_2_cells)] do
        if not crs[i]=0 then
            crs_2_cells[i]:=0;
        fi;
    od; # list of 2-cells where there's a 0-crossings
    v_cap_sep:=List([1..Length(v_cap)],x->[]);
    for i in [1..Length(v_cap)] do
        if v_cap[i]<>[] and
            Intersection(v_cap[i],crs_2_cells)<>[] then
                l:=[];
                for j in [1..Length(v_cap[i])] do
                    if not v_cap[i][j] in crs_2_cells then
                        Add(l,v_cap[i][j]);
                    else
                        Add(l,v_cap[i][j]);
                        Add(v_cap_sep[i],l);
                        l:=[v_cap[i][j]*1];
                    fi;
                od;
                Add(v_cap_sep[i],l);
        elif v_cap[i]<>[] then
            v_cap_sep[i]:=v_cap[i]*1;
        fi;
    od;
    # need to make a list similar to v_cap_sep for the 1-cells in v_ceiling
    # have to be careful because the boundary of the 2-cells doesn't necessarily
    # contain all the 1-cells we need (we've added extra corner loops etc.)
    v_ceiling_sep:=List([1..Length(v_ceiling)]);
    for i in [1..Length(v_cap)] do
        if v_cap_sep[i]=[] then
            v_ceiling_sep[i]:=[];
        elif v_cap_sep[i][1]*0=0 then
            v_ceiling_sep[i]:=[];
        else
            v_ceiling_sep[i]:=List(
                [1..Length(v_cap_sep[i])],
                x->[]
            );
        fi;
    od;
    for i in [1..Length(v_ceiling_sep)] do
        if v_ceiling_sep<>[] then
            if Length(v_ceiling_sep[i])>1 then
                for j in [1..Length(v_ceiling_sep[i])] do
                    k:=Set(
                        Concatenation(
                            List(
                                Filtered(v_cap_sep[i][j],x->not x in crs_2_cells),
                                x->Concatenation(
                                    List(
                                        bnd[3][x]{[2..bnd[3][x][1]+1]},
                                        y->bnd[2][y]{[2,3]}
                                    )
                                )
                            )
                        )
                    );
                    l:=List(
                        v_ceiling[i],
                        x->bnd[2][x]{[2,3]}  
                    );
                    for m in [1..Length(l)] do
                        if Length(Intersection(k,l[m]))=2 then
                            Add(v_ceiling_sep[i][j],v_ceiling[i][m]*1);
                        fi;
                    od;
                od;
            else
                v_ceiling_sep[i]:=v_ceiling[i]*1;
            fi;
        fi;
    od;
    # v_ceiling_sep contains enough information to obtain the boundary of our capping
    # vertical 2-cells. once those are added to bnd, add their indexing to v_cap_sep's
    # 3-cells. the 3-cells can then be formed by comparing intersections (at corners only)
    # with the horizontal 3-cells in h_cap
    for i in [1..Length(v_cap_sep)] do
        if v_cap_sep[i]<>[] then
            if IsList(v_cap_sep[i][1]) then # you forgot the else statement here...
                x:=0; # binary counter for entry/exit
                for j in [1..Length(v_cap_sep[i])] do
                    if Intersection(v_cap_sep[i][j],crs_2_cells)<>[] then
                        int:=List(crossings+l0)[
                            Position(
                                crs_2_cells,
                                Intersection(
                                    v_cap_sep[i][j],
                                    crs_2_cells
                                )[1]
                            )
                        ];
                        for k in [1..Length(v_2_bnd)] do
                            if Set(int)=Set(v_2_bnd[k][1]) then
                                Append(v_ceiling_sep[i][j],v_2_bnd[k][2][x+1]);
                                cell:=v_ceiling_sep[i][j]*1;
                                Add(cell,Length(cell),1);
                                Add(bnd[3],cell);
                                Add(sub[3],Length(bnd[3]));
                                Add(v_cap_sep[i][j],Length(bnd[3]));
                                Add(ucap,Length(bnd[3]));
                                Append(v_cap_sep[i][j],v_3_bnd[k][2][x+1]);
                                Unbind(
                                    v_cap_sep[i][j][
                                        Position(
                                            v_cap_sep[i][j],
                                            Intersection(
                                                v_cap_sep[i][j],
                                                crs_2_cells
                                            )[1]
                                        )
                                    ]
                                );
                                v_cap_sep[i][j]:=Set(v_cap_sep[i][j]);
                                x:=(x+1) mod 2;
                            fi;
                        od;
                    fi;
                od;
            else
                cell:=Set(v_ceiling_sep[i])*1;
                Add(cell,Length(cell),1);
                Add(bnd[3],cell);
                Add(sub[3],Length(bnd[3]));
                Add(v_cap_sep[i],Length(bnd[3]));
                Add(ucap,Length(bnd[3]));
            fi;
        fi;
    od;

    # final step, matching the boundaries of 3-cells in h_cap with those in v_cap_sep
    3cells:=[];
    for i in [1..Length(h_cap)] do
        if h_cap[i]<>[] then
            Add(3cells,Filtered(h_cap[i],x->not x in crs_2_cells));
        fi;
    od;
    for i in [1..Length(v_cap_sep)] do
        if v_cap_sep[i]<>[] then
            if IsList(v_cap_sep[i][1]) then
                for j in [1..Length(v_cap_sep[i])] do
                    Add(3cells,v_cap_sep[i][j]);
                od;
            else
                Add(3cells,v_cap_sep[i]);
            fi;
        fi;
    od;
    for i in [1..Length(3cells)] do
        for j in [i+1..Length(3cells)] do
            if IsBound(3cells[i]) and IsBound(3cells[j]) then
                int1:=List(3cells[i],x->bnd[3][x]{[2..bnd[3][x][1]+1]});
                int2:=List(3cells[j],x->bnd[3][x]{[2..bnd[3][x][1]+1]});
                int1:=Concatenation(List(int1,x->List(x,y->bnd[2][y]{[2,3]})));
                int2:=Concatenation(List(int2,x->List(x,y->bnd[2][y]{[2,3]})));
                int1:=Filtered(int1,x->x[1]<=2*l0 and x[2]<=2*l0);
                int2:=Filtered(int2,x->x[1]<=2*l0 and x[2]<=2*l0);
                int1:=Filtered(int1,x->Length(Positions(bnd[2],[2,x[1],x[2]]))=3);
                int2:=Filtered(int2,x->Length(Positions(bnd[2],[2,x[1],x[2]]))=3);
                int1:=Set(int1); int2:=Set(int2);
                if Intersection(int1,int2)<>[] then
                    Append(3cells[i],3cells[j]);
                    Unbind(3cells[j]);
                fi;
            fi;
        od;
    od;
    3cells:=Set(3cells,Set);
    for i in [1..Length(3cells)] do
        Add(3cells[i],Length(3cells[i]),1);
        Add(bnd[4],3cells[i]*1);
    od;
###################################################################################

# add a cap to both ends of D x [0,1] 
####################################################################################
    Add(
        bnd[3],
        [
            2,
            Positions(bnd[2],[2,1,l0])[1],
            Positions(bnd[2],[2,1,l0])[2]
        ]
    );
    Add(lcap,Length(bnd[3]));
    lcap:=Set(lcap);
    Add(lcap,Length(lcap),1);
    Add(bnd[4],lcap);
    Add(
        bnd[3],
        [
            2,
            Positions(bnd[2],[2,l0+1,2*l0])[1],
            Positions(bnd[2],[2,l0+1,2*l0])[2]
        ]
    );
    Add(ucap,Length(bnd[3]));
    ucap:=Set(ucap);
    Add(ucap,Length(ucap),1);
    Add(bnd[4],ucap);
####################################################################################

# add colour
####################################################################################
    for i in [1..Length(bnd[3])] do
        if not IsBound(colour[3][i]) then
            colour[3][i]:=[0];
        fi;
    od;
    colour_:={n,k}->colour[n+1][k];
####################################################################################
    x:=CWSubcomplexToRegularCWMap([RegularCWComplex(bnd),sub]);
    x!.colour:=colour_;
    return x;
end;
#gap> b:=HopfSatohSurface();
#Pure cubical complex of dimension 4.
#
#gap> b:=RegularCWComplex(b);
#Regular CW-complex of dimension 4
#
#gap> b:=SimplifiedComplex(b);
#Regular CW-complex of dimension 4
#
#gap> Homology(b,0);
#[ 0, 0 ]
#gap> Homology(b,1);
#[ 0, 0, 0, 0 ]
#gap> Homology(b,2);
#[ 0, 0 ]
#gap> Homology(b,3);
#[  ]
#gap> Homology(b,4);
#[  ]
