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
        check, l0, l1, l2, l1_, l2_, IsEdgeInDuplicate,
        copy1, hbars2, vbars2, copy2, 3cell, colour, lcap, reg,
        closure, ucap, l1__, l2__, path_comp, pipes, HorizontalOrVertical,
        l, AboveBelow0Cell, IntersectingCylinders, pos, colour_;

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

# (i) the 0-skeleton of the disk
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

# (ii) the 1-skeleton of the disk
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
            Add(sub[2],Length(bnd[2])); # loops always in subcomplex... for now
            Add(bnd[2],[2,loops1[Length(loops1)],loops2[Length(loops2)]]);
            Add(sub[2],Length(bnd[2]));
        od;
    od;

    # add the remaining four 1-cells to keep things regular
    Add(bnd[2],[2,1,2]); Add(bnd[2],[2,nr0cells-1,nr0cells]);
    Add(bnd[2],[2,1,nr0cells]); Add(bnd[2],[2,1,nr0cells]);
####################################################################################

# (iii) the 2-skeleton of the disk
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

# (iv) the direct product of the disk with [0,1]
####################################################################################
# make a duplicate of the disk
    l0:=Length(bnd[1]);
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
# sub contains all duplicate 0-cells except for those in the frame... for now

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
    l1__:=Length(sub[2]); l2__:=Length(sub[3]);

# 13-02-21 remove any path components of sub that consist of two 0-cells and two 1-cells only
    path_comp:=PathComponentsCWSubcomplex([RegularCWComplex(bnd),sub]);
    for i in [1..Length(path_comp)] do
        if Length(path_comp[i][2][1])=2 and
            Length(path_comp[i][2][2])=2 then
                for j in [1,2] do
                    sub[j]:=Difference(sub[j],path_comp[i][2][j]);
                od;
        fi;
    od;

# join the original disk to the copy
# each n-cell of the disk yields an (n+1)-cell which connects it
# to its duplicate
    for i in [1..l0] do
        Add(bnd[2],[2,i,i+l0]);
        if i in loops and
            i in sub[1] and
                i+l0 in sub[1] and
                    not i in [1,l0] then
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
            Length(Positions(bnd[2],bnd[2][i]))=2 and
                i+l1 in sub[2] then
                Add(sub[3],Length(bnd[3]));
        fi;       
    od;

    for i in [1..l2] do
        3cell:=[];
        # 13-02-21 plug the corner holes
            if bnd[3][i][1]=2 and
                bnd[3][i][2] in sub[2] and
                    not bnd[3][i][2]+l1 in sub[2] then
                        Add(sub[3],i);
            fi;
            if bnd[3][i][1]=2 and
                not bnd[3][i][2] in sub[2] and
                    bnd[3][i][2]+l1 in sub[2] then
                        Add(sub[3],i+l2);
            fi;
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

# (v) join the open ended 'pipes' at either end of B^2xI according to crs
####################################################################################
    colour:=List([1..4],x->[]);

    if not IsBound(crs) then
        crs:=List([1..Length(crossings)],x->1);
    fi;

# (a) start with the lower caps, they're more straight forward #####################
    lcap:=[];
    reg:=RegularCWComplex(bnd);
    for i in [1..l2] do
        if bnd[3][i][1]<>2 and not i in sub[3] then
            Add(lcap,i);
        elif bnd[3][i][1]=2 and not i in sub[3] then
            closure:=ClosureCWCell(reg,2,i);
            if not IsSubset(sub[2],closure[2][2]) then
                Add(lcap,i);
            fi;
        fi;
    od;

    pipes:=[]; # the 0,1 and 2-cells of each 'pipe' in the lower dome which will join
    # together to form the (intersecting) tubular surface
    for i in [1..l2] do
        if i in sub[3] then
            Add(
                pipes,
                [
                    [i],
                    bnd[3][i]{[2..bnd[3][i][1]+1]},
                    Set(
                        Concatenation(
                            List(
                                bnd[3][i]{[2..bnd[3][i][1]+1]},
                                x->bnd[2][x]{[2,3]}
                            )
                        )
                    )
                ]
            );
        fi;
    od;
    pipes:=Filtered(pipes,x->Length(x[3])>2);
    for i in [1..Length(pipes)] do
        for j in [i+1..Length(pipes)] do
            if Intersection(
                Intersection(
                    pipes[i][3],pipes[j][3]
                ),
                Concatenation(crossings)
            )<>[] then
                Append(pipes[i][1],pipes[j][1]);
                Append(pipes[i][2],pipes[j][2]);
                Append(pipes[i][3],pipes[j][3]);
                pipes[j][3]:=[];
            fi;
        od;
    od;
    Apply(pipes,x->[Set(x[1]),Set(x[2]), Set(x[3])]);
    pipes:=Filtered(pipes,x->x[3]<>[]);
    for i in [1..Length(pipes)] do # either end of each pipe has two 2-cells which
    # are currently absent from pipe[i][1], this step will add them
        for j in Filtered(bnd[3],x->x[1]=2) do
            int:=Intersection(pipes[i][2],j{[2,3]});
            if int<>[] then
                Add(pipes[i][1],Position(bnd[3],j));
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do # if two pipes' 0-skeletons intersect, add a
        for j in [i+1..Length(pipes)] do # new 1-cell whose boundary is their intersection
            int:=Intersection(pipes[i][3],pipes[j][3]);
            if int<>[] then
                Add(int,2,1);
                Add(bnd[2],int);
                Add(sub[2],Length(bnd[2]));
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do # replace the 1-cells that occur twice with
        for j in [1..Length(pipes[i][2])] do # either the new cell we just added or the
        # other occurance of that cell which is not in pipes[i][2]
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][pipes[i][2][j]][2],
                    bnd[2][pipes[i][2][j]][3]
                ]
            );
            if Length(pos)>=2 then
                cell:=Last(Filtered(pos,x->x<>pipes[i][2][j]));
                pipes[i][2][j]:=cell;
            fi;
        od;
    od;
    HorizontalOrVertical:=function(l)
        if l in hbars then
            return "horizontal";
        fi;
        return "vertical";
    end;
    for i in [1..Length(pipes)] do # if a pipe contains 1-cells from crossings then
        # filter out the horizontal/vertical 1-cells
        # should said pipe be vertical/horizontal itself
        if Intersection(pipes[i][3],Concatenation(crossings))<>[] then
            ori:=HorizontalOrVertical(pipes[i][3]);
            l:=Filtered(crossings,x->Intersection(x,pipes[i][3])<>[])[1];
            if ori="horizontal" then
                pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,l[1],l[2]]),
                        Position(bnd[2],[2,l[3],l[4]])
                    ]
                );
            else
                pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,l[1],l[3]]),
                        Position(bnd[2],[2,l[2],l[4]])
                    ]
                );
            fi;
        fi;
    od;
    for i in [1..Length(pipes)] do # add new 2-cells to bnd[3] whose boundary is
    # the pipes[i][2] 1-cells
        cell:=Set(pipes[i][2]);
        Add(cell,Length(cell),1);
        Add(bnd[3],cell);
        Add(sub[3],Length(bnd[3]));
        Add(lcap,Length(bnd[3]));
        Add(pipes[i][1],Length(bnd[3]));
    od;
    for i in [1..Length(pipes)] do # find where pipes intersect in their 2-skeleta
    # use this to form the 3-cells comprising the interiors of the pipes
        for j in [i+1..Length(pipes)] do
            if Intersection(pipes[i][1],pipes[j][1])<>[] then
                Append(pipes[i][1],pipes[j][1]);
                pipes[j][1]:=[];
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do
        if pipes[i][1]<>[] then
            cell:=Set(pipes[i][1]);
            Add(cell,Length(cell),1);
            Add(bnd[4],cell);
        fi;
    od;

# (b) now for the upper caps, 0 in crs leads to a very elaborate CW-structure #######
# 13-02-21 this function was added to help alter the intersection structure below
# after discovering it was incorrect
    AboveBelow0Cell:=function(n,str)
        local n_, i, j, l;

        if n>l0 then
            n_:=n-l0;
        else
            n_:=n*1;
        fi;

        i:=Position(List(GRD,x->n_ in x),true);
        j:=Position(GRD[i],n_);
        
        if str="above" then
            l:=Reversed([1..i-1]);
        else
            l:=[i+1..Length(GRD)];
        fi;

        for k in l do
            if GRD[k][j]<>0 then
                if n>l0 then
                    return GRD[k][j]+l0;
                else
                    return GRD[k][j];
                fi;
            fi;
        od;
    end;
    
    ucap:=[];
    reg:=RegularCWComplex(bnd);
    for i in [l2+1..2*l2] do
        if bnd[3][i][1]<>2 and not i in sub[3] then
            Add(ucap,i);
        elif bnd[3][i][1]=2 and not i in sub[3] then
            closure:=ClosureCWCell(reg,2,i);
            if not IsSubset(sub[2],closure[2][2]) then
                Add(ucap,i);
            fi;
        fi;
    od;

    pipes:=[]; # the 0,1 and 2-cells of each 'pipe' in the upper dome which will join
    # together to form the (intersecting) tubular surface
    for i in [l2+1..2*l2] do
        if i in sub[3] then
            Add(
                pipes,
                [
                    [i],
                    bnd[3][i]{[2..bnd[3][i][1]+1]},
                    Set(
                        Concatenation(
                            List(
                                bnd[3][i]{[2..bnd[3][i][1]+1]},
                                x->bnd[2][x]{[2,3]}
                            )
                        )
                    )
                ]
            );
        fi;
    od;
    pipes:=Filtered(pipes,x->Length(x[3])>2);
    
    IntersectingCylinders:=function(a,b,c,d) ######################################
        local
            n, i, m, j, l, a_abv,
            b_abv, c_blw, d_blw,
            del_2_cell_1, del_2_cell_2,
            del_2_cell_1_i, del_2_cell_2_i,
            dif, f_clr, f_clr1, f_clr2;
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

        # 13-02-21 need to remove cells from sub to accommodate the
        # new self-intersection structure
        a_abv:=AboveBelow0Cell(a,"above");
        b_abv:=AboveBelow0Cell(b,"above");
        c_blw:=AboveBelow0Cell(c,"below");
        d_blw:=AboveBelow0Cell(d,"below");
        del_2_cell_1:=[
            Position(bnd[2],[2,a_abv,a]),
            Position(bnd[2],[2,b_abv,b]),
            Position(bnd[2],[2,a,b])
        ];
        del_2_cell_1_i:=Position(
            List(
                bnd[3],
                x->Length(
                    Intersection(
                        x{[2..x[1]+1]},
                        del_2_cell_1
                    )
                )=3
            ),
            true
        );
        del_2_cell_2:=[
            Position(bnd[2],[2,c,c_blw]),
            Position(bnd[2],[2,d,d_blw]),
            Position(bnd[2],[2,c,d])
        ];
        del_2_cell_2_i:=Position(
            List(
                bnd[3],
                x->Length(
                    Intersection(
                        x{[2..x[1]+1]},
                        del_2_cell_2
                    )
                )=3
            ),
            true
        );
        sub[3]:=Difference(sub[3],[del_2_cell_1_i,del_2_cell_2_i]);
        sub[2]:=Difference(
            sub[2],
            [
                Position(bnd[2],[2,a_abv,a]),
                Position(bnd[2],[2,b_abv,b]),
                Position(bnd[2],[2,c,c_blw]),
                Position(bnd[2],[2,d,d_blw])
            ]
        );

        # 14-02-21 this new structure has 4 additional 1-cells
        Add(bnd[2],[2,a_abv,n]); Add(sub[2],Length(bnd[2])); # m+18
        Add(bnd[2],[2,b_abv,n+1]); Add(sub[2],Length(bnd[2])); # m+19
        Add(bnd[2],[2,n+6,c_blw]); Add(sub[2],Length(bnd[2])); # m+20
        Add(bnd[2],[2,n+7,d_blw]); Add(sub[2],Length(bnd[2])); # m+21

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
        Add(sub[3],Length(bnd[3]));
        # these 2-cells are those which should be coloured #########################
        if IsBound(clr) then                                                      ##
            pos:=Position(List(crossings,x->Set(x)+l0),Set([a,b,c,d]));           ##
            pos:=pos-Length(Filtered(crs{[1..pos-1]},x->x<>0));                   ##
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
        # 14-02-21 there are four more 2-cells                                    ##
        dif:=Difference(bnd[3][del_2_cell_1_i]{[2..5]},del_2_cell_1)[1];          ##
        Add( # l+8 bottom                                                         ##
            bnd[3],                                                               ##
            [                                                                     ##
                4,                                                                ##
                dif,                                                              ##
                m+5,                                                              ##
                m+18,                                                             ##
                m+19                                                              ##
            ]                                                                     ##
        );                                                                        ##
        Add(sub[3],Length(bnd[3]));                                               ##
        dif:=Difference(bnd[3][del_2_cell_2_i]{[2..5]},del_2_cell_2)[1];          ##
        Add( # l+9 bottom                                                         ##
            bnd[3],                                                               ##
            [                                                                     ##
                4,                                                                ##
                dif,                                                              ##
                m+15,                                                             ##
                m+20,                                                             ##
                m+21                                                              ##
            ]                                                                     ##
        );                                                                        ##
        Add(sub[3],Length(bnd[3]));                                               ##                                             ##
        if IsBound(clr) then                                                      ##
            if clr[pos]=1 then                                                    ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            elif clr[pos]=2 then                                                  ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            elif clr[pos]=3 then                                                  ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            else                                                                  ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            fi;                                                                   ##
        fi;                                                                       ##
        ############################################################################
        Add(bnd[3],[2,m+4,m+5]); # l+10
        Add(sub[3],Length(bnd[3]));
        Add(bnd[3],[2,m+14,m+15]); # l+11
        Add(sub[3],Length(bnd[3]));
        # from this point onwards, cells added are only present in bnd, not sub
        Add( # l+12
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
        Add( # l+13
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
        # 14-02-21 another few 2-cells to bnd only
        Add( # l+14
            bnd[3],
            [
                3,
                Position(bnd[2],[2,a_abv,a]),
                m,
                m+18
            ]
        );
        Add( # l+15
            bnd[3],
            [
                3,
                Position(bnd[2],[2,b_abv,b]),
                m+1,
                m+19
            ]
        );
        Add( # l+16
            bnd[3],
            [
                3,
                Position(bnd[2],[2,c,c_blw]),
                m+2,
                m+20
            ]
        );
        Add( # l+17
            bnd[3],
            [
                3,
                Position(bnd[2],[2,d,d_blw]),
                m+3,
                m+21
            ]
        );
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
                l+10,
                l+11
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
                l+12,
                l+13
            ]
        );
        # 14-02-21 need two more 3-cells to finish off the structure
        Add(
            bnd[4],
            [
                5,
                del_2_cell_1_i,
                l,
                l+8,
                l+14,
                l+15
            ]
        );
        Add(
            bnd[4],
            [
                5,
                del_2_cell_2_i,
                l+1,
                l+9,
                l+16,
                l+17
            ]
        );

        if IsBound(clr) then
            f_clr:=clr[pos]*1;
            if f_clr=1 then
                f_clr1:=-1;
                f_clr2:=-1;
            elif f_clr=2 then
                f_clr1:=-1;
                f_clr2:=1;
            elif f_clr=3 then
                f_clr1:=1;
                f_clr2:=-1;
            else
                f_clr1:=1;
                f_clr2:=1;
            fi;
        else
            f_clr1:=0;
            f_clr2:=0;
        fi;

        for i in [1..Length(pipes)] do
            if IsBound(pipes[i][1][1]) then
                if pipes[i][1][1]=del_2_cell_1_i then
                    pipes[i]:=[
                        [
                            Position(
                                bnd[3],
                                Concatenation(
                                    [2],
                                    Set(
                                        Positions(
                                            bnd[2],
                                            Concatenation(
                                                [2],
                                                Set(
                                                    [a_abv,b_abv]
                                                )
                                            )
                                        )
                                    )
                                )
                            ),
                            l+8,
                            l+10
                        ],
                        [
                            Filtered(
                                pipes[i][2],
                                x->Length(
                                    Positions(
                                        bnd[2],
                                        bnd[2][x]
                                    )
                                )=2
                            )[1],
                            m+4,
                            m+18,
                            m+19
                        ],
                        [
                            a_abv,
                            b_abv,
                            n,
                            n+1
                        ],
                        '*', # just to mark this pipe as the vertical part of an intersection
                        f_clr1
                    ];
                elif pipes[i][1][1]=del_2_cell_2_i then
                    pipes[i]:=[
                        [
                            Position(
                                bnd[3],
                                Concatenation(
                                    [2],
                                    Set(
                                        Positions(
                                            bnd[2],
                                            Concatenation(
                                                [2],
                                                Set(
                                                    [c_blw,d_blw]
                                                )
                                            )
                                        )
                                    )
                                )
                            ),
                            l+9,
                            l+11
                        ],
                        [
                            Filtered(
                                pipes[i][2],
                                x->Length(
                                    Positions(
                                        bnd[2],
                                        bnd[2][x]
                                    )
                                )=2
                            )[1],
                            m+14,
                            m+20,
                            m+21
                        ],
                        [
                            c_blw,
                            d_blw,
                            n+6,
                            n+7
                        ],
                        '*',
                        f_clr2
                    ];
                fi;
            fi;
        od;
        for i in [1..Length(pipes)] do
            if not IsBound(pipes[i][4]) then
                int:=Set(
                    Intersection(
                        pipes[i][3],
                        [a,b,c,d]
                    )
                );
                if Length(int)=4 then
                    if crs[Position(List(crossings+l0,Set),int)]=0 then
                        pipes[i]:=[
                            [
                                l+2,
                                l+4,
                                l+6,
                                l+12,
                                l+13
                            ],
                            [
                                m,
                                m+1,
                                m+2,
                                m+3,
                                m+4,
                                m+14
                            ],
                            [
                                a,
                                b,
                                c,
                                d,
                                n,
                                n+1,
                                n+6,
                                n+7
                            ],
                        ];
                    fi;
                fi;
            fi;
        od;
    end; ##########################################################################

    for i in [1..Length(crs)] do
        if crs[i]=0 then
            IntersectingCylinders(
                crossings[i][1]+l0,
                crossings[i][3]+l0,
                crossings[i][2]+l0,
                crossings[i][4]+l0
            );
        fi;
    od;
    for i in [1..Length(pipes)] do # ignoring all self-intersections, join all
    # pipes which intersect at a common crossing point
        for j in [i+1..Length(pipes)] do
            if not IsBound(pipes[i][4]) and not IsBound(pipes[j][4]) then
                int:=Intersection(pipes[i][3],pipes[j][3]);
                if Length(int)>=2 and Intersection(int,Concatenation(crossings)+l0)=int then
                    Append(pipes[i][1],pipes[j][1]);
                    pipes[i][1]:=Set(pipes[i][1]); pipes[j][1]:=[];
                    Append(pipes[i][2],pipes[j][2]);
                    pipes[i][2]:=Set(pipes[i][2]); pipes[j][2]:=[];
                    Append(pipes[i][3],pipes[j][3]);
                    pipes[i][3]:=Set(pipes[i][3]); pipes[j][3]:=[];
                fi;
            fi;
        od;
    od;
    pipes:=Filtered(pipes,x->x<>[[],[],[]]);
    for i in [1..Length(pipes)] do # add the degree two 2-cells to the ends of each pipe
    # except for the vertical pipes at each intersection (they already have them)
        if not IsBound(pipes[i][4]) then
            for j in Filtered(bnd[3]{[l2+1..2*l2]},x->x[1]=2) do
                int:=Intersection(pipes[i][2],j{[2,3]});
                if int<>[] then
                    Add(pipes[i][1],Position(bnd[3],j));
                fi;
            od;
        fi;
    od;
    for i in [1..Length(pipes)] do # if two pipes' 0-skeletons intersect (at somewhere other
    # than a crossing), add a new 1-cell whose boundary is their intersection
        for j in [i+1..Length(pipes)] do
            int:=Intersection(pipes[i][3],pipes[j][3]);
            if int<>[] then
                if Intersection(int,Concatenation(crossings)+l0)=[] and
                    int<[2*l0,2*l0] then
                        Add(int,2,1);
                        Add(bnd[2],int);
                        Add(sub[2],Length(bnd[2]));
                fi;
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do # replace the 1-cells (added pre IntersectingCylinders)
    # that occur twice with either the new cell we just added or the other occurance of
        for j in [1..Length(pipes[i][2])] do # that cell which is not in pipes[i][2]
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][pipes[i][2][j]][2],
                    bnd[2][pipes[i][2][j]][3]
                ]
            );
            if Length(pos)>=2 and
                [
                    bnd[2][pipes[i][2][j]][2],
                    bnd[2][pipes[i][2][j]][3]
                ]<[2*l0,2*l0] then
                    cell:=Last(Filtered(pos,x->x<>pipes[i][2][j]));
                    pipes[i][2][j]:=cell;
            fi;
        od;
    od;
    HorizontalOrVertical:=function(l)
        if l in hbars+l0 then
            return "horizontal";
        fi;
        return "vertical";
    end;
    for i in [1..Length(pipes)] do # if a pipe contains 1-cells from crossings then
        # filter out the horizontal/vertical 1-cells
        # should said pipe be vertical/horizontal itself
        # however if a pipe contains a self intersection remove both the horizontal
        # AND vertical 1-cells
        int:=Intersection(pipes[i][3],Concatenation(crossings)+l0);
        if Length(int)=4 then
            pos:=Position(List(crossings,Set)+l0,Set(int));
            if crs[pos]=0 then
                pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][2]+l0]),
                        Position(bnd[2],[2,crossings[pos][3]+l0,crossings[pos][4]+l0]),
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][3]+l0]),
                        Position(bnd[2],[2,crossings[pos][2]+l0,crossings[pos][4]+l0])
                    ]
                );
            else
                ori:=HorizontalOrVertical(pipes[i][3]);
                if ori="horizontal" then
                    pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][2]+l0]),
                        Position(bnd[2],[2,crossings[pos][3]+l0,crossings[pos][4]+l0]),
                    ]
                );
                else
                    pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][3]+l0]),
                        Position(bnd[2],[2,crossings[pos][2]+l0,crossings[pos][4]+l0])
                    ]
                );
                fi;
            fi;
        fi;
    od;
    for i in [1..Length(pipes)] do # add new 2-cells to bnd[3] whose boundary is
    # the pipes[i][2] 1-cells
        cell:=Set(pipes[i][2]);
        Add(cell,Length(cell),1);
        Add(bnd[3],cell);
        Add(sub[3],Length(bnd[3]));
        Add(ucap,Length(bnd[3]));
        Add(pipes[i][1],Length(bnd[3]));
        if IsBound(pipes[i][4]) then
            if IsBound(clr) then
                colour[Length(bnd[3])]:=pipes[i][5];
            fi;
        fi;
    od;
    for i in [1..Length(pipes)] do # find where pipes intersect in their 2-skeleta
    # use this to form the 3-cells comprising the interiors of the pipes
        for j in [i+1..Length(pipes)] do
            if Intersection(pipes[i][1],pipes[j][1])<>[] then
                Append(pipes[i][1],pipes[j][1]);
                pipes[j][1]:=[];
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do
        if pipes[i][1]<>[] then
            cell:=Set(pipes[i][1]);
            Add(cell,Length(cell),1);
            Add(bnd[4],cell);
        fi;
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
    #return List(Filtered(sub[2],x->RegularCWComplex(bnd)!.coboundaries[2][x][1]=5),y->bnd[2][y]);
    return x;
end;
