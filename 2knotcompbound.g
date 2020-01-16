2KnotComplementWithBoundary:=function(arc,lvls)
    local
        arc_presentation, levels, grid_number, grid,
        i, IsIntersection, CornerConfiguration, big_grid,
        j, k, nr0cells, bound, knot_boundary,
        inc_mapping, hbars, vbars, 1SkeletonOfDisk,
        2SkeletonOfDisk;

    arc_presentation:=ShallowCopy(arc);
    levels:=ShallowCopy(lvls);
    grid_number:=Length(arc);

    if Length(levels)<>2*grid_number then
        Error(
            "the specified list of levels is incorrect, ",
            2*grid_number,
            " entries expected"
        );
    fi;

    grid:=List([1..grid_number],x->[1..grid_number]*0);
    # forms a matrix from the arc presentation
    for i in [0..grid_number-1] do
        grid[grid_number-i][arc_presentation[i+1][1]]:=1;
        grid[grid_number-i][arc_presentation[i+1][2]]:=1;
    od;

    IsIntersection:=function(i,j)
        if grid[i][j]=0 then
            if 1 in grid[i]{[1..j]} then
                if 1 in grid[i]{[j..grid_number]} then
                    if 1 in List([1..i],x->grid[x][j]) then
                        if 1 in List([i..grid_number],x->grid[x][j]) then
                            return true;
                        fi;
                    fi;
                fi;
            fi;
        fi;

        return false;
    end;

    CornerConfiguration:=function(i,j);
        if grid[i][j]=1 then
            if Size(Positions(grid[i]{[j..grid_number]},1))=2 then
                if Size(Positions(List([i..grid_number],x->grid[x][j]),1))=2 then
                         # Corner type 1, i.e : __
                    return 1; #                |
                elif Size(Positions(List([1..i],x->grid[x][j]),1))=2 then
                         # Corner type 3, i.e :
                    return 3; #                |__
                fi;
            elif Size(Positions(grid[i]{[1..j]},1))=2 then
                if Size(Positions(List([i..grid_number],x->grid[x][j]),1))=2 then
                         # Corner type 2, i.e : __
                    return 2; #                   |
                elif Size(Positions(List([1..i],x->grid[x][j]),1))=2 then
                         # Corner type 4, i.e :
                    return 4; #                 __|
                fi;
            fi;
        fi;

        return 0;
    end;
    
    big_grid:=List([1..4*grid_number],x->[1..4*grid_number]*0);
# quadruple the size of the grid to allow for the 0-skeleton to be displayed
# nicely in an array
    for i in [1..grid_number] do
        for j in [1..grid_number] do
            if CornerConfiguration(i,j) in [1,4] then
                big_grid[4*i-3][4*j-3]:=1;
                big_grid[4*i][4*j]:=1;
            elif CornerConfiguration(i,j) in [2,3] then
                big_grid[4*i-3][4*j]:=1;
                big_grid[4*i][4*j-3]:=1;
            elif IsIntersection(i,j) then
                for k in [0..3] do
                    big_grid[4*i-3][4*j-3+k]:=1;
                    big_grid[4*i][4*j-3+k]:=1;
                od;
            fi;
        od;
    od;
# labels the 0-cells by a horizontal ordering
    nr0cells:=2;
    for i in [1..4*grid_number] do
        for j in [1..4*grid_number] do
            if big_grid[i][j]=1 then
                big_grid[i][j]:=nr0cells;
                nr0cells:=nr0cells+1;
            fi;
        od;
    od;
    big_grid:=FrameArray(big_grid);
    big_grid[1][1]:=1; big_grid[4*grid_number+2][4*grid_number+2]:=nr0cells;

    bound:=List([1..6],x->[]); # the beginnings of the boundary list

    knot_boundary:=ShallowCopy(bound)*1;
    inc_mapping:=ShallowCopy(bound)*1; # these will be used to create an
    # inclusion from the boundary of the knot to this complement space
    hbars:=List([1..grid_number],x->[]);
    vbars:=List([1..grid_number],x->[]);

################################################################################
    1SkeletonOfDisk:=function(bnd)
    # adds to bound a regular CW-decomposition of a solid disk
        local
            i, hslice, j, vslice, loops1, loops2, k;

        for i in [1..nr0cells] do
            Add(bnd[1],[1,0]);
        od;

        # add the horizontal 1-cells
        for i in [2..Length(big_grid)-1] do
            hslice:=Filtered(big_grid[i],x->x<>0);
            for j in [1..Length(hslice)-1] do
                Add(bnd[2],[2,hslice[j],hslice[j+1]]);
            od;
        od;

        # add the vertical 1-cells (almost identically so)
        for i in TransposedMat(big_grid){[2..Length(big_grid)-1]} do
            vslice:=Filtered(i,x->x<>0);
            for j in [1..Length(vslice)-1] do
                Add(bnd[2],[2,vslice[j],vslice[j+1]]);
            od;
        od;

        # add the loops
        for i in [2,6..Length(big_grid)-4] do
            loops1:=Filtered(big_grid[i],x->x<>0);
            loops2:=Filtered(big_grid[i+3],x->x<>0);
            for j in [1,2] do
                Add(bnd[2],[2,loops1[1],loops2[1]]);
                Add(bnd[2],[2,loops1[Length(loops1)],loops2[Length(loops2)]]);
            od;
        od;

        # add the remaining 4 1-cells to keep things regular
        Add(bnd[2],[2,1,2]); Add(bnd[2],[2,nr0cells-1,nr0cells]);
        Add(bnd[2],[2,1,nr0cells]); Add(bnd[2],[2,1,nr0cells]);

        ########################################################################
        # now to add the appropriate cells to knot_boundary to work towards an
        # inclusion
        for i in [1..grid_number] do
            for j in big_grid{[2..5]+4*(i-1)} do
                for k in Filtered(j,x->x<>0) do
                    Add(hbars[i],k);
                od;
            od;
        od;
        for i in [2..Length(big_grid)-1] do
            for j in [1..grid_number] do
                for k in big_grid[i]{[2..5]+4*(j-1)} do
                    if k<>0 then
                        Add(vbars[j],k);
                    fi;
                od;
            od;
        od;

        inc_mapping[1]:=[2..Length(bnd[1])-1];
        inc_mapping[1]:=Concatenation(inc_mapping[1],inc_mapping[1]+Length(bnd[1])-1);
        for i in [1..Length(inc_mapping[1])] do
            Add(knot_boundary[1],[1,0]);
        od;

        for i in [1..2*grid_number] do
            if i<=grid_number then
                for j in Filtered(bnd[2],x->(x[2] in hbars[i] and x[3] in hbars[i])) do
                    if levels[i]=-1 then
# if the level is + then the cells will be added later to avoid convolution
                        if not Length(Positions(knot_boundary[2],[2,j[2]-1,j[3]-1]))>2 then 
                            Add(knot_boundary[2],[2,j[2]-1,j[3]-1]);
                            if not Position(bnd[2],j) in inc_mapping[2] then
                                Add(inc_mapping[2],Position(bnd[2],j));
                            else
                                Add(inc_mapping[2],Position(bnd[2],j)+2);
                            fi;
                        fi;
                    fi;
                od;
            else
                for j in Filtered(bnd[2],x->(x[2] in vbars[i-grid_number] and x[3] in vbars[i-grid_number])) do
                    if levels[i]=-1 then
                        if not [2,j[2]-1,j[3]-1] in knot_boundary[2] then
                            Add(knot_boundary[2],[2,j[2]-1,j[3]-1]);
                            Add(inc_mapping[2],Position(bnd[2],j));
                        fi;
                    fi;
                od;
            fi;
        od;
        ########################################################################
        return bnd;
    end;
################################################################################
    1SkeletonOfDisk(bound);
################################################################################
# this is mostly self-plagiarised from knotcomp
    2SkeletonOfDisk:=function(bnd)
        local
            Orient, Clockwise, path, FaceTrace;

        Orient:=function()
            local 
                unchosen, neighbours, i, j,
                Clockwise;

            unchosen:=List(ShallowCopy(bnd[2]),x->[x[2],x[3]]);
            neighbours:=List(ShallowCopy(bnd[1]),x->[]);

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
                        x:=bound[2][neighbours[i][j]];
                        x:=Filtered(x{[2,3]},y->y<>i)[1];
                        for k in [1..Length(big_grid)] do
                            for l in [1..Length(big_grid[1])] do
                                if i=big_grid[k][l] then
                                    posi:=[k,l];
                                fi;
                                if x=big_grid[k][l]then
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

            return Clockwise(neighbours);
        end;

        path:=Orient();

        FaceTrace:=function(path)
            local
                unselectedEdges, sourceORtarget, faceloops,
                x, ClockwiseTurn, 2cell, sORt, ori, e1, e0, i,
                loops, present_loops, j;

            unselectedEdges:=List([1..Length(bnd[2])-2]);
            unselectedEdges:=Concatenation(unselectedEdges,unselectedEdges);
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
                if (not Set(2cell) in List(bnd[3],x->Set(x))) then
                    for i in Filtered(2cell{[2..Length(2cell)]},y->y in unselectedEdges) do
                        Unbind(unselectedEdges[Position(unselectedEdges,i)]);
                    od;
                    Add(bnd[3],2cell);
                fi;
            od;

            # if loops are present in the bottom layer,
            # then they're not filtered out below
            loops:=Filtered(bnd[3],x->x[1]=2);
            present_loops:=[];
            for i in [1..grid_number] do
                if levels[i]=-1 then
                    for j in Filtered(loops,x->bnd[2][x[2]][2] in hbars[i] and bnd[2][x[2]][3] in hbars[i]) do
                        Add(present_loops,j);
                    od;
                fi;
                if levels[grid_number+i]=-1 then
                    for j in Filtered(loops,x->bnd[2][x[2]][2] in vbars[i] and bnd[2][x[2]][3] in vbars[i]) do
                        Add(present_loops,j);
                    od;
                fi;
            od;
            bnd[3]:=Filtered(bnd[3],x->not x in Difference(loops,present_loops));
            bnd[3]:=List(bnd[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));

        end;
        FaceTrace(path);
    end;
################################################################################
    2SkeletonOfDisk(bound);
    return bound;#[bound,knot_boundary,inc_mapping];
end;
i:=[
    [ [ 2, 4 ], [ 1, 3 ], [ 2, 4 ], [ 1, 3 ] ],
    0*[1..8]+1
];;
