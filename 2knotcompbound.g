2KnotComplementWithBoundary:=function(arc,lvls)
    local
        arc_presentation, levels, grid_number, grid,
        i, IsIntersection, CornerConfiguration, big_grid,
        j, k, nr0cells, bound, knot_boundary,
        inc_mapping, hbars, vbars, inc, 1SkeletonOfDisk,
        2SkeletonOfDisk, cap_bound, 3SkeletonOfTube, CappedCylinder;

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

    inc:={n,k}->inc_mapping[n+1][k];

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
        for i in [1..Length(inc_mapping[1])] do
            Add(knot_boundary[1],[1,0]);
        od;

        for i in [1..2*grid_number] do
            if i<=grid_number then
                for j in Filtered(bnd[2],x->(x[2] in hbars[i] and x[3] in hbars[i])) do
                    #if levels[i]=-1 then
# if the level is + then the cells will be added later to avoid convolution
                        if not Length(Positions(knot_boundary[2],[2,j[2]-1,j[3]-1]))>2 then 
                            Add(knot_boundary[2],[2,j[2]-1,j[3]-1]);
                            if not Position(bnd[2],j) in inc_mapping[2] then
                                Add(inc_mapping[2],Position(bnd[2],j));
                            else
                                Add(inc_mapping[2],Position(bnd[2],j)+2);
                            fi;
                        fi;
                    #fi;
                od;
            else
                for j in Filtered(bnd[2],x->(x[2] in vbars[i-grid_number] and x[3] in vbars[i-grid_number])) do
                    #if levels[i]=-1 then
                        if not [2,j[2]-1,j[3]-1] in knot_boundary[2] then
                            Add(knot_boundary[2],[2,j[2]-1,j[3]-1]);
                            Add(inc_mapping[2],Position(bnd[2],j));
                        fi;
                    #fi;
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
                loops, present_loops, j, k, vertices, check;

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

            bnd[3]:=Filtered(bnd[3],x->x[1]<>2);
            bnd[3]:=List(bnd[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));

            # now to add the appropriate 2-cells to the boundary of the knot
            for i in bnd[3] do
                vertices:=[];
                for j in i{[2..5]} do
                    Add(vertices,bnd[2][j][2]);
                    Add(vertices,bnd[2][j][3]);
                od;
                vertices:=Set(vertices);
                check:=false;
                for j in [1..grid_number] do
                    if Intersection(vertices,hbars[j])=vertices then
                        check:=true;
                    fi;
                    if Intersection(vertices,vbars[j])=vertices then
                        check:=true;
                    fi;
                od;
                if check then
                    Add(
                        knot_boundary[3],
                        Concatenation(
                            [i[1]],
                            List(
                                i{[2..Length(i)]},
                                x->Position(inc_mapping[2],x)
                            )
                        )
                    );
                    Add(inc_mapping[3],Position(bnd[3],i));
                fi;
            od;
        end;
        FaceTrace(path);
    end;
################################################################################
    2SkeletonOfDisk(bound);

    cap_bound:=[[Difference(inc_mapping[3],[1..Length(bound[3])])],[]];
    cap_bound[2]:=cap_bound[1]+Length(bound[3]); 
# this list will be used in recording the boundary of the upper & lower 3-cells 
################################################################################
    3SkeletonOfTube:=function(bnd)
        local
            l0, l1, l2, k0, k1, k2,
            DuplicateDisk, JoinDisks;

        l0:=Length(bnd[1]); k0:=Length(knot_boundary[1]);
        l1:=Length(bnd[2]); k1:=Length(knot_boundary[2]);
        l2:=Length(bnd[3]); k2:=Length(knot_boundary[3]);

        DuplicateDisk:=function()
            local
                copy1, copy2, i,
                pos1, 2cell, pos2;

            bnd[1]:=Concatenation(bnd[1],bnd[1]);

            copy1:=List(bnd[2],x->[x[1],x[2]+l0,x[3]+l0]);
            bnd[2]:=Concatenation(bnd[2],copy1);

            copy2:=List(bnd[3],x->Concatenation([x[1]],List(x{[2..Length(x)]},y->y+l1)));
            bnd[3]:=Concatenation(bnd[3],copy2);

            ####################################################################
            # now to duplicate the knot_boundary and update inc_mapping
            knot_boundary[1]:=Concatenation(knot_boundary[1],knot_boundary[1]);
            inc_mapping[1]:=Concatenation(inc_mapping[1],inc_mapping[1]+l0);

            for i in [1..k1] do
                Add(knot_boundary[2],[2,knot_boundary[2][i][2]+k0,knot_boundary[2][i][3]+k0]);
                pos1:=Position(
                    bnd[2],
                    [2,inc_mapping[1][knot_boundary[2][i][2]+k0],
                    inc_mapping[1][knot_boundary[2][i][3]+k0]]
                );
                if pos1 in inc_mapping[2] then
                    Add(inc_mapping[2],pos1+2); # accounting for loops
                else
                    Add(inc_mapping[2],pos1);
                fi;
            od;

            knot_boundary[3]:=List(
                knot_boundary[3],
                x->Concatenation(
                    [x[1]],
                    Set(x{[2..Length(x)]})
                )
            );
            for i in [1..k2] do
                2cell:=knot_boundary[3][i];
                2cell:=2cell+2cell*0+k1-[k1];
                Add(knot_boundary[3],2cell);
                pos2:=Position(
                    bnd[3],
                    Concatenation(
                        [2cell[1]],
                        Set(
                            2cell{[2..Length(2cell)]},
                            x->inc_mapping[2][x]
                        )
                    )
                );
                Add(inc_mapping[3],pos2);
            od;
        end;
        DuplicateDisk();

        JoinDisks:=function()
            local
                i, j, 1cell, 2cell,
                i2cell, 3cell, i3cell;

            for i in [1..l0] do
                Add(bnd[2],[2,i,i+l0]);
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

            ####################################################################
            # \/ inclusion \/
            for i in [1..k0] do
                1cell:=[2,i,i+k0];
                Add(knot_boundary[2],1cell);
                Add(
                    inc_mapping[2],
                    Position(
                        bnd[2],
                        [2,inc_mapping[1][1cell[2]],inc_mapping[1][1cell[3]]]
                    )
                );
            od;

            for i in [1..k1] do
                2cell:=[
                    i,
                    Position(
                        knot_boundary[2],
                        [
                            2,
                            knot_boundary[2][i][2],
                            knot_boundary[2][i+k1][2]
                        ]
                    ),
                    Position(
                        knot_boundary[2],
                        [
                            2,
                            knot_boundary[2][i][3],
                            knot_boundary[2][i+k1][3]
                        ]
                    ),
                    i+k1
                ];
                2cell:=Concatenation([4],Set(2cell));
                Add(knot_boundary[3],2cell);

                i2cell:=List(2cell{[2..5]},x->inc_mapping[2][x]);
                i2cell:=Set(i2cell);
                Add(i2cell,4,1);
                Add(inc_mapping[3],Position(bnd[3],i2cell));
            od;

            for i in [1..k2] do
                3cell:=[i];
                for j in knot_boundary[3][i]{[2..Length(knot_boundary[3][i])]} do
                    Add(
                        3cell,
                        Position(
                            knot_boundary[3],
                            [
                                4,
                                j,
                                j+k1,
                                Position(
                                    knot_boundary[2],
                                    [
                                        2,
                                        knot_boundary[2][j][2],
                                        knot_boundary[2][j+k1][2]
                                    ]
                                ),
                                Position(
                                    knot_boundary[2],
                                    [
                                        2,
                                        knot_boundary[2][j][3],
                                        knot_boundary[2][j+k1][3]
                                    ]
                                ),
                            ]
                        )
                    );
                od;
                Add(3cell,i+k2);
                3cell:=Set(3cell);
                i3cell:=List(3cell,x->inc_mapping[3][x]);
                Add(3cell,Length(3cell),1); 
                i3cell:=Set(i3cell);
                Add(i3cell,Length(i3cell),1);

                Add(knot_boundary[4],3cell);
                Add(inc_mapping[4],Position(bnd[4],i3cell));
            od;

        end;
        JoinDisks();
    end;
################################################################################
    3SkeletonOfTube(bound);
################################################################################
    CappedCylinder:=function(bnd)
        local
            loops, lo_loops, hi_loops,
            caps, i, h, v, j, left, right,
            loop, 2cell;

        loops:=Filtered(bnd[2],x->Length(Positions(bnd[2],x))=2);
        lo_loops:=Set(loops{[1..(Length(loops)/2)-2]});
        hi_loops:=Set(loops{[(Length(loops)/2)+1..Length(loops)-2]});

        caps:=List([1..grid_number],x->[[],[]]);

        for i in [1..2*grid_number] do # plug the degree 2 faces
            h:=0;
            v:=0;
            for j in [1..grid_number] do
                if lo_loops[i][2] in hbars[j] and lo_loops[i][3] in hbars[j] then
                    h:=h+j;
                fi;
                if lo_loops[i][2] in vbars[j] and lo_loops[i][3] in vbars[j] then
                    v:=v+j;
                fi;
            od;
            if levels[h]=levels[grid_number+v] then
                if levels[h]=1 then
                    Add(caps[Int((i+1)/2)][((i-1) mod 2)+1],-1);

                    Add(
                        bnd[3],
                        [
                            2,
                            Positions(bnd[2],lo_loops[i])[1],
                            Positions(bnd[2],lo_loops[i])[2]
                        ]
                    );
                    Add(inc_mapping[3],Length(bnd[3]));
                    Add(
                        knot_boundary[3],
                        [
                            2,
                            Positions(knot_boundary[2],lo_loops[i]-[0,1,1])[1],
                            Positions(knot_boundary[2],lo_loops[i]-[0,1,1])[2]
                        ]
                    );
                elif levels[h]=-1 then
                    Add(caps[Int((i+1)/2)][((i-1) mod 2)+1],1);

                    Add(
                        bnd[3],
                        [
                            2,
                            Positions(bnd[2],hi_loops[i])[1],
                            Positions(bnd[2],hi_loops[i])[2]
                        ]
                    );
                    Add(inc_mapping[3],Length(bnd[3]));
                    Add(
                        knot_boundary[3],
                        [
                            2,
                            Positions(knot_boundary[2],hi_loops[i]-[0,3,3])[1],
                            Positions(knot_boundary[2],hi_loops[i]-[0,3,3])[2]
                        ]
                    );
                fi;
            fi;
        od;

    # tube time
        for i in [1..grid_number] do
            if levels[i]=-1 then
                if Length(caps[i][1])=0 then
                    left:=false;
                else
                    left:=true;
                    loop:=lo_loops[2*i-1];
                    Add(bnd[2],[2,loop[2],loop[3]]);
                    Add(inc_mapping[2],Length(inc_mapping[2]));
                    Add(knot_boundary[2],[2,loop[2]-1,loop[3]-1]);
                fi;
                if Length(caps[i][2])=0 then
                    right:=false;
                else
                    right:=true;
                    loop:=lo_loops[2*i];
                    Add(bnd[2],[2,loop[2],loop[3]]);
                    Add(inc_mapping[2],Length(inc_mapping[2]));
                    Add(knot_boundary[2],[2,loop[2]-1,loop[3]-1]);
                fi;
                if not left and not right then
                elif not left and right then
                elif left and not right then
                elif left and right then
                    2cell:=[];
                    tube:=List([1..Length(hbars[i])-4],x->x+Length(bnd[1]));
                    outer:=List(tube,x->x+Length(tube));
                    for j in [1..Length(tube)*2] do
                        Add(bnd[1],[1,0]);
                        Add(inc_mapping[1],Length(bnd[1]));
                        Add(knot_boundary[1],[1,0]);
                    od;

                    for j in [1..(Length(hbars[i])/2)-1] do
                        Add(2cell,Position(bnd[2],[2,hbars[i][j],hbars[i][j+1]]));
                    od;
                    for j in [(Length(hbars[i])/2)+1..Length(hbars[i])-1] do
                        Add(2cell,Position(bnd[2],[2,hbars[i][j],hbars[i][j+1]]));
                    od;
                    Add(2cell,Positions(bnd[2],lo_loops[2*i-1])[1]);
                    Add(2cell,Positions(bnd[2],lo_loops[2*i])[1]);
                    2cell:=Set(2cell);
                    Add(2cell,Length(2cell),1);
                    Add(bnd[3],2cell);
                fi;
            else

            fi;
        od;

        return caps;
    end;
################################################################################
    CappedCylinder(bound);

    return [bound,knot_boundary,inc_mapping];
end;
i:=[
    [ [ 2, 4 ], [ 1, 3 ], [ 2, 4 ], [ 1, 3 ] ],
    0*[1..8]+1
];;
