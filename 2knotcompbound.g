2KnotComplementWithBoundary:=function(arc,lvls)
    local
        arc_presentation, levels, grid_number, grid,
        i, IsIntersection, CornerConfiguration, big_grid,
        j, k, nr0cells, bound, knot_boundary,
        inc_mapping, 1SkeletonOfDisk, 2SkeletonOfDisk;

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
# detects the number of crossings in the grid
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
# assigns an orientation to the various corners in the arc presentation
        if grid[i][j]=1 then
            if Size(Positions(grid[i]{[j..grid_number]},1))=2 then
                if Size(Positions(List([i..grid_number],x->grid[x][j]),1))=2
                    then # Corner type 1, i.e : __
                    return 1; #                |
                elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                    then # Corner type 3, i.e :
                    return 3; #                |__
                fi;
            elif Size(Positions(grid[i]{[1..j]},1))=2 then
                if Size(Positions(List([i..grid_number],x->grid[x][j]),1))=2
                    then # Corner type 2, i.e : __
                    return 2; #                   |
                elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                    then # Corner type 4, i.e :
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
    inc_mapping:=ShallowCopy(bound)*1; # these will be used to create an inclusion
    # from the boundary of the knot to this complement space

    1SkeletonOfDisk:=function(bnd)
# adds to bound a regular CW-decomposition of a solid disk
        local
            i, nr, hslice, j,
            vslice, index, loops1, loops2;

        for i in [1..nr0cells] do
            Add(bnd[1],[1,0]);
        od;

        # add the horizontal 1-cells
        nr:=0; # this just counts which horizontal bar the loop is currently in
        for i in [2..Length(big_grid)-1] do
            hslice:=Filtered(big_grid[i],x->x<>0);
            if hslice<>[] then
                nr:=nr+1;
            fi;
            for j in [1..Length(hslice)-1] do
                Add(bnd[2],[2,hslice[j],hslice[j+1]]);
                if levels[Int((nr-1)/2)+1]=-1 then
                # if this bar is given -ve sign in the
                # input, then it also appars in the knot boundary
                    Add(knot_boundary[2],[2,hslice[j]-1,hslice[j+1]-1]);
                    Add(inc_mapping[2],Length(bnd[2]));
                fi;
            od;
        od;

        # add the vertical 1-cells (almost identically so)
        nr:=0;
        for i in TransposedMat(big_grid){[2..Length(big_grid)-1]} do
            vslice:=Filtered(i,x->x<>0);
            if vslice<>[] then
                nr:=nr+1;
            fi;
            for j in [1..Length(vslice)-1] do
                Add(bnd[2],[2,vslice[j],vslice[j+1]]);
                if levels[Int((nr-1)/2)+1]=-1 then
                    Add(knot_boundary[2],[2,vslice[j]-1,vslice[j+1]-1]);
                    Add(inc_mapping[2],Length(bnd[2]));
                fi;
            od;
        od;

        # add the loops (my brother)
        index:=1;
        for i in [2,6..Length(big_grid)-4] do
            loops1:=Filtered(big_grid[i],x->x<>0);
            loops2:=Filtered(big_grid[i+3],x->x<>0);
            for j in [1,2] do
                Add(bnd[2],[2,loops1[1],loops2[1]]);
                Add(bnd[2],[2,loops1[Length(loops1)],loops2[Length(loops2)]]);
                # these 'loops' also appear in the boundary of the knot
                Add(knot_boundary[2],[2,index,index+2]);
                Add(knot_boundary[2],[2,index+1,index+3]);
                Add(inc_mapping[2],Length(bnd[2])-1);
                Add(inc_mapping[2],Length(bnd[2]));
            od;
            index:=index+4;
        od;

        # add the remaining 4 1-cells to keep things regular
        Add(bound[2],[2,1,2]); Add(bound[2],[2,nr0cells-1,nr0cells]);
        Add(bound[2],[2,1,nr0cells]); Add(bound[2],[2,1,nr0cells]);

        # lastly, for the 1-skeleton, add the 0-cells to knot_boundary
        for i in [1..nr0cells-2] do
            Add(knot_boundary[1],[1,0]);
        od;

        return bnd;
    end;

    1SkeletonOfDisk(bound);

    2SkeletonOfDisk:=function(bnd)
        local
            orientation, neighbours,
            NeighbourMatch, i, j;

# each 0-cell will have the 1-cells connecting it to its neighbours oriented.
# there are 12 ways of doing this. it can be visualised as the time on a clock.
        orientation:=List([1..Length(bnd[1])],x->[1..12]*0);

        # the first and last 0-cells are exceptions, their orientations are
        # assigned separately
        orientation[1][3]:=Length(bnd[1]);
        orientation[1][4]:=2;
        orientation[1][6]:=Length(bnd[1]);
        
        orientation[Length(bnd[1])][9]:=1;
        orientation[Length(bnd[1])][10]:=Length(bnd[1])-1;
        orientation[Length(bnd[1])][12]:=1;

        # simply a list of each 0-cell's neighbours
        neighbours:=List([1..Length(bnd[1])],x->[]);
        for i in List(bnd[2],x->[x[2],x[3]]) do
            Add(neighbours[i[1]],i[2]);
            Add(neighbours[i[2]],i[1]);
        od;
        neighbours:=List(neighbours,x->Set(x));

################################################################################
        NeighbourMatch:=function(i,j)
            local
                x, pos_i, pos_j;
            
            for x in [1..Length(big_grid)] do
                if Position(big_grid[x],i)<>fail then
                    pos_i:=[Position(big_grid[x],i),x];
                fi;
                if Position(big_grid[x],j)<>fail then
                    pos_j:=[Position(big_grid[x],j),x];
                fi;
            od;
            
            if pos_i[2]=pos_j[2] then
                if pos_i[1]<pos_j[1] then
                    orientation[i][3]:=j;
                else
                    orientation[i][9]:=j;
                fi;
            elif pos_i[2]<pos_j[2] then
                if pos_i[1]=pos_j[1] then
                    orientation[i][6]:=j;
                elif pos_i[1]<pos_j[1] then
                    orientation[i][4]:=j;
                    orientation[i][5]:=j;
                else
                    orientation[i][7]:=j;
                    orientation[i][8]:=j;
                fi;
            else
                if pos_i[1]=pos_j[1] then
                    orientation[i][12]:=j;
                elif pos_i[1]<pos_j[1] then
                    orientation[i][1]:=j;
                    orientation[i][2]:=j;
                else
                    orientation[i][10]:=j;
                    orientation[i][11]:=j;
                fi;
            fi;
        end;
################################################################################

        for i in [2..Length(neighbours)-1] do
            for j in neighbours[i] do
                NeighbourMatch(i,j);
            od;
        od;
# orientations are done, just implement facetrace and finish the 2-skeleton of the disk
        return orientation;
    end;

    return 2SkeletonOfDisk(bound);
end;
i:=[
    [ [ 2, 4 ], [ 1, 3 ], [ 2, 4 ], [ 1, 3 ] ],
    0*[1..8]+1
];;
