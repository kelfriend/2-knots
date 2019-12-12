2KnotComplement:=function(arc,lvls)
    local
        arc_presentation, levels, grid_number, grid,
        i, IsIntersection, CornerConfiguration, big_grid,
        j, k, nr0cells, bound, 1SkeletonOfDisk;

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

    1SkeletonOfDisk:=function(bnd)
# adds to bound a regular CW-decomposition of a solid disk
        local
            i;

        for i in [1..nr0cells] do
            Add(bnd[1],[1,0]);
        od;

        return bnd;
    end;

    1SkeletonOfDisk(bound);

    return bound;
end;
i:=[
    [ [ 2, 4 ], [ 1, 3 ], [ 2, 4 ], [ 1, 3 ] ],
    0*[1..8]+1
];;
