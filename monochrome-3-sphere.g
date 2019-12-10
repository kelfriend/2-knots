# Input: a marked arc presentation together with a list of 0s, +1s & -1s
# detailing the type of crossing occuring at each point in the associated knot
# Output: a regular CW-decomposition of the 3-sphere towards a colouring
# which will be used to embed the space in 4-dimensional space (a 2-knot).
2knot:=function(arc,crs)
    local
        arc_presentation, crossings, grid_number, grid,
        i, IsIntersection, crossing_coords, j, crossing_type,
        bound, 2skeleton;

    arc_presentation:=ShallowCopy(arc);
    crossings:=ShallowCopy(crs);
    grid_number:=Length(arc);

    grid:=List([1..grid_number],x->[1..grid_number]*0);
    # forms a matrix from the arc presentation
    for i in [0..grid_number-1] do
        grid[grid_number-i][arc_presentation[i+1][1]]:=1;
        grid[grid_number-i][arc_presentation[i+1][2]]:=1;
    od;

    IsIntersection:=function(i,j)
    # detects the number of crossings in the grid
        if grid[i][j]=0
            then
            if 1 in grid[i]{[1..j]}
                then
                if 1 in grid[i]{[j..grid_number]}
                    then
                    if 1 in List([1..i],x->grid[x][j])
                        then
                        if 1 in List([i..grid_number],x->grid[x][j])
                            then
                            return true;
                        fi;
                    fi;
                fi;
            fi;
        fi;

        return false;
    end;

    crossing_coords:=[];
    for i in [1..grid_number] do
        for j in [1..grid_number] do
            if IsIntersection(i,j) then
                Add(crossing_coords,[i,j]);
            fi;
        od;
    od;

    if Length(crossing_coords)<>Length(crossings) then
        Error(
            "the specified number of crossings is incorrect, ",
            Length(crossing_coords),
            " crossings expected"
        );
    fi;

    crossing_type:=function(c)

        if c in crossing_coords then
            return crossings[Position(crossing_coords,c)];
        else
            return '*';
        fi;
    end;

    bound:=List([1..5],x->[]);

    2skeleton:=function(bound)
        local
            nr0cells, i, c,
            horizontal_lengths;

        nr0cells:=4*grid_number+8*(Length(crossings));
        # each crossing contributes the same number of 0-cells (eight)

        for i in [1..nr0cells] do
            Add(bound[1],[1,0]); # add the 0-cells
        od;

        c:=1;
        horizontal_lengths:=List([1..grid_number],x->[]);
        for i in [1..grid_number] do
            if i=1 then
                Add(horizontal_lengths[i],c);
                Add(bound[2],[2,c,c+2]); Add(bound[2],[2,c,c+2]); # loops
                Add(bound[2],[2,c,c+1]); # connectors
                Add(bound[2],[2,c+2,c+3]);
                Add(bound[2],[2,c+1,c+3]); Add(bound[2],[c+1,c+3]); # loops
                c:=c+4;
                Add(horizontal_lengths[i],c-1);
            else
                for j in Filtered(crossing_coords,x->x[1]=i) do
                    if crossing_type(j)<>'*' then
                        Add(horizontal_lengths,c);

                    fi;
                od;
            fi;

    end;

    2skeleton(bound);

    return bound;
end;
