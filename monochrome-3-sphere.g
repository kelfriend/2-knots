# Input: a marked arc presentation together with a list of 0s, +1s & -1s
# detailing the type of crossing occuring at each point in the associated knot
# Output: a regular CW-decomposition of the 3-sphere towards a colouring
# which will be used to embed the space in 4-dimensional space (a 2-knot).
2knot:=function(marked_arc)
    local
        arc_presentation, crossings, grid_number, grid,
        i, IsIntersection, crossing_coords, j, crossing_type;

    arc_presentation:=ShallowCopy(marked_arc[1]);
    crossings:=ShallowCopy(marked_arc[2]);
    grid_number:=Length(marked_arc[1]);

    grid:=List([1..grid_number],x->[1..grid_number]*0);
    for i in [0..grid_number-1] do
        grid[grid_number-i][arc_presentation[i+1][1]]:=1;
        grid[grid_number-i][arc_presentation[i+1][2]]:=1;
    od;

    IsIntersection:=function(i,j)

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
        Error("the specified number of crossings is incorrect");
    fi;

    crossing_type:={c}->crossings[Position(crossing_coords,c)];

    return "Stuff works";
end;