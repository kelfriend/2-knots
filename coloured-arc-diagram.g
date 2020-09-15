ViewColouredArcDiagram:=function(arc,cross,levels)
    local
        AppendTo, PrintTo, file, filetxt, ArcPresentationToArray,
        bin_arr, i, res, crs, crs_coord, j, colours, n, bar1, loc_3,
        bar2, nsew, clr;

    AppendTo:=HAP_AppendTo;
    PrintTo:=HAP_PrintTo;
    file:="/tmp/HAPtmpImage";
    filetxt:="/tmp/HAPtmpImage.txt";

    ArcPresentationToArray:=function(arc)
        local
            arr, i, pair, j, k;

        arr:=List([0..5*(Length(arc)-1)],x->List([0..5*(Length(arc)-1)])*0);

        for i in [1..Length(arc)] do
            pair:=arc[i]*1;
            arr[Length(arr)-5*(i-1)][1+5*(pair[1]-1)]:=3;
            arr[Length(arr)-5*(i-1)][1+5*(pair[2]-1)]:=3;
        od;

        for i in [1..Length(arr)] do
            for j in [1..Length(arr[i])] do
                if arr[i][j]=3 then
                    k:=j*1;
                    while arr[i][k+1]<>3 do
                        k:=k+1;
                        arr[i][k]:=1;
                    od;
                    break;
                fi;
            od;
        od;
        for i in [1..Length(arr)-1] do
            for j in [1..Length(arr[i])] do
                if arr[i][j]=3 then
                    if 3 in List(arr{[i+1..Length(arr)]},x->x[j]) then
                        k:=i*1;
                        while arr[k+1][j]<>3 do
                            k:=k+1;
                            arr[k][j]:=arr[k][j]+1;
                        od;
                    fi;
                fi;
            od;
        od;

        return arr;
    end;

    bin_arr:=FrameArray(ArcPresentationToArray(arc))-3;
    bin_arr:=List(bin_arr,x->Concatenation(x,[-3,-3]));
    for i in [1,2] do
        Add(bin_arr,bin_arr[1]*0-3);
    od;

    # resolution of the output, chosen to display on (my) laptop nicely
    res:=String(Minimum(50*Length(bin_arr),950));

    crs:=List(bin_arr,x->Positions(x,-1));
    crs_coord:=[]; # the co-ordinates of each crossing
    for i in [1..Length(crs)] do
        if crs[i]<>[] then
            for j in [1..Length(crs[i])] do
                Add(crs_coord,[i,crs[i][j]]);
            od;
        fi;
    od;

    #if Length(crs_coord)<>Length(cross) then
    #    Error("number of specified crossings is incorrect");
    #elif levels<>[] then
    #    if Length(Positions(cross,0))<Length(levels) then
    #        Error("all ",Length(Positions(cross,0))," 4-d crossing(s) must be accounted for");
    #    fi;
    #fi;

    PrintTo(
        filetxt,
        "# ImageMagick pixel enumeration: ",
        Length(bin_arr[1]),
        ",",
        Length(bin_arr),
        ",255,RGB\n"
    );

    colours:=NewDictionary([0,""],true);
    colours!.entries:=[
        [-3,"(255,255,255)"], # white
        [-2,"(0,255,255)"], # cyan
        [-1,"(0,0,255)"], # blue
        [0,"(0,255,0)"], # green
        [1,"(255,0,0)"], # red
        [2,"(255,255,0)"], # yellow
        [3,"(255,0,255)"] # magenta
    ];

    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            if bin_arr[i][j] in [-1,0] then
                bin_arr[i][j]:=3; # temp. colour for the corners/crossings
            elif bin_arr[i][j]=-2 then
                bin_arr[i][j]:=0; # everything should be green 
            fi; # unless there's a 4-d crossing involved
        od;
    od;

    n:=0;
    for i in [1..Length(crs_coord)] do
        if cross[i]=0 then
            n:=n+1;
            bar1:=List([1..crs_coord[i][1]-1]);
            loc_3:=Positions(List(bar1,x->bin_arr[x][crs_coord[i][2]]),3);
            bar1:=bar1{
                [loc_3[Length(loc_3)]+1..Length(bar1)]
            };
            bar2:=List([crs_coord[i][1]+1..Length(bin_arr)]);
            loc_3:=Positions(List(bar2,x->bin_arr[x][crs_coord[i][2]]),3);
            bar2:=bar2{
                [1..loc_3[1]-1]
            };
            for j in bar1 do
                if levels[n] in [1,2] then
                    bin_arr[j][crs_coord[i][2]]:=-1;
                else
                    bin_arr[j][crs_coord[i][2]]:=1;
                fi;
            od;
            for j in bar2 do
                if levels[n] in [1,3] then
                    bin_arr[j][crs_coord[i][2]]:=-1;
                else
                    bin_arr[j][crs_coord[i][2]]:=1;
                fi;
            od;
        elif cross[i]=-1 then
            j:=crs_coord[i]*1;
            bin_arr[j[1]][j[2]-1]:=-3;
            bin_arr[j[1]][j[2]+1]:=-3;
        fi;
    od;
    for i in [1..Length(crs_coord)] do
        if cross[i]=1 then
            j:=crs_coord[i]*1;
            bin_arr[j[1]-1][j[2]]:=-3;
            bin_arr[j[1]+1][j[2]]:=-3;
        fi;
    od;

    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            if bin_arr[i][j]=3 then
                nsew:=[
                    bin_arr[i-1][j],
                    bin_arr[i+1][j],
                    bin_arr[i][j+1],
                    bin_arr[i][j-1]
                ];
                if -3 in nsew then
                    if Length(Positions([nsew[1],nsew[2]],-3))=1 then # corner
                        if -1 in nsew then
                            bin_arr[i][j]:=-2;
                        elif 1 in nsew then
                            bin_arr[i][j]:=2;
                        else
                            bin_arr[i][j]:=0;
                        fi;
                    else # int. point
                        if nsew[1]=-3 then
                            bin_arr[i][j]:=0;
                        else
                            if nsew[1]=nsew[2] then
                                bin_arr[i][j]:=nsew[1]*1;
                            elif Set([nsew[1],nsew[2]])=Set([1,0]) then
                                bin_arr[i][j]:=2;
                            elif Set([nsew[1],nsew[2]])=Set([-1,0]) then
                                bin_arr[i][j]:=-2;
                            fi;
                        fi;
                    fi;
                else
                    bin_arr[i][j]:=0;
                fi;
            fi;
        od;
    od;

    # colour entries according to levels
    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            clr:=LookupDictionary(colours,bin_arr[i][j]);
            AppendTo(filetxt,j,",",i,": ",clr,"\n");
        od;
    od;

    # convert the binary array to an upscaled png
    Exec(
        Concatenation("convert ",filetxt," ","-scale ",res,"x",res," ",file,".png")
    );
    # delete the old txt file
    Exec(
        Concatenation("rm ",filetxt)
    );
    # display the image
    Exec(
        Concatenation(DISPLAY_PATH," ","/tmp/HAPtmpImage.png")
    );
    # delete the image on window close
    Exec(
        Concatenation("rm  ","/tmp/HAPtmpImage.png")
    );
end;
