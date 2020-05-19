function [Y,x,y] = ASVEDE_extract(total_time, nx, ny)
    
	addpath('./data');
    Y = [];
    x = [];
    y = [];
    for i=1:total_time
        % Setting up file name
        if i < 10
            string = "vort000" + i;
        elseif i<100
            string = "vort00" + i;
        else
            string ="vort0"+i;
        end

        table = readtable(string, 'HeaderLines', 0);
        Y(:,i) = table.Var3;
        disp('Reading vorticity from '+string+' success!');
        drawnow;
    end
    
    x = table.Var1;
    y = table.Var2;
    disp('Reading Coordinates Success!')
end
    
