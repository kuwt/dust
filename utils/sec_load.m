function sec_data = sec_load(file_out, type)
    switch type
        case 'l'
            field = {'vel_outplane_isolated';
                'vel_outplane';
                'vel_2d_isolated';
                'vel_2d';
                'Mo';
                'Fz';
                'Fy';
                'Fx';
                'Cm';
                'Cl';
                'Cd';
                'alpha_isolated';
                'alpha'};
        case 'v'
            field = {'vel_outplane_isolated';
                'vel_outplane';
                'vel_2d_isolated';
                'vel_2d';
                'Cm';
                'Cl';
                'Cd';
                'alpha_isolated';
                'alpha'};
    end
    
    % name file
    sec_data = struct('n_time',[],'n_sec',[], 'sec',[], 'l_sec',[],'vel_outplane',[],...
        'vel_outplane_isolated',[],'vel_2d_isolated',[],'vel_2d',[],'Mo',[],'Fx',[],...
        'Fy',[],'Fz',[],'Cm',[],'Cl',[],'Cd',[],'alpha_isolated',[],'alpha',[]);
    
    for i = 1:numel(field)
        file = [file_out '_' field{i} '.dat'];
        sec_data = parse_file_sectional(file, field{i}, sec_data);
    end

end


function sec_data = parse_file_sectional(file, field, sec_data)

    fid = fopen(file,'rt');
    it = 0;
    skip_line = false;
    sec_data.n_time = 0; % initialization
    long_field = {  'Mo';
        'Fz';
        'Fy';
        'Fx'};
    
    while skip_line || ~feof(fid)
    
        if ~skip_line
            line_new = fgets(fid);
        else
            skip_line = false;
        end
    
        it = it + 1;
    
    
        if it == 2
            line_2 = textscan(line_new, '# n_sec : %d ; n_time : %d. Next lines: y_cen , y_span');
            sec_data.n_sec = line_2{1};
            sec_data.n_time = line_2{2};
        end
    
        if it == 3
            line_3 = textscan(line_new,[repmat('%n',[1, sec_data.n_sec])]);
            sec_data.sec = cell2mat(line_3(1:end));
        end
    
        if it == 4
            line_4 = textscan(line_new,[repmat('%n',[1, sec_data.n_sec])]);
            sec_data.l_sec = cell2mat(line_4(1:end));
        end
    
        if it > 5
            while skip_line || ~feof(fid)
    
                if ~skip_line
                    line = fgets(fid);
                else
                    skip_line = false;
                end
    
                if any(strcmpi(field,long_field))
                    line = textscan(line, repmat('%n',[1, 1 + sec_data.n_sec + 9 + 3]));
                    sec_data.(field).time(it-5,1) = line{1};
                    sec_data.(field).value(it-5,:) = cell2mat(line(2:end-12));
                    sec_data.(field).R(:,:,it-5) = reshape(cell2mat(line(end-11:end-3)),3,3);
                    sec_data.(field).X(:,:,it-5) = line(end-2:end);
                    it = it + 1;
                else
                    line = textscan(line, repmat('%n',[1, sec_data.n_sec + 1]));
                    sec_data.(field).time(it-5,1) = line{1};
                    sec_data.(field).value(it-5,:) = cell2mat(line(2:end));
                    it = it + 1;
                end
    
            end
        elseif it > 5 + sec_data.n_time
            break
        end
    end
end




