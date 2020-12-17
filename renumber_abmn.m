function [ele, conf] = renumber_abmn(survey_all)

    conf_full = survey_all.data(:, 1:4);
    conf = conf_full;
    [ele, ~, map2full] = unique(survey_all.ele, 'rows');
    for ii = 1:length(map2full)
       if map2full == ii
           % skip
       else
           tmp = conf_full == ii;
           conf(tmp) = map2full(ii);
       end
    end

end
