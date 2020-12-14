function cl = plotSens(mesh, J_plot, f_num, ABMN_ele, ylim_min, ylim_max, cl)

    y_pl_var = -0.3;
    y_pl_sym = -1.3;

    Plot.plotMesh(mesh, J_plot, f_num);
    set(gca, 'Ydir', 'normal');
    xlim([min(ABMN_ele)-5, max(ABMN_ele)+5]);
    ylim([ylim_min, ylim_max]);
    colormap(Plot.bluewhitered)
    colorbar
    if nargin < 7
        cl = caxis;
    else
        caxis(cl);
    end
    hold on
        plot([ABMN_ele(1), ABMN_ele(2)]', zeros(2,1), 'v');
        if length(ABMN_ele)==4
                    plot([ABMN_ele(3), ABMN_ele(4)]', zeros(2,1), 'o');
                    plot(ABMN_ele(1), 0, 'v');
                    text(ABMN_ele(1),y_pl_var,'A');
                    text(ABMN_ele(1),y_pl_sym, num2str(ABMN_ele(1)));
                    plot(ABMN_ele(2), 0, 'v');
                    text(ABMN_ele(2),y_pl_var,'B');
                    text(ABMN_ele(2),y_pl_sym, num2str(ABMN_ele(2)));
                    plot(ABMN_ele(3), 0, 'o');
                    text(ABMN_ele(3),y_pl_var,'M');
                    text(ABMN_ele(3),y_pl_sym, num2str(ABMN_ele(3)));
                    plot(ABMN_ele(4), 0, 'o');
                    text(ABMN_ele(4),y_pl_var,'N');
                    text(ABMN_ele(4),y_pl_sym, num2str(ABMN_ele(4)));
        elseif length(ABMN_ele)==3
                    plot(ABMN_ele(1), 0, 'vy');
                    text(ABMN_ele(1),y_pl_var,'A');
                    text(ABMN_ele(1),y_pl_sym,num2str(ABMN_ele(1)));
                    plot(ABMN_ele(2), 0, 'vg');
                    text(ABMN_ele(2),y_pl_var,'M');
                    text(ABMN_ele(2),y_pl_sym,num2str(ABMN_ele(2)));
                    plot(ABMN_ele(3), 0, 'vg');
                    text(ABMN_ele(3),y_pl_var,'N');
                    text(ABMN_ele(3),y_pl_sym,num2str(ABMN_ele(3)));
        elseif length(ABMN_ele)==2
                    plot(ABMN_ele(1), 0, 'vy');
                    text(ABMN_ele(1),y_pl_var,'A');
                    text(ABMN_ele(1),y_pl_sym,num2str(ABMN_ele(1)));
                    plot(ABMN_ele(2), 0, 'vg');
                    text(ABMN_ele(2),y_pl_var,'M');
                    text(ABMN_ele(2),y_pl_sym,num2str(ABMN_ele(2)));
        else
                    error('Something went wrong with electrode plotting')
        end
    hold off
end
