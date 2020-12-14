function simple_taylor(m, dm, h)
    % Simplified Haber test
    %
    % clc();
    % close('all');
    % clearvars;

    err_f    = 0e-2;    % Fehler auf gesamte Lsg
    scale_f  = 1 + 0e-2;% const. Faktor   in Lsg
    err_df   = 0e-6;    % Fehler auf gesamte Ableitung der Lsg
    scale_df = 1 + 1e-2;% const. Faktor   in Ableitung der Lsg

    % Referenzlösung für homog. HR
    f = @(x, m, err, scale) err+(scale)/(x*m*2*pi*7.5);
    df = @(x, m, err, scale) err+(-scale)/(x*m^2*2*pi*7.5);

    x = -10:0.5:10;
    x = x(:);

    u = f(x, m, err_f, scale_f);
    J = df(x, m, err_df, scale_df);

    [e0, e1] = deal(zeros(length(h), 1));

    for k = 1:length(h)
        upert = f(x, m + h(k) * dm, err_f , scale_f);
    %     plot(x, f(x, m, 0), x, f(x, m + h(k) * dm, 0), ...
    %         x, f(x,m) + h(k) * J * dm, 0);

        e0(k) = norm(upert - u);
        e1(k) = norm(upert - u - h(k) * J * dm);

    end

    %%
    figure(101)
    loglog( ...
        h, e0 / e0(1), 'x-m', ...
        h, e1 / e1(1), 'r-o', ...
        h, h / h(1), 'k--', ...
        h, (h ./ h(1)).^2, 'k.')
    title(sprintf('Taylor simplified (fak = %d)', err_df));
    legend('e_0(h)', 'e_1(h)', 'O(h)', 'O(h^2)', 'Location', 'SouthWest');
    xlabel('h');
    ylabel('error norm');
    ylim([1e-17, 1e2]);
    xlim([min(h), max(h)]);
    set(gca, 'XDir', 'reverse');
end
