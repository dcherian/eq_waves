% run sanity checks on inferred structures
function [] = RunTests(opt)

    load([opt.name '-' opt.filt.window '.mat']);

    % if corr coef is insignificant,
    % regression slope bounds should include 0.
    for mm=1:11
        for nn=1:7
            r = modes.corr{mm,nn};
            sl = modes.InferredModeOLS{mm,nn};
            err = modes.InferredModeErrorOLS{mm,nn};

            mask = abs(r) < eps;

            assert( all((sl(mask)-err(mask)) .* (sl(mask)+err(mask)) < 0), ...
                    ['Regression slope does not bound 0 for (' ...
                     num2str(mm) ', ' num2str(nn) ')']);
        end
    end
    disp(['Test passed. When corr. coeff. is insignificant, slope ' ...
          'bounds includes 0.']);

    for mm=1:11
        for nn=1:7
            if isempty(modes.intercept{mm,nn})
                continue;
            end

            c = modes.intercept{mm,nn}(:,1);
            cerr = modes.intercept{mm,nn}(:,2);

            assert(all(abs(c) < 1e-3), ...
                   ['Intercept > 1e-3 for (' ...
                     num2str(mm) ', ' num2str(nn) ')']);
            % assert(all((c+cerr) .* (c-cerr) < 0), ...
            %        ['Intercept does not bound 0 for (' ...
            %          num2str(mm) ', ' num2str(nn) ')']);

        end
    end
    disp(['Test passed. Intercept is 0.']);
