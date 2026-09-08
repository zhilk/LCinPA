function behavioralData = addConditioningMeans(behavioralData, opts)
    
%  Append per-subject low/high conditioning-mean columns.
%  Adds AvHighCond and AvLowCond to each experiment table. Values are the
%  subject's mean VASRating during conditioning, mapped onto ALL their rows.
%  If either column already exists in a table, that table is skipped.
 
    arguments
        behavioralData cell
        opts.ExcludeNoResponse (1,1) logical = true
        opts.ExcludeCatch      (1,1) logical = true
    end

    newCols = {'AvLowCond','AvHighCond'};

    for e = 1:numel(behavioralData)
        T = behavioralData{e};

        exists = ismember(newCols, T.Properties.VariableNames);
        if any(exists)
            warning('Exp %d: %s already present — skipping.', ...
                    e, strjoin(newCols(exists), ', '));
            continue
        end

        % ---- per-subject conditioning means ----
        keep = strcmp(T.Phase,'conditioning');
        if opts.ExcludeNoResponse && ismember('VASResponse', T.Properties.VariableNames)
            keep = keep & T.VASResponse == 1;
        end
        if opts.ExcludeCatch && ismember('CatchTrial', T.Properties.VariableNames)
            keep = keep & T.CatchTrial == 0;
        end
        cond = T(keep, :);

        g = groupsummary(cond, {'SubID','CueAssociation'}, 'mean', 'VASRating');
        w = unstack(g, 'mean_VASRating', 'CueAssociation', 'GroupingVariables', 'SubID');

        highMean = local_getcol(w, 'high');
        lowMean  = local_getcol(w, 'low');

        % ---- map subject means onto every row ----
        [tf, loc] = ismember(T.SubID, w.SubID);
        T.AvHighCond = nan(height(T),1);
        T.AvLowCond  = nan(height(T),1);
        T.AvHighCond(tf) = highMean(loc(tf));
        T.AvLowCond(tf)  = lowMean(loc(tf));

        behavioralData{e} = T;
    end
end

function col = local_getcol(w, name)
    if ismember(name, w.Properties.VariableNames)
        col = w.(name);
    else
        col = nan(height(w),1);   % category absent this experiment
    end
end
