function ph = phenotype_simple(sub, sus)
%PHENOTYPE_SIMPLE  Assignment phenotype of the simple model at sigma_us.
%
%  Per block: run the forward pass and ask whether the majority of morph
%  trials is assigned to a cause founded during conditioning ('learned')
%  or to a new cause ('novel'). The subject is "novel" only if every block
%  is novel-dominant, otherwise "categorizer".

    doms = strings(0);
    for bi = 1:numel(sub.ub)
        m  = sub.blocks == sub.ub(bi);
        ix = find(m);
        r  = lcm_infer_simple(sub.X(ix, :), sus);
        ci = sub.isCond(ix);
        mi = find(sub.isMix(ix));
        cc = [];                     % causes founded during conditioning
        for k = 1:3
            fk = find(r.z == k, 1);
            if ~isempty(fk) && ci(fk)
                cc(end+1) = k; %#ok<AGROW>
            end
        end
        if ~isempty(mi)
            if isempty(cc)
                fNov = 1.0;
            else
                fNov = mean(~ismember(r.z(mi), cc));
            end
            if fNov >= 0.5
                doms(end+1) = "novel";   %#ok<AGROW>
            else
                doms(end+1) = "learned"; %#ok<AGROW>
            end
        end
    end
    if all(doms == "novel")
        ph = "novel";
    else
        ph = "categorizer";
    end
end
