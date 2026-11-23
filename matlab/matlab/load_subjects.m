function subs = load_subjects(root)
%LOAD_SUBJECTS  Load all subjects of studies 1-3 (N = 156).
%
%  Reads LCPA_STUDY{1,2,3}_behavioral.csv from the repository root and
%  returns a cell array of per-subject structs (see prep_subject.m), each
%  with .study and .SubID added. Study 3 subject 36 is excluded (not part
%  of the analyzed sample).
%
%  root: repository root; defaults to the parent of this folder.

    if nargin < 1
        root = fileparts(fileparts(mfilename('fullpath')));
    end
    subs = {};
    for si = 1:3
        T = readtable(fullfile(root, sprintf('LCPA_STUDY%d_behavioral.csv', si)));
        if si == 3
            T(T.SubID == 36, :) = [];
        end
        ids = unique(T.SubID);
        for u = 1:numel(ids)
            sub = prep_subject(T(T.SubID == ids(u), :));
            sub.study = si;
            sub.SubID = ids(u);
            subs{end+1} = sub; %#ok<AGROW>
        end
    end
end
