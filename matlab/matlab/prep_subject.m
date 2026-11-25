function sub = prep_subject(subT)
%PREP_SUBJECT  Build the per-subject struct used by all fitters and analyses.
%
%  Sorts one subject's trials by TrialGlobal and extracts:
%    US      delivered heat, TargetVAS/100 (stimulus features live in [0,1])
%    CS      cue morph coordinates [x_face, x_house]
%    X       stimulus feature matrix [US, x_face] (latent-cause models;   
%            x_house = 1 - x_face is redundant and would double-count cues)
%    CR      raw VAS ratings (0-100)
%    blocks  block label per trial; ub: block order of appearance
%    cat     visual category per trial ('face' / 'house' / 'mix')
%  and the trial masks:
%    valid   rating given on the trial
%    isCond  conditioning phase
%    isMix   morph (mixed-cue) test trials
%    vt      valid test trials - the scoring set of all primary fits

    subT = sortrows(subT, 'TrialGlobal');
    sub.T      = subT;
    sub.n      = height(subT);
    sub.US     = subT.TargetVAS / 100;
    sub.CS     = [subT.x_face, subT.x_house];
    sub.X      = [sub.US, sub.CS(:, 1)];
    sub.CR     = subT.VASRating;
    sub.blocks = string(subT.Block);
    sub.ub     = unique(sub.blocks, 'stable');
    sub.cat    = string(subT.VisualCategory);
    sub.valid  = subT.VASResponse == 1;
    sub.isCond = strcmp(subT.Phase, 'conditioning');
    sub.isMix  = strcmp(subT.Phase, 'test') & sub.cat == "mix";
    sub.vt     = sub.valid & ~sub.isCond;
end
